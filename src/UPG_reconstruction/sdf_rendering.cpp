#include "sdf_rendering.h"
#include "graphics_utils/image_metrics.h"
#include "tinyEngine/engine.h"
#include "common_utils/distribution.h"
#include "preprocessing.h"
#include "LiteMath_conv.h"
#include <chrono>
#include <omp.h>
#include <thread>

namespace upg
{
  inline glm::vec3 EyeRayDirNormalized(float x/*in [0,1]*/, float y/*in [0,1]*/, glm::mat4 projInv)
  {
    glm::vec4 pos = glm::vec4(2.0f*x - 1.0f, -2.0f*y + 1.0f, 0.0f, 1.0f);
    pos = projInv * pos;
    pos /= pos.w;
    return glm::normalize(glm::vec3(pos));
  }

  inline glm::vec3 transformRay(glm::vec3 ray, glm::mat4 viewInv)
  {
    glm::vec3 p1 = glm::vec3(viewInv*glm::vec4(0,0,0,1));
    glm::vec3 p2 = glm::vec3(viewInv*glm::vec4(ray.x,ray.y,ray.z,1));
    return glm::normalize(p2-p1);
  }

  //intersects SDF with ray start_pos + t*dir. Returns if ray intersected with SDF and optionally - point of intersection
  bool sdf_sphere_tracing(const ProceduralSdf &sdf, const AABB &sdf_bbox, const glm::vec3 &start_pos, const glm::vec3 &dir, glm::vec3 *surface_pos = nullptr)
  {
    constexpr float EPS = 3*1e-6;
    glm::vec3 p0 = start_pos;
    if (!sdf_bbox.contains(p0))
    {
      float t = 0;
      if (sdf_bbox.intersects(start_pos, dir, &t))
        p0 = start_pos + t*dir;
      else //ray won't intersect SDF
        return false;
    }
    int iter = 0;
    float d = sdf.get_distance(p0);
    while (iter < 1000 && d > EPS && d < 1e6)
    {
      p0 += d * dir;
      d = sdf.get_distance(p0);
      iter++;
    }
    if (surface_pos)
      *surface_pos = p0;
    return d <= EPS;
  }

  Texture render_sdf(const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, SDFRenderMode mode)
  {
    unsigned char *data = new unsigned char[4*image_w*image_h];
    std::vector<float> res(image_w*image_h, 0);
    render_sdf_to_array(res, sdf, camera, image_w, image_h, spp, mode);
  
    for (int i=0;i<res.size(); i++)
    {
      data[4*i+0] = 255*CLAMP(res[i],0,1);
      data[4*i+1] = 255*CLAMP(res[i],0,1);
      data[4*i+2] = 255*CLAMP(res[i],0,1);
      data[4*i+3] = 255;
    }

    Texture t = engine::textureManager->create_texture(image_w, image_h, GL_RGBA8, 1, data, GL_RGBA);
    delete[] data;
    return t;
  }

  void evaluate_job(const ProceduralSdf &sdf, std::span<float> positions, std::span<float> distances, 
                    unsigned batch_size, unsigned offset, unsigned end)
  {
    while (offset < end)
    {
      unsigned size = MIN(batch_size, end - offset);
      sdf.get_distance_batch(size, positions.data() + 3*offset, distances.data() + offset, nullptr, nullptr);
      offset += size;
    }
  }

  void evaluate(const ProceduralSdf *sdfs, unsigned max_threads, std::span<float> positions, std::span<float> distances, unsigned count)
  {
    assert(positions.size() >= 3*count);
    assert(distances.size() >= count);

    unsigned batch_size = 512;
    unsigned batches = (count+batch_size-1)/batch_size;
    unsigned num_threads = MIN(max_threads, batches/10 + 1); //do not spawn too many tiny threads

    unsigned batches_per_thread = (batches+num_threads-1)/num_threads;

    if (num_threads == 1)
      evaluate_job(sdfs[0], positions, distances, batch_size, 0, count);
    else
    {
      std::vector<std::thread> threads;
      for (int i = 0; i < num_threads; i++)
      {
        unsigned offset = batch_size*batches_per_thread*i;
        unsigned end = (i==num_threads-1) ? count : batch_size*batches_per_thread*(i+1);
        threads.push_back(std::thread(&evaluate_job, sdfs[i], positions, distances, batch_size, offset, end));
      }

      for (auto &t : threads)
        t.join();
    }
  }

  void render_sdf_to_array_batched(std::span<float> out_array, const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, 
                                   int spp, SDFRenderMode mode)
  {
    constexpr unsigned max_threads = 1; //for some reason multithreading is actually slower
    std::vector<ProceduralSdf> sdfs;
    if (max_threads > 1)
    {
      for (int i=0;i<max_threads;i++)
        sdfs.emplace_back(sdf);
    }

    AABB bbox = sdf.get_bbox();
    AABB inflated_bbox = bbox.expand(1.1f);

    glm::mat4 projInv = glm::inverse(camera.get_proj(false));
    glm::mat4 viewInv = glm::inverse(camera.get_view());
    glm::vec3 light_dir = normalize(glm::vec3(1,1,1));

    int spp_a = MAX(1,floor(sqrtf(spp)));
    float mult = 1.0f/(spp_a*spp_a);

    std::vector<int> hits(image_w*image_h, 0); //one per pixel
    std::vector<int> pixel_indices(image_w*image_h*spp_a*spp_a);
    std::vector<glm::vec3> dirs(image_w*image_h*spp_a*spp_a);
    std::vector<float> points(3*image_w*image_h*spp_a*spp_a);
    std::vector<float> distances(image_w*image_h*spp_a*spp_a, 0);
    std::vector<float> fd_points, fd_distances, fd_pixel_indices;
    if (mode == SDFRenderMode::LAMBERT)
    {
      fd_points.resize(6*3*image_w*image_h*spp_a*spp_a);
      fd_distances.resize(6*image_w*image_h*spp_a*spp_a);
      fd_pixel_indices.resize(image_w*image_h*spp_a*spp_a);
    }
    
    int index = 0;
    for (int yi=0;yi<image_h;yi++)
    {
      for (int xi=0;xi<image_w;xi++)
      {
        for (int yp=0;yp<spp_a;yp++)
        {
          for (int xp=0;xp<spp_a;xp++)
          {
            float y = (float)(yi*spp_a+yp)/(image_h*spp_a);
            float x = (float)(xi*spp_a+xp)/(image_w*spp_a);
            glm::vec3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            glm::vec3 p0 = camera.origin;
            float t = 0;
            if (bbox.contains(p0) || bbox.intersects(p0, dir, &t))
            {
              p0 += t*dir;
              pixel_indices[index] = yi*image_w+xi;
              dirs[index] = dir;
              points[3*index + 0] = p0.x;
              points[3*index + 1] = p0.y;
              points[3*index + 2] = p0.z;
              index++;
            }
            //else ray won't intersect SDF
          }
        }
      }
    }

    constexpr int max_steps = 1000;
    constexpr float EPS = 3e-6;
    constexpr float h = 0.001;

    int points_left = index;
    int points_found = 0;
    int iteration = 0;
    while (iteration < max_steps && points_left > 0)
    {
      evaluate(sdfs.empty() ? &sdf : sdfs.data(), max_threads, points, distances, points_left);

      int new_points_left = points_left;
      int free_pos = 0;
      for (int i=0;i<points_left;i++)
      {
        float d = distances[i];
        glm::vec3 p0 = glm::vec3(points[3*i+0], points[3*i+1], points[3*i+2]);
        if (!inflated_bbox.contains(p0))//ray went too far, no hit
        {
          distances[i] = 2e6;
          new_points_left--;
        }
        else if (d < EPS) //found hit!
        {
          distances[i] = 2e6;
          new_points_left--;
            if (mode == SDFRenderMode::LAMBERT)
            {
              std::vector<float> np = {points[3*i+0]+h, points[3*i+1]  , points[3*i+2],
                                       points[3*i+0]-h, points[3*i+1]  , points[3*i+2],
                                       points[3*i+0]  , points[3*i+1]+h, points[3*i+2],
                                       points[3*i+0]  , points[3*i+1]-h, points[3*i+2],
                                       points[3*i+0]  , points[3*i+1]  , points[3*i+2]+h,
                                       points[3*i+0]  , points[3*i+1]  , points[3*i+2]-h};
              for (int k=0;k<3*6;k++)
                fd_points[3*6*points_found + k] = np[k];
              fd_pixel_indices[points_found] = pixel_indices[i];
            }
            else if (mode == SDFRenderMode::MASK)
            {
              out_array[pixel_indices[i]] += mult;
            }
            else if (mode == SDFRenderMode::DEPTH)
            {
              float z = glm::length(p0 - camera.origin);
              float d = (1 / z - 1 / camera.z_near) / (1 / camera.z_far - 1 / camera.z_near);
              out_array[pixel_indices[i]] += mult*d;
            }
            else if (mode == SDFRenderMode::LINEAR_DEPTH)
            {
              float z = glm::length(p0 - camera.origin);
              out_array[pixel_indices[i]] += mult*((z - camera.z_near) / (camera.z_far - camera.z_near));
            }
            else if (mode == SDFRenderMode::INVERSE_LINEAR_DEPTH)
            {
              float z = glm::length(p0 - camera.origin);
              out_array[pixel_indices[i]] += mult*(1 - (z - camera.z_near) / (camera.z_far - camera.z_near));
            }
            points_found++;
        }
        else //ray marching should continue
        {
          while (free_pos < i && distances[free_pos] < 1e6)
            free_pos++;
          
          points[3*free_pos + 0] = points[3*i + 0] + d*dirs[i].x;
          points[3*free_pos + 1] = points[3*i + 1] + d*dirs[i].y;
          points[3*free_pos + 2] = points[3*i + 2] + d*dirs[i].z;
          if (free_pos != i)
          {
            pixel_indices[free_pos] = pixel_indices[i];
            dirs[free_pos] = dirs[i];
            distances[free_pos] = distances[i];
            distances[i] = 2e6;
          }
        }
      }
      iteration++;
      points_left = new_points_left;
      //logerr("iteration %d, points left %d/%d", iteration, points_left, index);
    }

    if (mode == SDFRenderMode::LAMBERT)
    {
      evaluate(sdfs.empty() ? &sdf : sdfs.data(), max_threads, fd_points, fd_distances, 6*points_found);
      for (int i=0;i<points_found;i++)
      {
        float ddx = (fd_distances[6*i+0] - fd_distances[6*i+1])/(2*h);
        float ddy = (fd_distances[6*i+2] - fd_distances[6*i+3])/(2*h);
        float ddz = (fd_distances[6*i+4] - fd_distances[6*i+5])/(2*h);
        glm::vec3 n = glm::normalize(glm::vec3(ddx, ddy, ddz));
        out_array[fd_pixel_indices[i]] += mult*MAX(0.1f, dot(n, light_dir));
      }
    }
  }

  void render_sdf_to_array(std::span<float> out_array, const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, SDFRenderMode mode)
  {
    render_sdf_to_array_batched(out_array, sdf, camera, image_w, image_h, spp, mode); 
  }

  void sdf_to_outer_point_cloud(const ProceduralSdf &sdf, int points_count, std::vector<glm::vec3> *points, 
                                std::vector<glm::vec3> *outside_points)
  {
    //cast rays from large sphere outside of model to (0,0,0)
    //This method can't be an accurate representation of surface
    //as points won't be uniformly distributed.
    assert(points != nullptr);
    float r = 10000;

    AABB sdf_bbox = sdf.get_bbox();
    *points = {};
    points->reserve(points_count);
    for (int i=0;i<points_count;i++)
    {
      float phi = urand(0, 2*PI);
      float psi = acos(1 - 2*urand());
      if (urand() < 0.5)
        psi = -psi;
      glm::vec3 start_pos = r*glm::vec3{sin(psi)*cos(phi), cos(psi), sin(psi)*sin(phi)};
      glm::vec3 dir = glm::normalize(-start_pos); //tracing rays from random points to (0,0,0)
      glm::vec3 p0;
      if (sdf_sphere_tracing(sdf, sdf_bbox, start_pos, dir, &p0))
        points->push_back(p0);
    }

    if (outside_points)
    {
      AABB bbox = sdf.get_bbox();
      AABB inflated_bbox = AABB(bbox.min_pos - glm::vec3(0.01,0.01,0.01), bbox.max_pos + glm::vec3(0.01,0.01,0.01));
      
      *outside_points = {};
      outside_points->reserve(points_count);
      while (outside_points->size() < points_count)
      {
        glm::vec3 p = glm::vec3(urand(inflated_bbox.min_pos.x, inflated_bbox.max_pos.x),
                                urand(inflated_bbox.min_pos.y, inflated_bbox.max_pos.y),
                                urand(inflated_bbox.min_pos.z, inflated_bbox.max_pos.z));
        if (sdf.get_distance(p) > 0.01)
          outside_points->push_back(p);
      }
    }
  }

  void sdf_to_point_cloud_with_dist(const ProceduralSdf &sdf, int points_count, std::vector<glm::vec3> *points, 
                                    std::vector<float> *distances)
  {
    *points = {};
    points->resize(points_count);
    *distances = {};
    distances->resize(points_count);
    
    AABB bbox = sdf.get_bbox();
    glm::vec3 center = 0.5f*(bbox.max_pos + bbox.min_pos);
    glm::vec3 size = 0.5f*(bbox.max_pos - bbox.min_pos);
    AABB inflated_bbox = AABB(center - 2.0f*size, center + 2.0f*size);

    for (int i=0;i<points_count;i++)
    {
              glm::vec3 p = glm::vec3(urand(inflated_bbox.min_pos.x, inflated_bbox.max_pos.x),
                                urand(inflated_bbox.min_pos.y, inflated_bbox.max_pos.y),
                                urand(inflated_bbox.min_pos.z, inflated_bbox.max_pos.z));
      (*points)[i] = p;
      (*distances)[i] = sdf.get_distance(p);
    }
  }

  void sdf_to_point_cloud(const ProceduralSdf &sdf, int points_count, std::vector<glm::vec3> *points, 
                          std::vector<glm::vec3> *outside_points)
  {
    //More accurate ay to get point clouds from an arbitrary SDF.
    //It has two steps: first we estimate size of SDF (as bounding box)
    //Then we are tracing rays from random points outside of model in random
    //direction and store all hits. I haven't tried to prove it, but 
    //likely it will give more or less random distribution of points on 
    //relatively simple surfaces.
    assert(points != nullptr);
    float r = 10000;

    AABB sdf_bbox = sdf.get_bbox();
    *points = {};
    points->reserve(points_count);

    std::vector<glm::vec3> estimate_points;
    sdf_to_outer_point_cloud(sdf, 5000, &estimate_points, nullptr);
    AABB bbox = get_point_cloud_bbox(estimate_points);
    glm::vec3 center = 0.5f*(bbox.max_pos + bbox.min_pos);
    glm::vec3 size = 0.5f*(bbox.max_pos - bbox.min_pos);
    AABB inflated_bbox = AABB(center - 2.0f*size, center + 2.0f*size);

    int tries = 100*points_count;
    int iter = 0;
    while (iter<tries && points->size()<points_count)
    {
      glm::vec3 p = glm::vec3(urand(inflated_bbox.min_pos.x, inflated_bbox.max_pos.x),
                              urand(inflated_bbox.min_pos.y, inflated_bbox.max_pos.y),
                              urand(inflated_bbox.min_pos.z, inflated_bbox.max_pos.z));
      if (sdf.get_distance(p) > 0.01)
      {
        float phi = urand(0, 2*PI);
        float psi = acos(1 - 2*urand());
        if (urand() < 0.5)
          psi = -psi;
        glm::vec3 dir = glm::vec3{sin(psi)*cos(phi), cos(psi), sin(psi)*sin(phi)};
        glm::vec3 p0;
        if (sdf_sphere_tracing(sdf, sdf_bbox, p, dir, &p0))
          points->push_back(p0);
      }
      iter++;
    }

    if (outside_points)
    {
      *outside_points = {};
      outside_points->reserve(points_count);
      while (outside_points->size() < points_count)
      {
        glm::vec3 p = glm::vec3(urand(inflated_bbox.min_pos.x, inflated_bbox.max_pos.x),
                                urand(inflated_bbox.min_pos.y, inflated_bbox.max_pos.y),
                                urand(inflated_bbox.min_pos.z, inflated_bbox.max_pos.z));
        if (sdf.get_distance(p) > 0.01)
          outside_points->push_back(p);
      }
    }
  }

  AABB get_point_cloud_bbox(const std::vector<glm::vec3> &points)
  {
    glm::vec3 minv(1e9,1e9,1e9);
    glm::vec3 maxv(-1e9,-1e9,-1e9);
    for (auto &p : points)
    {
      minv = glm::min(minv, p);
      maxv = glm::max(maxv, p);
    }

    return AABB(minv, maxv);
  }

  float get_sdf_image_based_quality(const ProceduralSdf &reference_sdf, const ProceduralSdf &sdf)
  {
    int image_size = 512;
    int points_count = 1000;
    std::vector<glm::vec3> points;
    sdf_to_point_cloud(reference_sdf, points_count, &points);
    AABB bbox = get_point_cloud_bbox(points);
    float d = 1.5*length(bbox.max_pos - bbox.min_pos);
    CameraSettings cam;
    cam.target = 0.5f*(bbox.min_pos + bbox.max_pos);
    cam.origin = cam.target + glm::vec3(0,0,-d);
    cam.up = glm::vec3(0,1,0);

    auto cameras = get_cameras_uniform_sphere(cam, 64, d);

    float diff = 0.0;
    float mse = 0.0;
    for (auto &cam : cameras)
    {
      Texture t1 = render_sdf(reference_sdf, cam, image_size, image_size, 4, SDFRenderMode::MASK);
      Texture t2 = render_sdf(sdf, cam, image_size, image_size, 4, SDFRenderMode::MASK);
      mse += ImageMetric::get(t1, t2);
    }
    return -10*log10(MAX(1e-9f,mse/cameras.size()));
  }

  float get_sdf_similarity_MSE(const ProceduralSdf &reference_sdf, const ProceduralSdf &sdf)
  {
    int points = 50000;
    AABB bbox = AABB(min(reference_sdf.get_bbox().min_pos, sdf.get_bbox().min_pos),
                     max(reference_sdf.get_bbox().max_pos, sdf.get_bbox().max_pos));
    double dist_sum = 0;
    for (int i=0;i<points;i++)
    {
      glm::vec3 p = glm::vec3(urand(bbox.min_pos.x, bbox.max_pos.x),
                              urand(bbox.min_pos.y, bbox.max_pos.y),
                              urand(bbox.min_pos.z, bbox.max_pos.z));
      dist_sum += SQR(sdf.get_distance(p) - reference_sdf.get_distance(p));
    }
    return dist_sum/points;
  }

  float render_bbox(const AABB &bbox, glm::vec3 p, float line_q)
  {
    glm::vec3 pb = (p-bbox.min_pos)/bbox.size();
    glm::vec3 pc1 = max(abs(1.0f - line_q*pb), 0.0f);
    glm::vec3 pc2 = max(abs(1.0f - line_q*(1.0f-pb)), 0.0f);
    glm::vec3 pm = max(pc1, pc2);
    return std::max(pm.x, std::max(pm.y, pm.z));
  }

  Texture render_sdf(const SdfScene &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, SDFRenderMode mode)
  {
    unsigned char *data = new unsigned char[4*image_w*image_h];
    std::vector<float> res(image_w*image_h, 0);
    render_sdf_to_array(res, sdf, camera, image_w, image_h, spp, mode);
  
    for (int i=0;i<res.size(); i++)
    {
      data[4*i+0] = 255*CLAMP(res[i],0,1);
      data[4*i+1] = 255*CLAMP(res[i],0,1);
      data[4*i+2] = 255*CLAMP(res[i],0,1);
      data[4*i+3] = 255;
    }

    Texture t = engine::textureManager->create_texture(image_w, image_h, GL_RGBA8, 1, data, GL_RGBA);
    delete[] data;
    return t;
  }

  void render_sdf_to_array(std::span<float> out_array, const SdfScene &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, SDFRenderMode mode)
  {
    assert(sdf.conjunctions.size()>0 && sdf.objects.size()>0);
    AABB sdf_bbox = conv(sdf.conjunctions[0].bbox);
    for (auto &obj : sdf.conjunctions)
    {
      sdf_bbox.min_pos = min(sdf_bbox.min_pos, conv(obj.bbox.min_pos));
      sdf_bbox.max_pos = max(sdf_bbox.max_pos, conv(obj.bbox.max_pos));
    }
    std::vector<glm::vec3> edges =
    {
      {0,0,0},{0,0,1},
      {0,0,0},{0,1,0},
      {0,0,0},{1,0,0},
      {0,0,1},{0,1,1},

      {0,0,1},{1,0,1},
      {0,1,0},{0,1,1},
      {0,1,0},{1,1,0},
      {0,1,1},{1,1,1},

      {1,0,0},{1,1,0},
      {1,0,0},{1,0,1},
      {1,0,1},{1,1,1},
      {1,1,0},{1,1,1},
    };
    std::vector<AABB> boxes = {sdf_bbox};
    for (auto &obj : sdf.conjunctions)
     boxes.push_back(conv(obj.bbox));
    std::vector<AABB> debug_lines;
    float th = 0.01;
    for (auto &box : boxes)
    {
      for (int i=0;i<edges.size();i+=2)
      {
        glm::vec3 p1 = edges[i]*box.min_pos + (1.0f-edges[i])*box.max_pos;
        glm::vec3 p2 = edges[i+1]*box.min_pos + (1.0f-edges[i+1])*box.max_pos;
        debug_lines.push_back(AABB(min(p1,p2)-th,max(p1,p2)+th));
      }
    }
    ///sdf_bbox = AABB({-5,-5,-5},{5,5,5});
    glm::mat4 projInv = glm::inverse(glm::perspective(camera.fov_rad, 1.0f, camera.z_near, camera.z_far));
    glm::mat4 viewInv = glm::inverse(camera.get_view());
    //set light somewhere to the side 
    glm::vec3 light_dir = normalize(glm::vec3(1,1,1));
    int spp_a = MAX(1,floor(sqrtf(spp)));

    #pragma omp parallel for
    for (int yi=0;yi<image_h;yi++)
    {
      std::vector<float> cur_grad;
      std::vector<float> ddist_dpos = {0,0,0};
      for (int xi=0;xi<image_w;xi++)
      {
        float color = 0;
        for (int yp=0;yp<spp_a;yp++)
        {
          for (int xp=0;xp<spp_a;xp++)
          {
            float y = (float)(yi*spp_a+yp)/(image_h*spp_a);
            float x = (float)(xi*spp_a+xp)/(image_w*spp_a);
            glm::vec3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            glm::vec3 p0;

            bool debug_box_found = false;
            if (false)
            {
            for (auto &box : debug_lines)
            {
              if (box.intersects(camera.origin, dir))
              {
                debug_box_found = true;
                color += 1;
                break;
              }
            }
            }
            
            LiteMath::float3 lmp0;
            if (!debug_box_found && sdf_sphere_tracing(sdf, conv(sdf_bbox), conv(camera.origin), conv(dir), &lmp0))
            {
              p0 = conv(lmp0);
              if (mode == SDFRenderMode::MASK)
                color += 1;
              else if (mode == SDFRenderMode::LAMBERT)
              {
                constexpr float h = 0.001;
                cur_grad.clear();
                float ddx = (eval_dist_scene(sdf, lmp0 + LiteMath::float3(h,0,0)) - eval_dist_scene(sdf, lmp0 + LiteMath::float3(-h,0,0)))/(2*h);
                float ddy = (eval_dist_scene(sdf, lmp0 + LiteMath::float3(0,h,0)) - eval_dist_scene(sdf, lmp0 + LiteMath::float3(0,-h,0)))/(2*h);
                float ddz = (eval_dist_scene(sdf, lmp0 + LiteMath::float3(0,0,h)) - eval_dist_scene(sdf, lmp0 + LiteMath::float3(0,0,-h)))/(2*h);
                glm::vec3 n = glm::normalize(glm::vec3(ddx, ddy, ddz));
                color += MAX(0.1f, dot(n, light_dir));
              }
              else if (mode == SDFRenderMode::DEPTH)
              {
                float z = glm::length(p0 - camera.origin);
                float d = (1/z - 1/camera.z_near)/(1/camera.z_far - 1/camera.z_near);
                color += d;
              }
              else if (mode == SDFRenderMode::LINEAR_DEPTH)
              {
                float z = glm::length(p0 - camera.origin);
                color += (z - camera.z_near)/(camera.z_far - camera.z_near);
              }
              else if (mode == SDFRenderMode::INVERSE_LINEAR_DEPTH)
              {
                float z = glm::length(p0 - camera.origin);
                color += 1 - (z - camera.z_near)/(camera.z_far - camera.z_near);
              }
            }
          }
        }
        out_array[yi*image_w+xi] = color/SQR(spp_a);
      }
    }
  }
}