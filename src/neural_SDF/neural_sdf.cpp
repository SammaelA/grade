#include "neural_sdf.h"
#include "private_camera.h"
#include "common_utils/bbox.h"
#include "tinyEngine/engine.h"
#include "common_utils/distribution.h"
#include "neuralCore/siren.h"

#include <vector>
#include <fstream>

namespace nsdf
{
  struct PrimitiveSDF
  {
    float get_distance(glm::vec3 pos) const
    {
      float dist = 1e9;
      for (auto &c : circles)
      {
        glm::vec3 d = pos - glm::vec3(c.x, c.y, c.z);
        dist = std::min(dist, length(d) - c.w);
      }
      return dist;
    }

    AABB get_bbox() const
    {
      if (circles.empty())
        return AABB({-1,-1,-1},{1,1,1});
      
      glm::vec3 min_v(circles[0].x, circles[0].y, circles[0].z);
      glm::vec3 max_v(circles[0].x, circles[0].y, circles[0].z);
      for (auto &c : circles)
      {
        min_v = min(min_v, glm::vec3(c.x-c.w, c.y-c.w, c.z-c.w));
        max_v = max(max_v, glm::vec3(c.x+c.w, c.y+c.w, c.z+c.w));
      }
      return AABB(min_v,max_v);
    }

    std::vector<glm::vec4> circles;
  };

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
  bool primitive_sdf_sphere_tracing(const PrimitiveSDF &sdf, const AABB &sdf_bbox, const glm::vec3 &start_pos, const glm::vec3 &dir, glm::vec3 *surface_pos = nullptr)
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

  Texture render_primitive_sdf(const PrimitiveSDF &sdf, const CameraSettings &camera, glm::vec3 light_dir,
                               int image_w, int image_h, int spp, bool lambert)
  {
    AABB sdf_bbox = sdf.get_bbox();
    glm::mat4 projInv = glm::inverse(camera.get_proj());
    glm::mat4 viewInv = glm::inverse(camera.get_view());
    //set light somewhere to the side 
    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    #pragma omp parallel for
    for (int yi=0;yi<image_h;yi++)
    {
      for (int xi=0;xi<image_w;xi++)
      {
        glm::vec3 color = {0,0,0};
        for (int yp=0;yp<spp_a;yp++)
        {
          for (int xp=0;xp<spp_a;xp++)
          {
            float y = (float)(yi*spp_a+yp)/(image_h*spp_a);
            float x = (float)(xi*spp_a+xp)/(image_w*spp_a);
            glm::vec3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            glm::vec3 p0;
            
            if (primitive_sdf_sphere_tracing(sdf, sdf_bbox, camera.origin, dir, &p0))
            {
              if (lambert)
              {
                constexpr float h = 0.001;
                float ddx = (sdf.get_distance(p0 + glm::vec3(h,0,0)) - sdf.get_distance(p0 + glm::vec3(-h,0,0)))/(2*h);
                float ddy = (sdf.get_distance(p0 + glm::vec3(0,h,0)) - sdf.get_distance(p0 + glm::vec3(0,-h,0)))/(2*h);
                float ddz = (sdf.get_distance(p0 + glm::vec3(0,0,h)) - sdf.get_distance(p0 + glm::vec3(0,0,-h)))/(2*h);
                glm::vec3 n = glm::normalize(glm::vec3(ddx, ddy, ddz));
                color += glm::vec3(1,1,1) * MAX(0.1f, dot(n, light_dir));
              }
              else
                color += glm::vec3(1,1,1);
            }
          }
        }
        data[4*(yi*image_w+xi)+0] = 255*(color.x/SQR(spp_a));
        data[4*(yi*image_w+xi)+1] = 255*(color.y/SQR(spp_a));
        data[4*(yi*image_w+xi)+2] = 255*(color.z/SQR(spp_a));
        data[4*(yi*image_w+xi)+3] = 255;
      }
    }

    Texture t = engine::textureManager->create_texture(image_w, image_h, GL_RGBA8, 1, data, GL_RGBA);
    delete[] data;
    return t;
  }

  void task_1_create_references()
  {
    CameraSettings cam;
    cam.origin = glm::vec3(0,0,3);
    cam.target = glm::vec3(0,0,0);
    cam.up = glm::vec3(0,1,0);

    glm::vec3 light_dir = normalize(cam.origin + glm::vec3(cam.origin.z, cam.origin.y, cam.origin.x) - cam.target);
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};
    l.to_file("saves/task1_references/light.txt");

    CameraSettings cam1 = cam;
    cam1.origin = glm::vec3(0,3*sin(0),3*cos(0));
    convert(cam1).to_file("saves/task1_references/cam1.txt");

    CameraSettings cam2 = cam;
    cam2.origin = glm::vec3(0,3*sin(2*PI/3),3*cos(2*PI/3));
    convert(cam2).to_file("saves/task1_references/cam2.txt");

    CameraSettings cam3 = cam;
    cam3.origin = glm::vec3(0,3*sin(4*PI/3),3*cos(4*PI/3));
    convert(cam3).to_file("saves/task1_references/cam3.txt");

    PrimitiveSDF sdf1;
    sdf1.circles.push_back({0,0,0,0.7f});

    PrimitiveSDF sdf2;
    sdf2.circles.push_back({0,0.3f,0.4f,0.4f});
    sdf2.circles.push_back({-0.3f,0,0.1f,0.3f});
    sdf2.circles.push_back({0,-0.2f,0,0.25f});
    Texture t;
    
    t = render_primitive_sdf(sdf1, cam1, light_dir, 512, 512, 16, false);
    engine::textureManager->save_png(t, "task1_references/sdf1_cam1_reference");
    t = render_primitive_sdf(sdf1, cam2, light_dir, 512, 512, 16, false);
    engine::textureManager->save_png(t, "task1_references/sdf1_cam2_reference");
    t = render_primitive_sdf(sdf1, cam3, light_dir, 512, 512, 16, false);
    engine::textureManager->save_png(t, "task1_references/sdf1_cam3_reference");

    t = render_primitive_sdf(sdf2, cam1, light_dir, 512, 512, 16, false);
    engine::textureManager->save_png(t, "task1_references/sdf2_cam1_reference");
    t = render_primitive_sdf(sdf2, cam2, light_dir, 512, 512, 16, false);
    engine::textureManager->save_png(t, "task1_references/sdf2_cam2_reference");
    t = render_primitive_sdf(sdf2, cam3, light_dir, 512, 512, 16, false);
    engine::textureManager->save_png(t, "task1_references/sdf2_cam3_reference");
  
    t = render_primitive_sdf(sdf1, cam1, light_dir, 512, 512, 16, true);
    engine::textureManager->save_png(t, "task1_references/sdf1_cam1_lambert_reference");
    t = render_primitive_sdf(sdf1, cam2, light_dir, 512, 512, 16, true);
    engine::textureManager->save_png(t, "task1_references/sdf1_cam2_lambert_reference");
    t = render_primitive_sdf(sdf1, cam3, light_dir, 512, 512, 16, true);
    engine::textureManager->save_png(t, "task1_references/sdf1_cam3_lambert_reference");

    t = render_primitive_sdf(sdf2, cam1, light_dir, 512, 512, 16, true);
    engine::textureManager->save_png(t, "task1_references/sdf2_cam1_lambert_reference");
    t = render_primitive_sdf(sdf2, cam2, light_dir, 512, 512, 16, true);
    engine::textureManager->save_png(t, "task1_references/sdf2_cam2_lambert_reference");
    t = render_primitive_sdf(sdf2, cam3, light_dir, 512, 512, 16, true);
    engine::textureManager->save_png(t, "task1_references/sdf2_cam3_lambert_reference");
  }

  void load_points_cloud(std::string filename, std::vector<float> *points, std::vector<float> *distances)
  {
    std::ifstream in(filename, std::ios_base::binary);
    assert(in.is_open());
    int count = -1;
    in.read(reinterpret_cast<char*>(&count), sizeof(int));
    assert(count > 0 && count < (1<<16));
    points->resize(3*count);
    distances->resize(count);
    in.read(reinterpret_cast<char*>(points->data()), sizeof(float)*points->size());
    in.read(reinterpret_cast<char*>(distances->data()), sizeof(float)*distances->size());
    in.close();
  }

  void save_points_cloud(std::string filename, const std::vector<float> &points, const std::vector<float> &distances)
  {
    assert(points.size() > 0 && points.size() == 3*distances.size());
    std::ofstream out(filename, std::ios_base::binary);
    assert(out.is_open());
    int count = distances.size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(int));
    out.write(reinterpret_cast<const char*>(points.data()), sizeof(float)*points.size());
    out.write(reinterpret_cast<const char*>(distances.data()), sizeof(float)*distances.size());
    out.close();
  }

  void create_point_cloud_spheres(int count, AABB bbox, std::vector<float> *points, std::vector<float> *distances)
  {
    PrimitiveSDF sdf2;
    sdf2.circles.push_back({0,0.3f,0.4f,0.4f});
    sdf2.circles.push_back({-0.3f,0,0.1f,0.3f});
    sdf2.circles.push_back({0,-0.2f,0,0.25f});

    points->resize(3*count);
    distances->resize(count);

    for (int i=0;i<count;i++)
    {
      glm::vec3 p = glm::vec3(urand(bbox.min_pos.x, bbox.max_pos.x),
                              urand(bbox.min_pos.y, bbox.max_pos.y),
                              urand(bbox.min_pos.z, bbox.max_pos.z));
      (*points)[3*i+0] = p.x;
      (*points)[3*i+1] = p.y;
      (*points)[3*i+2] = p.z;
      (*distances)[i] = sdf2.get_distance(p);
    }
  }

  Texture render_neural_sdf(nn::Siren &sdf, AABB bbox, const CameraSettings &camera, 
                            int image_w, int image_h, int spp, bool lambert, glm::vec3 light_dir)
  {
    glm::vec3 center = 0.5f*(bbox.max_pos + bbox.min_pos);
    glm::vec3 size = 0.5f*(bbox.max_pos - bbox.min_pos);
    AABB inflated_bbox = AABB(center - 1.1f*size, center + 1.1f*size);

    glm::mat4 projInv = glm::inverse(camera.get_proj());
    glm::mat4 viewInv = glm::inverse(camera.get_view());

    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    std::vector<int> hits(image_w*image_h, 0); //one per pixel
    std::vector<glm::vec3> colors(image_w*image_h, glm::vec3(0,0,0));
    std::vector<int> pixel_indices(image_w*image_h*spp_a*spp_a);
    std::vector<glm::vec3> dirs(image_w*image_h*spp_a*spp_a);
    std::vector<float> points(3*image_w*image_h*spp_a*spp_a);
    std::vector<float> distances(image_w*image_h*spp_a*spp_a, 0);
    std::vector<float> fd_points, fd_distances, fd_pixel_indices;
    if (lambert)
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
      sdf.evaluate(points, distances, points_left);

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
            hits[pixel_indices[i]]++;
            if (lambert)
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
    }

    if (!lambert)
    {
      for (int i=0;i<image_w*image_h;i++)
        colors[i] += glm::vec3(hits[i], hits[i], hits[i]);
    }
    else
    {
      sdf.evaluate(fd_points, fd_distances, 6*points_found);
      for (int i=0;i<points_found;i++)
      {
        float ddx = (fd_distances[6*i+0] - fd_distances[6*i+1])/(2*h);
        float ddy = (fd_distances[6*i+2] - fd_distances[6*i+3])/(2*h);
        float ddz = (fd_distances[6*i+4] - fd_distances[6*i+5])/(2*h);
        glm::vec3 n = glm::normalize(glm::vec3(ddx, ddy, ddz));
        colors[fd_pixel_indices[i]] += glm::vec3(1,1,1) * MAX(0.1f, dot(n, light_dir));
      }
    }

    for (int i = 0; i < image_w * image_h; i++)
    {
      glm::vec3 color = colors[i];
      data[4 * i + 0] = 255 * (color.x / SQR(spp_a));
      data[4 * i + 1] = 255 * (color.y / SQR(spp_a));
      data[4 * i + 2] = 255 * (color.z / SQR(spp_a));
      data[4 * i + 3] = 255;
    }

    Texture t = engine::textureManager->create_texture(image_w, image_h, GL_RGBA8, 1, data, GL_RGBA);
    delete[] data;
    return t;
  }

  void task_2_create_references()
  {
    AABB bbox({-1,-1,-1},{1,1,1});

    CameraSettings cam;
    cam.origin = glm::vec3(0,0,3);
    cam.target = glm::vec3(0,0,0);
    cam.up = glm::vec3(0,1,0);

    glm::vec3 light_dir = normalize(cam.origin + glm::vec3(cam.origin.z, cam.origin.y, cam.origin.x) - cam.target);
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};
    l.to_file("saves/task2_references/light.txt");

    CameraSettings cam1 = cam;
    cam1.origin = glm::vec3(0,3*sin(0),3*cos(0));
    convert(cam1).to_file("saves/task2_references/cam1.txt");

    CameraSettings cam2 = cam;
    cam2.origin = glm::vec3(0,3*sin(2*PI/3),3*cos(2*PI/3));
    convert(cam2).to_file("saves/task2_references/cam2.txt");

    CameraSettings cam3 = cam;
    cam3.origin = glm::vec3(0,3*sin(4*PI/3),3*cos(4*PI/3));
    convert(cam3).to_file("saves/task2_references/cam3.txt");

    {
      nn::TensorProcessor::init("GPU");
      std::vector<float> points, distances;
      create_point_cloud_spheres(5000, bbox, &points, &distances);
      save_points_cloud("saves/task2_references/points.bin", points, distances);
      nn::Siren network(nn::Siren::Type::SDF, 3, 64);
      network.train(points, distances, 512, 15000);
      network.save_weights_to_file("saves/task2_references/weights.bin");
    }
    {
      std::vector<float> points, distances;
      load_points_cloud("saves/task2_references/points.bin", &points, &distances);
      nn::Siren network(nn::Siren::Type::SDF, 3, 64);
      network.initialize_from_file("saves/task2_references/weights.bin");
      
      Texture t;
      t = render_neural_sdf(network, bbox, cam1, 512, 512, 16, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf1_cam1_reference");
      t = render_neural_sdf(network, bbox, cam2, 512, 512, 16, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf1_cam2_reference");
      t = render_neural_sdf(network, bbox, cam3, 512, 512, 16, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf1_cam3_reference");
    }
  }

  void neural_SDF_test()
  {
    //task_1_create_references();
    task_2_create_references();
  }
}