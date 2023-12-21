#include "sdf_rendering.h"
#include "graphics_utils/image_metrics.h"
#include "tinyEngine/engine.h"
#include "common_utils/distribution.h"
#include "preprocessing.h"

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
  bool sdf_sphere_tracing(const ProceduralSdf &sdf, const glm::vec3 &start_pos, const glm::vec3 &dir, glm::vec3 *surface_pos = nullptr)
  {
    constexpr float EPS = 3*1e-6;
    glm::vec3 p0 = start_pos;
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

  struct SDFPartInfo
  {
    UPGPart part;
    AABB bbox;
    const SdfNode *part_root;
  };

  std::vector<SDFPartInfo> get_sdf_parts_info(const ProceduralSdf &sdf)
  {
    auto parts = get_sdf_parts(sdf.structure);
    std::vector<SDFPartInfo> parts_info(parts.size());

    for (int i=0;i<parts.size();i++)
    {
      parts_info[i].part = parts[i];
      parts_info[i].part_root = sdf.all_nodes[parts_info[i].part.s_range.first].get();
      parts_info[i].bbox = parts_info[i].part_root->get_bbox();
    }

    return parts_info;
  }

  Texture render_sdf(const ProceduralSdf &sdf, const CameraSettings &camera, int image_w, int image_h, int spp, bool lambert)
  {
    glm::mat4 projInv = glm::inverse(camera.get_proj());
    glm::mat4 viewInv = glm::inverse(camera.get_view());
    //set light somewhere to the side 
    glm::vec3 light_dir = normalize(camera.origin + glm::vec3(camera.origin.z, camera.origin.y, camera.origin.x) - camera.target);
    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    #pragma omp parallel for
    for (int yi=0;yi<image_h;yi++)
    {
      std::vector<float> cur_grad;
      std::vector<float> ddist_dpos = {0,0,0};
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
            
            if (sdf_sphere_tracing(sdf, camera.origin, dir, &p0))
            {
              if (lambert)
              {
                cur_grad.clear();
                sdf.get_distance(p0, &cur_grad, &ddist_dpos);
                glm::vec3 n = glm::normalize(glm::vec3(ddist_dpos[0], ddist_dpos[1], ddist_dpos[2]));
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

  void sdf_to_outer_point_cloud(const ProceduralSdf &sdf, int points_count, std::vector<glm::vec3> *points, 
                                std::vector<glm::vec3> *outside_points)
  {
    //cast rays from large sphere outside of model to (0,0,0)
    //This method can't be an accurate representation of surface
    //as points won't be uniformly distributed.
    assert(points != nullptr);
    float r = 10000;

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
      if (sdf_sphere_tracing(sdf, start_pos, dir, &p0))
        points->push_back(p0);
    }

    if (outside_points)
    {
      AABB bbox = sdf.root->get_bbox();
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
    
    AABB bbox = sdf.root->get_bbox();
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
        if (sdf_sphere_tracing(sdf, p, dir, &p0))
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

  float get_sdf_image_based_quality(ProceduralSdf reference_sdf, ProceduralSdf sdf)
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
      Texture t1 = render_sdf(reference_sdf, cam, image_size, image_size, 4, false);
      Texture t2 = render_sdf(sdf, cam, image_size, image_size, 4, false);
      mse += ImageMetric::get(t1, t2);
    }
    return -10*log10(MAX(1e-9f,mse/cameras.size()));
  }
}