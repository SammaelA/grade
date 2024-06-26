#include "neural_sdf.h"
#include "private_camera.h"
#include "public_scene.h"
#include "common_utils/bbox.h"
#include "tinyEngine/engine.h"
#include "tinyEngine/model.h"
#include "graphics_utils/modeling.h"
#include "common_utils/bvh.h"
#include "common_utils/distribution.h"
#include "neuralCore/siren.h"
#include "graphics_utils/render_point_cloud.h"
#include "../UPG_reconstruction/sdf_node.h"

#include <vector>
#include <fstream>

namespace nsdf
{
  float sdSphere(float3 p)
  {
    return length(p) - 1.0f;
  }

  float sdBox(float3 p)
  {
    float3 q = abs(p) - float3(1.0f);
    return length(LiteMath::max(q, 0.0f)) + LiteMath::min(LiteMath::max(q.x, LiteMath::max(q.y, q.z)), 0.0f);
  }

  float sdMengerSponge(float3 p)
  {
    float d = sdBox(p);
    float3 res = float3(d, 1.0f, 0.0f);

    int Iterations = 4;
    float s = 1.0f;
    for (int m = 0; m < Iterations; m++)
    {
      float3 a = mod(p * s, float3(2.0f)) - 1.0f;
      s *= 3.0f;
      float3 r = abs(1.0f - 3.0f * abs(a));

      float da = LiteMath::max(r.x, r.y);
      float db = LiteMath::max(r.y, r.z);
      float dc = LiteMath::max(r.z, r.x);
      float c = (LiteMath::min(da, LiteMath::min(db, dc)) - 1.0f) / s;

      d = std::max(d, c);
    }

    return d;
  }

  float sdMandelbulb(float3 pos) 
  {
    int Iterations = 4;
    float Power = 8;
    float Bailout = 1.15;

    float3 z = pos;
    float dr = 1.0;
    float r = 0.0;
    for (int i = 0; i < Iterations ; i++) 
    {
      r = length(z);
      if (r>Bailout) break;
      
      // convert to polar coordinates
      float theta = acos(z.z/r);
      float phi = atan2(z.y,z.x);
      dr =  pow( r, Power-1.0)*Power*dr + 1.0;
      
      // scale and rotate the point
      float zr = pow( r,Power);
      theta = theta*Power;
      phi = phi*Power;
      
      // convert back to cartesian coordinates
      z = zr*float3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
      z+=pos;
    }
    return 0.5*log(r)*r/dr;
  }

  float sdRecursiveTetrahedron(float3 z)
  {
    int Iterations = 6;
    float Scale = 2.0f;
    float Or = 0.01;
    float3 a1 = float3( 1, 1, 1);
    float3 a2 = float3(-1,-1, 1);
    float3 a3 = float3( 1,-1,-1);
    float3 a4 = float3(-1, 1,-1);
    float3 c;
    int n = 0;
    float dist, d;
    while (n < Iterations) {
      c = a1; dist = length(z-a1);
            d = length(z-a2); if (d < dist) { c = a2; dist=d; }
      d = length(z-a3); if (d < dist) { c = a3; dist=d; }
      d = length(z-a4); if (d < dist) { c = a4; dist=d; }
      z = Scale*z-c*(Scale-1.0f);
      n++;
    }

    return length(z) * pow(Scale, float(-n)) - Or;
  }

    enum PrimitiveType
    {
      SPHERE,
      BOX,
      SPONGE,
      MANDELBULB,
      REC_TETRAHEDRON
    };
    struct Primitive
    {
      Primitive() {};
      Primitive(PrimitiveType _t, float3 _shift, float3 _scale):
      type(_t),
      shift(_shift),
      scale(_scale) {}

      PrimitiveType type;
      float3 shift = float3(0,0,0);
      float3 scale = float3(1,1,1);
    };

  struct PrimitiveSDF
  {
    float get_distance(float3 pos) const
    {
      float dist = 1e9;
      for (auto &p : primitives)
      {
        float3 pt = (pos-p.shift)/p.scale;
        float d = 1e9;
        switch (p.type)
        {
        case SPHERE:
          d = sdSphere(pt);
          break;
        case BOX:
          d = sdBox(pt);
          break;
        case SPONGE:
          d = sdMengerSponge(pt);
          break;      
        case MANDELBULB:
          d = sdMandelbulb(pt);
          break;  
        case REC_TETRAHEDRON:
          d = sdRecursiveTetrahedron(pt);
          break;
        default:
          break;
        }
        d *= p.scale.x;
        dist = std::min(dist, d);
      }
      return dist;
    }

    AABB get_bbox() const
    {
      if (primitives.empty())
        return AABB({-1,-1,-1},{1,1,1});
      
      float3 min_v(1e9,1e9,1e9);
      float3 max_v(-1e9,-1e9,-1e9);
      for (auto &p : primitives)
      {
        min_v = min(min_v, float3(-1,-1,-1)*p.scale + p.shift);
        max_v = max(max_v, float3( 1, 1, 1)*p.scale + p.shift);
      }
      return AABB(min_v,max_v).expand(1.25f);
    }

    std::vector<Primitive> primitives;
  };

  inline float3 EyeRayDirNormalized(float x/*in [0,1]*/, float y/*in [0,1]*/, float4x4 projInv)
  {
    float4 pos = float4(2.0f*x - 1.0f, -2.0f*y + 1.0f, 0.0f, 1.0f);
    pos = projInv * pos;
    pos /= pos.w;
    return normalize(to_float3(pos));
  }

  inline float3 transformRay(float3 ray, float4x4 viewInv)
  {
    float3 p1 = to_float3(viewInv*float4(0,0,0,1));
    float3 p2 = to_float3(viewInv*float4(ray.x,ray.y,ray.z,1));
    return normalize(p2-p1);
  }

  //intersects SDF with ray start_pos + t*dir. Returns if ray intersected with SDF and optionally - point of intersection
  bool primitive_sdf_sphere_tracing(const PrimitiveSDF &sdf, const AABB &sdf_bbox, const float3 &start_pos, const float3 &dir, float3 *surface_pos = nullptr)
  {
    constexpr float EPS = 3*1e-6;
    float3 p0 = start_pos;
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

  enum ColorMode
  {
    GRAY,
    SINPOS,
    NORMALS
  };

  Texture render_primitive_sdf(const PrimitiveSDF &sdf, const CameraSettings &camera, float3 light_dir,
                               int image_w, int image_h, int spp, ColorMode color_mode, bool ao = false, bool soft_shadows = false)
  {
    AABB sdf_bbox = sdf.get_bbox();
    float4x4 projInv = LiteMath::inverse4x4(camera.get_proj(false));
    float4x4 viewInv = LiteMath::inverse4x4(camera.get_view());
    //set light somewhere to the side 
    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    #pragma omp parallel for
    for (int yi=0;yi<image_h;yi++)
    {
      for (int xi=0;xi<image_w;xi++)
      {
        float3 color = {0,0,0};
        for (int yp=0;yp<spp_a;yp++)
        {
          for (int xp=0;xp<spp_a;xp++)
          {
            float y = (float)(yi*spp_a+yp)/(image_h*spp_a);
            float x = (float)(xi*spp_a+xp)/(image_w*spp_a);
            float3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            float3 p0;
            
            if (primitive_sdf_sphere_tracing(sdf, sdf_bbox, camera.origin, dir, &p0))
            {
              constexpr float h = 0.001;
              float ddx = (sdf.get_distance(p0 + float3(h, 0, 0)) - sdf.get_distance(p0 + float3(-h, 0, 0))) / (2 * h);
              float ddy = (sdf.get_distance(p0 + float3(0, h, 0)) - sdf.get_distance(p0 + float3(0, -h, 0))) / (2 * h);
              float ddz = (sdf.get_distance(p0 + float3(0, 0, h)) - sdf.get_distance(p0 + float3(0, 0, -h))) / (2 * h);
              float3 n = normalize(float3(ddx, ddy, ddz));
              float shadow = primitive_sdf_sphere_tracing(sdf, sdf_bbox, p0 + 10.0f * h * light_dir, light_dir);
              if (soft_shadows)
              {
                unsigned steps = 64;
                float res = 1.0f;
                float t = h;
                for (int i=0;i<steps;i++)    
                {
                  float3 p = p0 + t*light_dir;
                  if (!sdf_bbox.contains(p))
                    break;
                  float d = sdf.get_distance(p);
                  if (d < 1e-6)
                  {
                    res = 0.0f;
                    break;
                  }
                  res = fminf(res, 10*d/t);
                  t += d;
                }   
                shadow = 1.0f - fminf(res, 1.0f);     
              }
              
              float l_val = (1 - shadow) * dot(n, light_dir);
              
              float ao_value = 0.0f;
              if (ao)
              {
                unsigned rays = 16;
                unsigned hits = 0;
                float step_size = 0.01;
                /* slow but correct ao
                float3 t1 = cross(n, abs(l_val) > 0.99 ? float3(-light_dir.y, light_dir.x, light_dir.z) : float3(light_dir.x, light_dir.y, light_dir.z));
                float3 t2 = cross(n, t1);

                for (int i = 0; i < rays; i++)
                {
                  float phi = urand(-PI, PI);
                  float psi = urand(0, PI / 2);
                  float3 ray = cos(phi) * cos(psi) * t1 + sin(psi) * n + sin(phi) * cos(psi) * t2;
                  hits += primitive_sdf_sphere_tracing(sdf, sdf_bbox, p0 + h * n, ray);
                }
                ao_value = hits / (float)rays;
                */
                //Distance Field Ambient Occlusion https://zephyrl.github.io/SDF/
                for (int i = 1; i <= rays; i++)
                {
                  float dist = step_size*i;
                  float i_intensity = 0.33f/rays;
                  ao_value += fmaxf((dist - sdf.get_distance(p0 + dist*n))/dist, 0.0f)*i_intensity;
                }
                ao_value = std::min(1.0f, ao_value);
              }
              float l = MAX(0.1f, l_val) * (1.0f - ao_value);
              switch (color_mode)
              {
              case GRAY:
                color += l*float3(1, 1, 1);
                break;
              case SINPOS:
                color += l*LiteMath::abs(float3(sin(8*PI*p0.x), sin(8*PI*p0.y), sin(8*PI*p0.z)));
                break;
              case NORMALS:
                color += l*LiteMath::abs(n);
                break;              
              default:
                break;
              }
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
    cam.origin = float3(0,3,3);
    cam.target = float3(0,0,0);
    cam.up = float3(0,1,0);

    float3 light_dir = normalize(float3(0.7,0.8,0.5));
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};
    l.to_file("saves/task1_references/light.txt");

    CameraSettings cam1 = cam;
    cam1.origin = float3(2*sin(0),2,2*cos(0));
    convert(cam1).to_file("saves/task1_references/cam1.txt");

    CameraSettings cam2 = cam;
    cam2.origin = float3(2*sin(2*PI/3),2,2*cos(2*PI/3));
    convert(cam2).to_file("saves/task1_references/cam2.txt");

    CameraSettings cam3 = cam;
    cam3.origin = float3(2*sin(4*PI/3),2,2*cos(4*PI/3));
    convert(cam3).to_file("saves/task1_references/cam3.txt");

    PrimitiveSDF sdf1;
    sdf1.primitives.push_back(Primitive(MANDELBULB, {0,0,0},{0.7,0.7,0.7}));

    {
      PrimitiveSDFScene psdf;
      for (auto &p : sdf1.primitives)
      {
        psdf.scene_data.push_back(p.type);
        psdf.scene_data.push_back(p.shift.x);
        psdf.scene_data.push_back(p.shift.y);
        psdf.scene_data.push_back(p.shift.z);
        psdf.scene_data.push_back(p.scale.x);
        psdf.scene_data.push_back(p.scale.y);
        psdf.scene_data.push_back(p.scale.z);
      }
      psdf.to_file("saves/task1_references/scene1.txt");
    }
    {
      PrimitiveSDFScene psdf;
      psdf.from_file("saves/task1_references/scene1.txt");
      sdf1.primitives = std::vector<Primitive>(psdf.scene_data.size()/7);
      for (int i=0;i<sdf1.primitives.size();i++)
        sdf1.primitives[i] = Primitive((PrimitiveType)((int)psdf.scene_data[7*i+0]),
                                       float3(psdf.scene_data[7*i+1], psdf.scene_data[7*i+2], psdf.scene_data[7*i+3]),
                                       float3(psdf.scene_data[7*i+4], psdf.scene_data[7*i+5], psdf.scene_data[7*i+6]));
    }
    
    PrimitiveSDF sdf2;
    sdf2.primitives.push_back(Primitive(SPONGE, {0,0.6,0},{0.6,0.6,0.6}));
    sdf2.primitives.push_back(Primitive(BOX, {0,-5,0.0f},{5,5,5}));

    {
      PrimitiveSDFScene psdf;
      for (auto &p : sdf2.primitives)
      {
        psdf.scene_data.push_back(p.type);
        psdf.scene_data.push_back(p.shift.x);
        psdf.scene_data.push_back(p.shift.y);
        psdf.scene_data.push_back(p.shift.z);
        psdf.scene_data.push_back(p.scale.x);
        psdf.scene_data.push_back(p.scale.y);
        psdf.scene_data.push_back(p.scale.z);
      }
      psdf.to_file("saves/task1_references/scene2.txt");
    }
    {
      PrimitiveSDFScene psdf;
      psdf.from_file("saves/task1_references/scene2.txt");
      sdf2.primitives = std::vector<Primitive>(psdf.scene_data.size()/7);
      for (int i=0;i<sdf2.primitives.size();i++)
        sdf2.primitives[i] = Primitive((PrimitiveType)((int)psdf.scene_data[7*i+0]),
                                       float3(psdf.scene_data[7*i+1], psdf.scene_data[7*i+2], psdf.scene_data[7*i+3]),
                                       float3(psdf.scene_data[7*i+4], psdf.scene_data[7*i+5], psdf.scene_data[7*i+6]));
    }

    PrimitiveSDF sdf3;
    sdf3.primitives.push_back(Primitive(SPHERE, {0.9,0.3,0.9},{0.1,0.1,0.1}));
    sdf3.primitives.push_back(Primitive(SPHERE, {0.76,0.3,0.35},{0.1,0.1,0.1}));
    sdf3.primitives.push_back(Primitive(SPHERE, {0.54,0.3,0.19},{0.1,0.1,0.1}));
    sdf3.primitives.push_back(Primitive(SPHERE, {0.32,0.3,0.57},{0.1,0.1,0.1}));
    sdf3.primitives.push_back(Primitive(SPHERE, {0.5,0.35,0.79},{0.15,0.15,0.15}));
    sdf3.primitives.push_back(Primitive(BOX, {0,-0.8,0},{1,1,1}));
    sdf3.primitives.push_back(Primitive(BOX, {0,-5,0.0f},{5,5,5}));

    Texture t;

    //int steps = 16;
    //for (int i=0;i<steps;i++)
    //{
    //  CameraSettings tc = cam;
    //  tc.origin = float3(3*cos(2*PI*i/steps), 2, 3*sin(2*PI*i/steps));
    //  t = render_primitive_sdf(sdf2, tc, light_dir, 1024, 1024, 16, GRAY, true, true);
    //  engine::textureManager->save_png(t, "task1_references/sdf1_turntable_"+std::to_string(i));
    //}

    t = render_primitive_sdf(sdf1, cam1, light_dir, 2048, 2048, 32, NORMALS);
    engine::textureManager->save_png(t, "task1_references/ref_1");
    t = render_primitive_sdf(sdf1, cam1, light_dir, 2048, 2048, 32, GRAY);
    engine::textureManager->save_png(t, "task1_references/ref_2");
    t = render_primitive_sdf(sdf1, cam1, light_dir, 2048, 2048, 32, SINPOS);
    engine::textureManager->save_png(t, "task1_references/ref_3");

    t = render_primitive_sdf(sdf2, cam1, light_dir, 2048, 2048, 32, GRAY);
    engine::textureManager->save_png(t, "task1_references/ref_4");
    t = render_primitive_sdf(sdf2, cam1, light_dir, 2048, 2048, 32, GRAY, true);
    engine::textureManager->save_png(t, "task1_references/ref_5");
    t = render_primitive_sdf(sdf2, cam1, light_dir, 2048, 2048, 32, GRAY, true, true);
    engine::textureManager->save_png(t, "task1_references/ref_6");

    t = render_primitive_sdf(sdf3, cam1, light_dir, 2048, 2048, 32, GRAY);
    engine::textureManager->save_png(t, "task1_references/ref_7");
    t = render_primitive_sdf(sdf3, cam1, light_dir, 2048, 2048, 32, GRAY, true);
    engine::textureManager->save_png(t, "task1_references/ref_8");
    t = render_primitive_sdf(sdf3, cam1, light_dir, 2048, 2048, 32, GRAY, true, true);
    engine::textureManager->save_png(t, "task1_references/ref_9");
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

  static float3 closest_point_triangle(const float3& p, const float3& a, const float3& b, const float3& c)
  {
    //implementation taken from Embree library
    const float3 ab = b - a;
    const float3 ac = c - a;
    const float3 ap = p - a;

    const float d1 = dot(ab, ap);
    const float d2 = dot(ac, ap);
    if (d1 <= 0.f && d2 <= 0.f) return a; //#1

    const float3 bp = p - b;
    const float d3 = dot(ab, bp);
    const float d4 = dot(ac, bp);
    if (d3 >= 0.f && d4 <= d3) return b; //#2

    const float3 cp = p - c;
    const float d5 = dot(ab, cp);
    const float d6 = dot(ac, cp);
    if (d6 >= 0.f && d5 <= d6) return c; //#3

    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
    {
        const float v = d1 / (d1 - d3);
        return a + v * ab; //#4
    }
      
    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
    {
        const float v = d2 / (d2 - d6);
        return a + v * ac; //#5
    }
      
    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
    {
        const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + v * (c - b); //#6
    }

    const float denom = 1.f / (va + vb + vc);
    const float v = vb * denom;
    const float w = vc * denom;
    return a + v * ab + w * ac; //#0
  }

  void model_to_point_cloud(const Mesh *m, int count, AABB bbox, std::vector<float> *points, std::vector<float> *distances)
  {
    std::vector<float3> face_normals(m->indices.size()/3);
    for (int i=0;i<m->indices.size();i+=3)
    {
      float3 p0(m->positions[3*m->indices[i+0]+0], m->positions[3*m->indices[i+0]+1], m->positions[3*m->indices[i+0]+2]);
      float3 p1(m->positions[3*m->indices[i+1]+0], m->positions[3*m->indices[i+1]+1], m->positions[3*m->indices[i+1]+2]);
      float3 p2(m->positions[3*m->indices[i+2]+0], m->positions[3*m->indices[i+2]+1], m->positions[3*m->indices[i+2]+2]);
      face_normals[i/3] = normalize(cross(p1-p0, p2-p0));
    }

    //BVH bvh(m->positions, m->indices);

    points->resize(3*count);
    distances->resize(count);
    std::vector<float3> surface_points(count);
    for (int j=0;j<count;j++)
    {
      float3 p = float3(urand(bbox.min_pos.x, bbox.max_pos.x),
                              urand(bbox.min_pos.y, bbox.max_pos.y),
                              urand(bbox.min_pos.z, bbox.max_pos.z));
      (*points)[3*j+0] = p.x;
      (*points)[3*j+1] = p.y;
      (*points)[3*j+2] = p.z;

      float min_dist = 1e9;
      float3 closest_pos(0,0,0);
      float3 closest_norm(1,0,0);
      for (int i=0;i<m->indices.size();i+=3)
      {
        float3 p0(m->positions[3*m->indices[i+0]+0], m->positions[3*m->indices[i+0]+1], m->positions[3*m->indices[i+0]+2]);
        float3 p1(m->positions[3*m->indices[i+1]+0], m->positions[3*m->indices[i+1]+1], m->positions[3*m->indices[i+1]+2]);
        float3 p2(m->positions[3*m->indices[i+2]+0], m->positions[3*m->indices[i+2]+1], m->positions[3*m->indices[i+2]+2]);
        
        float3 tp = closest_point_triangle(p, p0, p1, p2);
        float d = length(p-tp);
        if (d < min_dist)
        {
          closest_pos = tp;
          closest_norm = face_normals[i/3];
          min_dist = d;
        }
      }
      if (dot(p-closest_pos, closest_norm) < 0)
        min_dist = -min_dist;

      (*distances)[j] = min_dist;
      surface_points[j] = closest_pos;
      //printf("%f %f %f d= %f\n", p.x, p.y, p.z, min_dist);
    }
    /*
    for (int i=0;i<10;i++)
    {
      CameraSettings cam;
      cam.origin = float3(3*sin(0.2*PI*i),0, 3*cos(0.2*PI*i));
      cam.target = float3(0,0,0);
      cam.up = float3(0,1,0);
      PointCloudRenderer pcr;
      Texture t = pcr.render(surface_points, cam.get_viewProj(), 1024, 1024);
      engine::textureManager->save_png(t, "task2_references/point_cloud_"+std::to_string(i));
    }
    */
  }

  void create_point_cloud_house(int count, AABB bbox, std::vector<float> *points, std::vector<float> *distances)
  {
    Model *m = model_loader::load_model_from_obj_directly("resources/models/cup_1n.obj");
    model_loader::normalize_model(m);
    model_to_point_cloud(m, count, bbox, points, distances);
    delete m;
  }

  void create_point_cloud_spheres(int count, AABB bbox, std::vector<float> *points, std::vector<float> *distances)
  {
    PrimitiveSDF sdf2;
    sdf2.primitives.push_back(Primitive(SPHERE, {0,0.3f,0.4f},{0.4,0.4,0.4}));
    sdf2.primitives.push_back(Primitive(SPHERE, {-0.3f,0,0.1f},{0.3,0.3,0.3}));
    sdf2.primitives.push_back(Primitive(SPHERE, {0,-0.2f,0.0f},{0.25,0.25,0.25}));

    points->resize(3*count);
    distances->resize(count);

    for (int i=0;i<count;i++)
    {
      float3 p = float3(urand(bbox.min_pos.x, bbox.max_pos.x),
                              urand(bbox.min_pos.y, bbox.max_pos.y),
                              urand(bbox.min_pos.z, bbox.max_pos.z));
      (*points)[3*i+0] = p.x;
      (*points)[3*i+1] = p.y;
      (*points)[3*i+2] = p.z;
      (*distances)[i] = sdf2.get_distance(p);
    }
  }

  void create_real_params_from_chair_depended(std::vector<float> chair, std::vector<float> *params)
  {
    params->push_back(0);//first box to sit
    params->push_back(-chair[3]);
    params->push_back(0);
    params->push_back(chair[4]);
    params->push_back(chair[3]);
    params->push_back(chair[2]);

    params->push_back(chair[3] - chair[4]);//second another box
    params->push_back(-chair[5]);
    params->push_back(0);
    params->push_back(chair[3]);
    params->push_back(chair[5]);
    params->push_back(chair[2]);

    for (int j = 0; j < 4; ++j)//legs
    {
      float mult1 = 1, mult2 = 1;
      if (j % 2 != 0) mult1 = -1;
      if (j / 2 != 0) mult2 = -1;
      params->push_back((chair[4] - chair[0]) * mult1);
      params->push_back(chair[1]);
      params->push_back((chair[2] - chair[0]) * mult2);
      params->push_back(chair[1]);
      params->push_back(chair[0]);
    }
  }

  void create_point_and_params_cloud_chair(int count_params, int count_points, AABB bbox, std::vector<float> *points, std::vector<float> *distances)
  {
    upg::ProceduralSdf sdf({std::vector<uint16_t>{3, 3, 2, 4, 2, 4, 3, 3, 2, 5, 2, 5, 3, 2, 5, 2, 5}});
    unsigned int sz = 3 + 6;
    points->resize(sz*count_params*count_points);
    distances->resize(count_params*count_points);
    for (int i=0;i<count_params;i++)
    {
      std::vector<float> params = {}, par = {};//3mov 3box 3mov 3box 3mov 2cyl 3mov 2cyl 3mov 2cyl 3mov 2cyl
      params.push_back(urand(0.001, 0.4));//radius
      params.push_back(urand(0.001, 0.5));//height1
      params.push_back(urand(0.001, 1));//width
      params.push_back(urand(0.001, 0.3));//thick
      params.push_back(urand(0.001, 1));//length
      params.push_back(urand(0.001, 0.5));//height2

      create_real_params_from_chair_depended(params, &par);
      sdf.set_parameters(par);
      for (int j = 0; j < count_points; ++j)
      {
        float3 p = float3(urand(bbox.min_pos.x, bbox.max_pos.x),
                                urand(bbox.min_pos.y, bbox.max_pos.y),
                                urand(bbox.min_pos.z, bbox.max_pos.z));
        (*points)[sz*i*count_points+sz*j+0] = p.x;
        (*points)[sz*i*count_points+sz*j+1] = p.y;
        (*points)[sz*i*count_points+sz*j+2] = p.z;
        for (int k = 0; k < params.size(); ++k) (*points)[sz*i*count_points+sz*j+3+k] = params[k];
        (*distances)[i*count_points+j] = sdf.get_distance(p);
      }
    }
  }

  void create_point_and_params_cloud_32124(int count_params, int count_points, AABB bbox, std::vector<float> *points, std::vector<float> *distances)
  {
    upg::ProceduralSdf sdf({std::vector<uint16_t>{3, 2, 1, 2, 4}});
    unsigned int sz = 3 + sdf.desc.get_total_params_count();
    points->resize(sz*count_params*count_points);
    distances->resize(count_params*count_points);
    for (int i=0;i<count_params;i++)
    {
      std::vector<float> params = {};
      for (auto &par : sdf.desc.get_block_params())
      {
        for (auto &pr : par.second.p)
        {
          params.push_back(urand(pr.min_val, pr.max_val));
          //debug("%f, %f\n", pr.min_val, pr.max_val);
        }
      }
      //debug("\n");
      for (int j = 0; j < count_points; ++j)
      {
        float3 p = float3(urand(bbox.min_pos.x, bbox.max_pos.x),
                                urand(bbox.min_pos.y, bbox.max_pos.y),
                                urand(bbox.min_pos.z, bbox.max_pos.z));
        sdf.set_parameters(params);
        (*points)[sz*i*count_points+sz*j+0] = p.x;
        (*points)[sz*i*count_points+sz*j+1] = p.y;
        (*points)[sz*i*count_points+sz*j+2] = p.z;
        for (int k = 0; k < params.size(); ++k) (*points)[sz*i*count_points+sz*j+3+k] = params[k];
        (*distances)[i*count_points+j] = sdf.get_distance(p);
      }
    }
  }

  std::vector<float> add_params_to_points(std::vector<float> points, std::vector<float> params)
  {
    std::vector<float> res(0);
    for (int i = 0; i < points.size() / 3; ++i)
    {
      res.push_back(points[i * 3 + 0]);
      res.push_back(points[i * 3 + 1]);
      res.push_back(points[i * 3 + 2]);
      res.insert(res.end(), params.begin(), params.end());
    }
    return res;
  }

  Texture render_neural_sdf_with_params(nn::Siren &sdf, std::vector<float> params, AABB bbox, const CameraSettings &camera, 
                            int image_w, int image_h, int spp, bool lambert, float3 light_dir)
  {
    float3 center = 0.5f*(bbox.max_pos + bbox.min_pos);
    float3 size = 0.5f*(bbox.max_pos - bbox.min_pos);
    AABB inflated_bbox = AABB(center - 1.1f*size, center + 1.1f*size);

    float4x4 projInv = LiteMath::inverse4x4(camera.get_proj());
    float4x4 viewInv = LiteMath::inverse4x4(camera.get_view());

    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    std::vector<int> hits(image_w*image_h, 0); //one per pixel
    std::vector<float3> colors(image_w*image_h, float3(0,0,0));
    std::vector<int> pixel_indices(image_w*image_h*spp_a*spp_a);
    std::vector<float3> dirs(image_w*image_h*spp_a*spp_a);
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
            float3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            float3 p0 = camera.origin;
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
      std::vector<float> tmp = add_params_to_points(points, params);
      sdf.evaluate(tmp, distances, points_left);

      int new_points_left = points_left;
      int free_pos = 0;
      for (int i=0;i<points_left;i++)
      {
        float d = distances[i];
        float3 p0 = float3(points[3*i+0], points[3*i+1], points[3*i+2]);
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
        colors[i] += float3(hits[i], hits[i], hits[i]);
    }
    else
    {
      std::vector<float> tmp = add_params_to_points(fd_points, params);
      sdf.evaluate(tmp, fd_distances, 6*points_found);
      for (int i=0;i<points_found;i++)
      {
        float ddx = (fd_distances[6*i+0] - fd_distances[6*i+1])/(2*h);
        float ddy = (fd_distances[6*i+2] - fd_distances[6*i+3])/(2*h);
        float ddz = (fd_distances[6*i+4] - fd_distances[6*i+5])/(2*h);
        float3 n = normalize(float3(ddx, ddy, ddz));
        colors[fd_pixel_indices[i]] += float3(1,1,1) * MAX(0.1f, dot(n, light_dir));
      }
    }

    for (int i = 0; i < image_w * image_h; i++)
    {
      float3 color = colors[i];
      data[4 * i + 0] = 255 * (color.x / SQR(spp_a));
      data[4 * i + 1] = 255 * (color.y / SQR(spp_a));
      data[4 * i + 2] = 255 * (color.z / SQR(spp_a));
      data[4 * i + 3] = 255;
    }

    Texture t = engine::textureManager->create_texture(image_w, image_h, GL_RGBA8, 1, data, GL_RGBA);
    delete[] data;
    return t;
  }

  Texture render_neural_sdf(nn::Siren &sdf, AABB bbox, const CameraSettings &camera, 
                            int image_w, int image_h, int spp, bool lambert, float3 light_dir)
  {
    sdf.set_batch_size_for_evaluate(4096);

    float3 center = 0.5f*(bbox.max_pos + bbox.min_pos);
    float3 size = 0.5f*(bbox.max_pos - bbox.min_pos);
    AABB inflated_bbox = AABB(center - 1.1f*size, center + 1.1f*size);

    float4x4 projInv = LiteMath::inverse4x4(camera.get_proj());
    float4x4 viewInv = LiteMath::inverse4x4(camera.get_view());

    int spp_a = MAX(1,floor(sqrtf(spp)));
    unsigned char *data = new unsigned char[4*image_w*image_h];

    std::vector<int> hits(image_w*image_h, 0); //one per pixel
    std::vector<float3> colors(image_w*image_h, float3(0,0,0));
    std::vector<int> pixel_indices(image_w*image_h*spp_a*spp_a);
    std::vector<float3> dirs(image_w*image_h*spp_a*spp_a);
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
            float3 dir = transformRay(EyeRayDirNormalized(x,y,projInv), viewInv);
            float3 p0 = camera.origin;
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
        float3 p0 = float3(points[3*i+0], points[3*i+1], points[3*i+2]);
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
      //logerr("iteration %d, points left %d/%d", iteration, points_left, index);
    }

    if (!lambert)
    {
      for (int i=0;i<image_w*image_h;i++)
        colors[i] += float3(hits[i], hits[i], hits[i]);
    }
    else
    {
      sdf.evaluate(fd_points, fd_distances, 6*points_found);
      for (int i=0;i<points_found;i++)
      {
        float ddx = (fd_distances[6*i+0] - fd_distances[6*i+1])/(2*h);
        float ddy = (fd_distances[6*i+2] - fd_distances[6*i+3])/(2*h);
        float ddz = (fd_distances[6*i+4] - fd_distances[6*i+5])/(2*h);
        float3 n = normalize(float3(ddx, ddy, ddz));
        colors[fd_pixel_indices[i]] += float3(1,1,1) * MAX(0.1f, dot(n, light_dir));
      }
    }

    for (int i = 0; i < image_w * image_h; i++)
    {
      float3 color = colors[i];
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
    cam.origin = float3(0,0,3);
    cam.target = float3(0,0,0);
    cam.up = float3(0,1,0);

    float3 light_dir = normalize(cam.origin + float3(cam.origin.z, cam.origin.y, cam.origin.x) - cam.target);
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};
    l.to_file("saves/task2_references/light.txt");

    CameraSettings cam1 = cam;
    cam1.origin = float3(0,3*sin(0),3*cos(0));
    convert(cam1).to_file("saves/task2_references/cam1.txt");

    CameraSettings cam2 = cam;
    cam2.origin = float3(0,3*sin(2*PI/3),3*cos(2*PI/3));
    convert(cam2).to_file("saves/task2_references/cam2.txt");

    CameraSettings cam3 = cam;
    cam3.origin = float3(0,3*sin(4*PI/3),3*cos(4*PI/3));
    convert(cam3).to_file("saves/task2_references/cam3.txt");

    {
      std::vector<float> points, distances;
      create_point_cloud_spheres(5000, bbox, &points, &distances);
      save_points_cloud("saves/task2_references/sdf1_points.bin", points, distances);
      nn::Siren network(nn::Siren::Type::SDF, 3, 64);
      network.train(points, distances, 512, 15000, true);
      network.save_weights_to_file("saves/task2_references/sdf1_weights.bin");
      network.set_arch_to_file("saves/task2_references/sdf1_arch.txt");
      
      std::vector<float> estimated_distances = distances;
      network.evaluate(points, estimated_distances);
      save_points_cloud("saves/task2_references/sdf1_test.bin", points, estimated_distances);
    }
    {
      std::vector<float> points, distances;
      load_points_cloud("saves/task2_references/sdf1_points.bin", &points, &distances);
      nn::Siren network(nn::Siren::Type::SDF, 3, 64);
      network.initialize_from_file("saves/task2_references/sdf1_weights.bin");
      
      Texture t;
      t = render_neural_sdf(network, bbox, cam1, 64, 64, 1, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf1_cam1_reference");
      t = render_neural_sdf(network, bbox, cam2, 64, 64, 1, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf1_cam2_reference");
      t = render_neural_sdf(network, bbox, cam3, 64, 64, 1, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf1_cam3_reference");
    }

    {
      std::vector<float> points, distances;
      create_point_cloud_house(50000, bbox, &points, &distances);
      save_points_cloud("saves/task2_references/sdf2_points.bin", points, distances);
      nn::Siren network(nn::Siren::Type::SDF, 5, 256);
      network.train(points, distances, 512, 15000, true);
      network.save_weights_to_file("saves/task2_references/sdf2_weights.bin");
      network.set_arch_to_file("saves/task2_references/sdf2_arch.txt");
            
      std::vector<float> estimated_distances = distances;
      network.evaluate(points, estimated_distances);
      save_points_cloud("saves/task2_references/sdf2_test.bin", points, estimated_distances);
    }
    {
      std::vector<float> points, distances;
      load_points_cloud("saves/task2_references/sdf2_points.bin", &points, &distances);
      nn::Siren network(nn::Siren::Type::SDF, 5, 256);
      network.initialize_from_file("saves/task2_references/sdf2_weights.bin");
      
      Texture t;
      t = render_neural_sdf(network, bbox, cam1, 512, 512, 9, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf2_cam1_reference");
      t = render_neural_sdf(network, bbox, cam2, 512, 512, 9, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf2_cam2_reference");
      t = render_neural_sdf(network, bbox, cam3, 512, 512, 9, true, light_dir);
      engine::textureManager->save_png(t, "task2_references/sdf2_cam3_reference");
    }
  }

  void task_3_create_references()
  {
    AABB bbox({-1,-1,-1},{1,1,1});

    CameraSettings cam;
    cam.origin = float3(0,0,3);
    cam.target = float3(0,0,0);
    cam.up = float3(0,1,0);

    float3 light_dir = normalize(cam.origin + float3(cam.origin.z, cam.origin.y, cam.origin.x) - cam.target);
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};
    l.to_file("saves/task3_references/light.txt");

    CameraSettings cam1 = cam;
    cam1.origin = float3(0,3*sin(0),3*cos(0));
    convert(cam1).to_file("saves/task3_references/cam1.txt");

    CameraSettings cam2 = cam;
    cam2.origin = float3(0,3*sin(2*PI/3),3*cos(2*PI/3));
    convert(cam2).to_file("saves/task3_references/cam2.txt");

    CameraSettings cam3 = cam;
    cam3.origin = float3(0,3*sin(4*PI/3),3*cos(4*PI/3));
    convert(cam3).to_file("saves/task3_references/cam3.txt");

    {
      std::vector<float> points_and_params, distances;
      //create_point_and_params_cloud_1(1, 100000000, bbox, &points_and_params, &distances);
      //save_points_cloud("saves/task2_references/sdf1_points.bin", points, distances);
      nn::Siren network(nn::Siren::Type::Gen_SDF_32124, 5, 256);
      debug("START TRAIN\n");
      create_point_and_params_cloud_32124(100000, 100, bbox, &points_and_params, &distances);
      network.train(points_and_params, distances, 512, 80000, true);
      debug("END TRAIN\n");
      network.save_weights_to_file("saves/task3_references/sdf_and_params_weights.bin");
      network.set_arch_to_file("saves/task3_references/sdf_and_params_arch.txt");
      debug("END SAVE WEIGHTS\n");
      //std::vector<float> estimated_distances = distances;
      //network.evaluate(points_and_params, estimated_distances);
      //save_points_cloud("saves/task2_references/sdf1_test.bin", points, estimated_distances);
      //debug("END EVALUATE\n");
      nn::Siren network2(nn::Siren::Type::Gen_SDF_32124, 5, 256);
      network2.initialize_from_file("saves/task3_references/sdf_and_params_weights.bin");
      //std::vector<float> p = {0, 0, 0, 0.2, 0, 0, -0.2, 0.2};
      //std::vector<float> out = {0.2, 0};
      //std::vector<float> out2 = {0.2, 0};
      //network2.evaluate(p, out2);
      //debug("OUTS %f, %f - %f, %f", out[0], out2[0], out[1], out2[1]);
      debug("START RENDER 1\n");
      Texture t;
      t = render_neural_sdf_with_params(network2, {0.5, 0.0, 0.0, 0.2, 0.0, -0.5, -0.1, 0.3, 0.1, 0.4}, bbox, cam1, 256, 256, 1, true, light_dir);
      engine::textureManager->save_png(t, "task3_references/sdf_and_params_cam1_reference");
      debug("START RENDER 2\n");
      t = render_neural_sdf_with_params(network2, {0.5, 0.0, 0.0, 0.2, 0.0, -0.5, -0.1, 0.3, 0.1, 0.4}, bbox, cam2, 256, 256, 1, true, light_dir);
      engine::textureManager->save_png(t, "task3_references/sdf_and_params_cam2_reference");
      debug("START RENDER 3\n");
      t = render_neural_sdf_with_params(network2, {0.5, 0.0, 0.0, 0.2, 0.0, -0.5, -0.1, 0.3, 0.1, 0.4}, bbox, cam3, 256, 256, 1, true, light_dir);
      engine::textureManager->save_png(t, "task3_references/sdf_and_params_cam3_reference");
      debug("END RENDER\n");
    }
  }

  void task_4_chair()
  {
    AABB bbox({-1,-1,-1},{1,1,1});

    CameraSettings cam;
    cam.origin = float3(0,0,3);
    cam.target = float3(0,0,0);
    cam.up = float3(0,1,0);

    float3 light_dir = normalize(cam.origin + float3(cam.origin.z, cam.origin.y, cam.origin.x) - cam.target);
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};
    l.to_file("saves/task3_references/light.txt");

    CameraSettings cam1 = cam;
    cam1.origin = float3(3*sin(0),0,3*cos(0));
    convert(cam1).to_file("saves/task3_references/cam1.txt");

    CameraSettings cam2 = cam;
    cam2.origin = float3(3*sin(2*PI/3),0,3*cos(2*PI/3));
    convert(cam2).to_file("saves/task3_references/cam2.txt");

    CameraSettings cam3 = cam;
    cam3.origin = float3(3*sin(4*PI/3),0,3*cos(4*PI/3));
    convert(cam3).to_file("saves/task3_references/cam3.txt");

    {
      std::vector<float> points_and_params, distances;
      //create_point_and_params_cloud_1(1, 100000000, bbox, &points_and_params, &distances);
      //save_points_cloud("saves/task2_references/sdf1_points.bin", points, distances);
      nn::Siren network(nn::Siren::Type::Chair, 5, 256);
      debug("START TRAIN\n");
      create_point_and_params_cloud_chair(100000, 100, bbox, &points_and_params, &distances);
      network.train(points_and_params, distances, 512, 80000, true);
      debug("END TRAIN\n");
      network.save_weights_to_file("saves/task3_references/sdf_and_params_weights.bin");
      network.set_arch_to_file("saves/task3_references/sdf_and_params_arch.txt");
      debug("END SAVE WEIGHTS\n");
      //std::vector<float> estimated_distances = distances;
      //network.evaluate(points_and_params, estimated_distances);
      //save_points_cloud("saves/task2_references/sdf1_test.bin", points, estimated_distances);
      //debug("END EVALUATE\n");
      nn::Siren network2(nn::Siren::Type::Chair, 5, 256);
      network2.initialize_from_file("saves/task3_references/sdf_and_params_weights.bin");
      //std::vector<float> p = {0, 0, 0, 0.2, 0, 0, -0.2, 0.2};
      //std::vector<float> out = {0.2, 0};
      //std::vector<float> out2 = {0.2, 0};
      //network2.evaluate(p, out2);
      //debug("OUTS %f, %f - %f, %f", out[0], out2[0], out[1], out2[1]);
      debug("START RENDER 1\n");
      Texture t;
      t = render_neural_sdf_with_params(network2, {0.2, 0.4, 0.5, 0.05, 0.7, 0.4}, bbox, cam1, 256, 256, 1, true, light_dir);
      engine::textureManager->save_png(t, "task3_references/sdf_and_params_cam1_reference");
      debug("START RENDER 2\n");
      t = render_neural_sdf_with_params(network2, {0.2, 0.4, 0.5, 0.05, 0.7, 0.4}, bbox, cam2, 256, 256, 1, true, light_dir);
      engine::textureManager->save_png(t, "task3_references/sdf_and_params_cam2_reference");
      debug("START RENDER 3\n");
      t = render_neural_sdf_with_params(network2, {0.2, 0.4, 0.5, 0.05, 0.7, 0.4}, bbox, cam3, 256, 256, 1, true, light_dir);
      engine::textureManager->save_png(t, "task3_references/sdf_and_params_cam3_reference");
      debug("END RENDER\n");
    }
  }

  void test_complex_model()
  {
    AABB bbox({-1,-1,-1},{1,1,1});

    CameraSettings cam;
    cam.origin = float3(0,0,3);
    cam.target = float3(0,0,0);
    cam.up = float3(0,1,0);
    CameraSettings cam1 = cam; cam1.origin = float3(0,3*sin(0),3*cos(0));
    CameraSettings cam2 = cam; cam2.origin = float3(0,3*sin(2*PI/3),3*cos(2*PI/3));
    CameraSettings cam3 = cam; cam3.origin = float3(0,3*sin(4*PI/3),3*cos(4*PI/3));

    float3 light_dir = normalize(cam.origin + float3(cam.origin.z, cam.origin.y, cam.origin.x) - cam.target);
    DirectedLight l{light_dir.x, light_dir.y, light_dir.z, 1.0f};

    std::vector<float> points, distances;
    Model *m = model_loader::load_model_from_obj_directly("/home/sammael/grade_resources/Huawei_models/h1_norm.obj");
    model_loader::normalize_model(m);
    model_to_point_cloud(m, 100000, bbox, &points, &distances);
    delete m;

    nn::Siren network(nn::Siren::Type::SDF, 5, 512);
    network.train(points, distances, 1024, 5000, true);
    Texture t;
    t = render_neural_sdf(network, bbox, cam1, 512, 512, 1, true, light_dir);
    engine::textureManager->save_png(t, "detail_nSDF_cam1");
    t = render_neural_sdf(network, bbox, cam2, 512, 512, 1, true, light_dir);
    engine::textureManager->save_png(t, "detail_nSDF_cam2");
    t = render_neural_sdf(network, bbox, cam3, 512, 512, 1, true, light_dir);
    engine::textureManager->save_png(t, "detail_nSDF_cam3");
  }

  void neural_SDF_test()
  {
    nn::TensorProcessor::init("GPU");
    task_1_create_references();
    //task_2_create_references();
    //task_3_create_references();
    //task_4_chair();
  }
}