#include "sdf_scene.h"
#include <cassert>
#include <fstream>

using namespace LiteMath;

  float dist_prim(const SdfScene &sdf, const SdfObject &prim, float3 p)
  {
    float3 pos = prim.transform * p;
    //printf("%f %f %f -- %f %f %f\n", p.x, p.y, p.z, pos.x, pos.y, pos.z);
    //printf("%f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n", 
    //       prim.transform(0, 0), prim.transform(0, 1), prim.transform(0, 2), prim.transform(0, 3),
    //       prim.transform(1, 0), prim.transform(1, 1), prim.transform(1, 2), prim.transform(1, 3),
    //       prim.transform(2, 0), prim.transform(2, 1), prim.transform(2, 2), prim.transform(2, 3),
    //       prim.transform(3, 0), prim.transform(3, 1), prim.transform(3, 2), prim.transform(3, 3));

    switch (prim.type)
    {
    case 1:
    {
      float r = sdf.parameters[prim.params_offset + 0];
      //fprintf(stderr, "sphere %f %f %f - %f",pos.x, pos.y, pos.z, r);
      return length(pos) - r;
    }
    case 4:
    {
      float3 size(sdf.parameters[prim.params_offset + 0],
                     sdf.parameters[prim.params_offset + 1],
                     sdf.parameters[prim.params_offset + 2]);
      //fprintf(stderr, "box %f %f %f - %f %f %f - %f %f %f",p.x, p.y, p.z, pos.x, pos.y, pos.z, size.x, size.y, size.z);
      float3 q = abs(pos) - size;
      return length(max(q,float3(0.0f))) + min(max(q.x,max(q.y,q.z)),0.0f);
    }    
    case 5:
    {
      float h = sdf.parameters[prim.params_offset + 0];
      float r = sdf.parameters[prim.params_offset + 1];
      float2 d = abs(float2(sqrt(pos.x*pos.x + pos.z*pos.z),pos.y)) - float2(r,h);
      return min(max(d.x,d.y),0.0f) + length(max(d,float2(0.0f)));
    }
    default:
      fprintf(stderr, "unknown type %u", prim.type);
      assert(false);
      break;
    }
    return -1000;
  }

  float get_dist(const SdfScene &sdf, float3 p)
  {
    float d = 1e6;
    for (auto &conj : sdf.conjunctions)
    {
      float conj_d = -1e6;
      for (unsigned pid = conj.offset; pid < conj.offset + conj.size; pid++)
      {
        float prim_d = sdf.objects[pid].distance_mult*dist_prim(sdf, sdf.objects[pid], p) + sdf.objects[pid].distance_add;
        //fprintf(stderr, "prim %f", prim_d);
        conj_d = max(conj_d, sdf.objects[pid].complement ? -prim_d : prim_d);
      }
      //fprintf(stderr, "%f %f\n", d, conj_d);
      d = min(d, conj_d);
    }
    //fprintf(stderr, "\n");
    return d;
  }

  bool sdf_sphere_tracing(const SdfScene &sdf, const AABB &sdf_bbox, const float3 &pos, const float3 &dir, 
                          float3 *surface_pos)
  {
    constexpr float EPS = 1e-5;
    float t = 0;
    float tFar = 1e4;
    if (!sdf_bbox.contains(pos))
    {
      if (!sdf_bbox.intersects(pos, dir, &t, &tFar))
        return false;
    }
    int iter = 0;
    float d = get_dist(sdf, pos + t*dir);
    while (iter < 1000 && d > EPS && t < tFar)
    {
      t += d + EPS;
      d = get_dist(sdf, pos + t*dir);
      iter++;
    }
    if (surface_pos)
      *surface_pos = pos + t*dir;
    //fprintf(stderr, "st %d (%f %f %f)", iter, p0.x, p0.y, p0.z);
    return d <= EPS;
  }

  void save_sdf_scene(const SdfScene &scene, const std::string &path)
  {
    std::ofstream fs(path, std::ios::binary); 
    unsigned c_count = scene.conjunctions.size();
    unsigned o_count = scene.objects.size();
    unsigned p_count = scene.parameters.size();

    fs.write((const char *)(&c_count), sizeof(unsigned));
    fs.write((const char *)(&o_count), sizeof(unsigned));
    fs.write((const char *)(&p_count), sizeof(unsigned));

    fs.write((const char *)scene.conjunctions.data(), c_count*sizeof(SdfConjunction));
    fs.write((const char *)scene.objects.data(), o_count*sizeof(SdfObject));
    fs.write((const char *)scene.parameters.data(), p_count*sizeof(float));
    fs.flush();
    fs.close();
  }

  void load_sdf_scene(SdfScene &scene, const std::string &path)
  {
    std::ifstream fs(path, std::ios::binary); 
    unsigned c_count = 0;
    unsigned o_count = 0;
    unsigned p_count = 0;

    fs.read((char *)(&c_count), sizeof(unsigned));
    fs.read((char *)(&o_count), sizeof(unsigned));
    fs.read((char *)(&p_count), sizeof(unsigned));

    assert(c_count > 0);
    assert(o_count > 0);
    assert(p_count > 0);
    scene.conjunctions.resize(c_count);
    scene.objects.resize(o_count);
    scene.parameters.resize(p_count);

    fs.read((char *)scene.conjunctions.data(), c_count*sizeof(SdfConjunction));
    fs.read((char *)scene.objects.data(), o_count*sizeof(SdfObject));
    fs.read((char *)scene.parameters.data(), p_count*sizeof(float));
    fs.close();
  }