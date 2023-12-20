#pragma once

#include <vector>
#include <string>
#include <cstring> // for memset
#include <random>  // for random gen 
#include <cassert> 

#include "LiteMath.h"
#include "Image2d.h"
#include "utils.h"
#include "scene.h"
namespace diff_render
{
enum class SHADING_MODEL {UNDEFINED = 0,
                          SILHOUETTE = 1,
                          VERTEX_COLOR = 2,
                          TEXTURE_COLOR = 3, 
                          LAMBERT = 4,
                          PHONG = 5,
                          GGX = 6,
                          PATH_TEST = 7};

typedef float GradReal;

class DScene;
extern void PrintAndCompareGradients(const DScene& grad1, const DScene& grad2);
class DMesh
{
public:
  friend class DScene;
  friend void PrintAndCompareGradients(const DScene& grad1, const DScene& grad2);
  DMesh(){};
  DMesh(SHADING_MODEL mode, const TriangleMesh &mesh, int _instances, Scene::MeshInstancingType instance_type, int m_id)
  {
    assert(mode != SHADING_MODEL::UNDEFINED);
    mesh_id = m_id;

    bool dpos = true;
    bool dcolor = mode == SHADING_MODEL::VERTEX_COLOR;
    bool dtransform = instance_type != Scene::RESTRICTED_TRANSFORM;
    bool dresttransform = instance_type == Scene::RESTRICTED_TRANSFORM;
    bool dtextures = !(mode == SHADING_MODEL::SILHOUETTE || mode == SHADING_MODEL::VERTEX_COLOR);
    
    int tex_size = 0;
    for (auto &t : mesh.textures)
      tex_size += t.w*t.h*t.channels;
    int sz = 3*mesh.vertices.size()*dpos + 3*mesh.vertices.size()*dcolor + dtextures*tex_size + 
             dtransform*TRANSFORM_SIZE*_instances + dresttransform*RESTRICTED_TRANSFORM_SIZE*_instances;
    data.resize(sz, 0);

    vertices = mesh.vertices.size();
    instances = _instances;
    int offset = 0;
    if (dpos)
    {
      pos_ptr = data.data() + offset;
      offset += 3*mesh.vertices.size();
    }
    if (dcolor)
    {
      color_ptr = data.data() + offset;
      offset += 3*mesh.vertices.size();
    }
    if (dtextures)
    {
      for (auto &t : mesh.textures)
      {
        tex_ptrs.push_back(data.data() + offset);
        tex_params.push_back({t.w, t.h, t.channels});
        offset += t.w*t.h*t.channels;
      }
    }
    if (dtransform)
    {
      transform_ptr = data.data() + offset;
      offset += TRANSFORM_SIZE*_instances;
    }
    if (dresttransform)
    {
      restricted_transform_ptr = data.data() + offset;
      offset += RESTRICTED_TRANSFORM_SIZE*_instances;      
    }
  }
  ~DMesh() = default;
 
  DMesh& operator=(const DMesh& dmesh)
  {
    copy(dmesh);
    return *this;
  }
  DMesh(const DMesh& dmesh)
  {
    copy(dmesh);
  }
  void copy(const DMesh& dmesh)
  {   
    if (this == &dmesh)
      return;
    vertices = dmesh.vertices;
    instances = dmesh.instances;
    data = dmesh.data;
    tex_params = dmesh.tex_params;
    mesh_id = dmesh.mesh_id;

    //restore points from offsets
    if (dmesh.pos_ptr)
      pos_ptr = data.data() + (dmesh.pos_ptr - dmesh.data.data());
    
    if (dmesh.color_ptr)
      color_ptr = data.data() + (dmesh.color_ptr - dmesh.data.data());
    
    tex_ptrs.resize(dmesh.tex_ptrs.size(), nullptr);
    for (int i=0; i<dmesh.tex_ptrs.size();i++)
      tex_ptrs[i] = data.data() + (dmesh.tex_ptrs[i] - dmesh.data.data());
    
    if (dmesh.transform_ptr)
      transform_ptr = data.data() + (dmesh.transform_ptr - dmesh.data.data());

    if (dmesh.restricted_transform_ptr)
      restricted_transform_ptr = data.data() + (dmesh.restricted_transform_ptr - dmesh.data.data());
  }

  inline GradReal *pos(int index) { return pos_ptr + 3*index; }
  inline GradReal *color(int index) { return color_ptr + 3*index; }
  inline GradReal &tex(int tex_n, int px, int py, int channel) 
  { 
    return *(tex_ptrs[tex_n] + tex_params[tex_n].z*(tex_params[tex_n].x*py + px) + channel);
  }
  inline GradReal *tex(int tex_n, int offset) { return tex_ptrs[tex_n] + offset; }
  inline int3 get_tex_info(int tex_n) const { return tex_params[tex_n]; } //(x,y,channels)
  inline GradReal *transform_mat(int instance_n) { return transform_ptr ? transform_ptr + TRANSFORM_SIZE*instance_n : nullptr; }
  inline GradReal *restricted_transform(int instance_n) { return restricted_transform_ptr ? restricted_transform_ptr + RESTRICTED_TRANSFORM_SIZE*instance_n : nullptr; }

  inline GradReal *full_data() { return data.data(); }
  inline int full_size() { return data.size(); }
  inline int vertex_count() { return vertices; }
  inline int instance_count() { return instances; }
  inline int tex_count() { return tex_params.size(); }
  inline int get_mesh_id() { return mesh_id; }

  static constexpr int TRANSFORM_SIZE = 12;
  static constexpr int RESTRICTED_TRANSFORM_SIZE = 7;
protected:
  ::std::vector<GradReal> data;
  GradReal *pos_ptr = nullptr;
  GradReal *color_ptr = nullptr;
  ::std::vector<int3> tex_params;
  ::std::vector<GradReal *> tex_ptrs;
  GradReal *transform_ptr = nullptr;
  GradReal *restricted_transform_ptr = nullptr;
  int vertices = 0;
  int instances = 0;
  int mesh_id = -1;
};
class DScene
{
public:
  DScene(){};
  DScene(const Scene &scene, SHADING_MODEL mode, const ::std::vector<int> &dmesh_ids)
  {
    reset(scene, mode, dmesh_ids);
  }
  ~DScene() = default;
 
  DScene& operator=(const DScene& dscene)
  {
    copy(dscene);
    return *this;
  }
  DScene(const DScene& dscene)
  {
    copy(dscene);
  }
  void copy(const DScene& dscene)
  {
    if (this == &dscene)
      return;
    
    shading_model = dscene.shading_model;
    dmeshes.resize(dscene.dmeshes.size());
    dmeshes_by_id.resize(dscene.dmeshes_by_id.size(), nullptr);
    int k = 0;
    for (int i=0;i<dmeshes_by_id.size();i++)
    {
      if (dscene.dmeshes_by_id[i])
      {
        dmeshes[k] = *(dscene.dmeshes_by_id[i]);
        dmeshes_by_id[i] = dmeshes.data() + k;
        k++;
      }
    }
  }
  void reset(const Scene &scene, SHADING_MODEL mode, const ::std::vector<int> &dmesh_ids)
  {
    bool need_reset = (mode != shading_model) || (dmesh_ids.size() != dmeshes.size());
    if (!need_reset)
    {
      for (auto mesh_id : dmesh_ids)
        if (get_dmesh(mesh_id)->vertices != scene.get_mesh(mesh_id).vertices.size())
          need_reset = true;
    }
    if (need_reset)
    {
      assert(dmesh_ids.size()>0);
      scene.prepare_for_render();
      shading_model = mode;
      int max_mesh_id = scene.get_meshes().size();

      for (auto mesh_id : dmesh_ids)
      {
        max_mesh_id = ::std::max(max_mesh_id, mesh_id);
        dmeshes.push_back(DMesh(mode, scene.get_mesh(mesh_id), scene.get_transform(mesh_id).size(), scene.get_instancing_type(mesh_id), mesh_id));
      }
      dmeshes_by_id.resize(max_mesh_id+1, nullptr);
      for (int i=0;i<dmesh_ids.size();i++)
        dmeshes_by_id[dmesh_ids[i]] = dmeshes.data()+i;
    }
    else
    {
      clear();
    }
  }
  DMesh *get_dmesh(int mesh_id) { return dmeshes_by_id[mesh_id]; }
  const ::std::vector<DMesh> &get_dmeshes() const { return dmeshes; }
  ::std::vector<DMesh> &get_dmeshes() { return dmeshes; }
  void clear()
  {
    for (auto &dm : dmeshes)
      ::std::fill(dm.data.begin(), dm.data.end(), 0.0f);
  }
  int full_size()
  {
    int sz = 0;
    for (auto &m : dmeshes)
      sz += m.full_size();
    return sz;
  }
  SHADING_MODEL get_shading_model() { return shading_model; }
private:
  ::std::vector<DMesh> dmeshes;
  ::std::vector<DMesh*> dmeshes_by_id;
  SHADING_MODEL shading_model = SHADING_MODEL::UNDEFINED;
};

static inline float3 SummOfPixels(const Img& a_image) 
{
  const auto& color = a_image.vector();
  double summ[3] = {0.0, 0.0, 0.0};
  for(size_t i=0;i<color.size();i++) {
    summ[0] += double(color[i].x);
    summ[1] += double(color[i].y);
    summ[2] += double(color[i].z); 
  }
  return float3(summ[0], summ[1], summ[2]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}