#pragma once
#include <vector>
#include <string>
#include <cstring> 
#include <cassert> 
#include <map>
#include "../src/utils.h"
#include "common_utils/LiteMath_ext.h"

namespace diff_render
{
  using LiteMath::float2;
  using LiteMath::float3;
  using LiteMath::float4;
  using LiteMath::uint4;
  using LiteMath::int4;
  using LiteMath::uint3;
  using LiteMath::int3;
  using LiteMath::int2;
  using LiteMath::uint2;
  using LiteMath::ushort4;
  using LiteMath::uchar4;
  using LiteMath::clamp;
  using LiteMath::float4x4;

struct CPUTexture
{
  CPUTexture() = default;
  CPUTexture(const LiteImage::Image2D<float3> &img)
  {
    w = img.width();
    h = img.height();
    channels = 3;
    data = ::std::vector<float>((float*)img.data(), (float*)img.data() + w*h*channels);
  }

  inline int pixel_to_offset(int x, int y) const { return channels*(y*w + x); }
  inline int pixel_to_offset(int2 pixel) const { return channels*(pixel.y*w + pixel.x); }
  /**
  \brief UNSAFE ACCESS!!!

  */
  const float *get(int x, int y) const
  {
    return data.data() + pixel_to_offset(x,y); 
  }
  ::std::vector<float> data;
  int w,h,channels;
};

/**
\brief input/output mesh
*/
struct TriangleMesh 
{
  TriangleMesh() = default;
  TriangleMesh(const ::std::vector<float3> &_vertices, const ::std::vector<float3> &_colors, const ::std::vector<unsigned> &_indices = {})
  {
    vertices = _vertices;
    colors = _colors;
    indices = _indices;
  }
  
  TriangleMesh(const ::std::vector<float3> &_vertices, const ::std::vector<float2> &_tc, const ::std::vector<unsigned> &_indices = {})
  {
    vertices = _vertices;
    tc = _tc;
    indices = _indices;
  }

  inline int vertex_count() const { return vertices.size(); }
  inline int face_count() const { return indices.size()/3; }

  //vertex attributes, some of them might be empty
  ::std::vector<float3>     vertices;
  ::std::vector<float3>     colors;
  ::std::vector<float2>     tc;
  ::std::vector<float3>     normals;
  ::std::vector<float3>     tangents;

  ::std::vector<unsigned>   indices;

  ::std::vector<CPUTexture> textures; // an arbitrary number of textures
};

void transform(TriangleMesh &mesh, const LiteMath::float4x4 &transform);

struct PointLight
{
  PointLight() = default;
  PointLight(const float3 &col, float inten, const float3 &_pos)
  {
    color = normalize(col);
    intensity = inten;
    pos = _pos;
  }

  float3 color;
  float intensity;
  float3 pos;
};



struct TransformR
{
  TransformR():
    TransformR(float3(0,0,0), float3(0,0,0), 1.0f)
  { };
  TransformR(float3 tr, float3 rot, float sc)
  {
    data.s.translate = tr;
    data.s.rotate = rot;
    data.s.scale = sc;
  };
  float4x4 to_mat() const
  {
    return LiteMath::translate4x4(data.s.translate)*
           LiteMath::rotate4x4Z(data.s.rotate.z)*LiteMath::rotate4x4Y(data.s.rotate.y)*LiteMath::rotate4x4X(data.s.rotate.x)*
           LiteMath::scale4x4({data.s.scale, data.s.scale, data.s.scale});
  }
  struct TRS
  {
    float3 translate;
    float3 rotate;
    float scale;
  };

  union TRU
  {
    TRU():s(){};
    TRS s;
    float M[7];
  } data;
};

/**
\brief base scene description without auxilary data for differentiable rendering
*/
struct Scene
{
public:
  enum MeshInstancingType
  {
    SINGLE,
    BASE_TRANSFORM,
    RESTRICTED_TRANSFORM
  };
  void add_mesh(const TriangleMesh &mesh, const ::std::vector<float4x4> &transform = {float4x4()}, ::std::string name = "")
  {
    meshes.push_back(mesh);
    mesh_instancing_types.push_back(transform.size() == 1 ? SINGLE : BASE_TRANSFORM);
    transforms.push_back(transform);
    restricted_transforms.emplace_back();
    meshes_by_name.emplace(name, meshes.size()-1);
    prepared = false;
  }

  void add_mesh(const TriangleMesh &mesh, const ::std::vector<TransformR> &transform, ::std::string name = "")
  {
    meshes.push_back(mesh);
    mesh_instancing_types.push_back(RESTRICTED_TRANSFORM);
    transforms.emplace_back();
    restricted_transforms.push_back(transform);
    meshes_by_name.emplace(name, meshes.size()-1);
    prepared = false;
  }

  void set_mesh(const TriangleMesh &mesh, int id, const ::std::vector<TransformR> &transform)
  {
    assert(id >= 0 && id < meshes.size());
    meshes[id] = mesh;
    mesh_instancing_types[id] = RESTRICTED_TRANSFORM;
    restricted_transforms[id] = transform;
    prepared = false;
  }

  void set_mesh(const TriangleMesh &mesh, int id, const ::std::vector<float4x4> &transform = {float4x4()})
  {
    assert(id >= 0 && id < meshes.size());
    meshes[id] = mesh;
    mesh_instancing_types[id] = transform.size() == 1 ? SINGLE : BASE_TRANSFORM;
    transforms[id] = transform;
    prepared = false;
  }

  inline unsigned get_index(unsigned mesh_id, unsigned instance_id, unsigned vertex_id) const 
  { 
    return preparedData.indices[mesh_id][instance_id_mapping[instance_id]][vertex_id]; 
  }
  inline unsigned get_vertex_n(unsigned mesh_id, unsigned vertex_id) const
  {
    return meshes[mesh_id].indices[vertex_id];
  }
  //id is an index returned by get_index() function
  inline float3 get_pos(unsigned id)     const { return preparedData.vertices[id]; }
  inline float3 get_pos_orig(unsigned id)const { return preparedData.orig_vertices[id]; }
  inline float3 get_color(unsigned id)   const { return preparedData.colors[id]; }
  inline float3 get_norm(unsigned id)    const { return preparedData.normals[id]; }
  inline float3 get_tang(unsigned id)    const { return preparedData.tangents[id]; }
  inline float2 get_tc(unsigned id)      const { return preparedData.tc[id]; }
  
  inline const CPUTexture &get_tex(unsigned mesh_id, unsigned n) const { return meshes[mesh_id].textures[n]; }
  
  inline const TriangleMesh &get_mesh(unsigned n) const { return meshes[n]; }
  inline const ::std::vector<TriangleMesh> &get_meshes() const { return meshes; }
  inline TriangleMesh &get_mesh_modify(unsigned n) { prepared = false; return meshes[n]; }
  inline ::std::vector<TriangleMesh> &get_meshes_modify() { prepared = false; return meshes; }
  
  inline const ::std::vector<float4x4>  &get_transform(unsigned n) const { return transforms[n]; }
  inline const ::std::vector<::std::vector<float4x4>>  &get_transforms() const { return transforms; }
  inline const ::std::vector<TransformR>  &get_restricted_transform(unsigned n) const { return restricted_transforms[n]; }
  inline const ::std::vector<::std::vector<TransformR>>  &get_restricted_transforms() const { return restricted_transforms; }
  inline const ::std::vector<float4x4>  &get_transform_inv(unsigned n) const { return transforms_inv[n]; }
  inline const ::std::vector<::std::vector<float4x4>>  &get_transforms_inv() const { return transforms_inv; }
  inline ::std::vector<float4x4>  &get_transform_modify(unsigned n) { prepared = false; return transforms[n];  }
  inline ::std::vector<::std::vector<float4x4>>  &get_transforms_modify() { prepared = false; return transforms; }
  inline ::std::vector<TransformR>  &get_restricted_transform_modify(unsigned n) { prepared = false; return restricted_transforms[n];  }
  inline ::std::vector<::std::vector<TransformR>>  &get_restricted_transforms_modify() { prepared = false; return restricted_transforms; }
  inline MeshInstancingType get_instancing_type(unsigned n) const { return mesh_instancing_types[n]; }
  inline unsigned indices_size() const 
  {
    unsigned c = 0;
    for (auto &m : meshes)
      c+=m.indices.size();
    return c;
  }

  void restore_meshes(bool restore_normals, bool restore_tangents, bool transform_to_unindexed_mesh);

  //applies all transorms to meshes and puts them into one big structure
  void prepare_for_render() const;

  void get_prepared_mesh(TriangleMesh &mesh) const;
  void set_instance_id_mapping(const ::std::vector<int> &v) const { instance_id_mapping = v; }
  ::std::vector<int> get_instance_id_mapping() const { return instance_id_mapping; }
  void invalidate_prepared_scene() const { prepared = false; }

protected:
  ::std::vector<TriangleMesh> meshes;
  ::std::vector<MeshInstancingType> mesh_instancing_types;
  ::std::vector<::std::vector<TransformR>> restricted_transforms;
  mutable ::std::vector<::std::vector<float4x4>> transforms;
  mutable ::std::vector<::std::vector<float4x4>> transforms_inv;
  ::std::map<::std::string, int> meshes_by_name; //position in meshes vector
  //index in this vector is the instanceId from ray tracer
  //value is the instance id to take values from indices
  mutable ::std::vector<int> instance_id_mapping;

  float3 ambient_light_color = float3(0,0,0);
  float3 environment_light_mult = float3(1,1,1);
  CPUTexture environment_light_texture;
  ::std::vector<PointLight> point_lights;


  mutable bool prepared = false;
  mutable struct PreparedData
  {
    //3-dimentional array id = indices[mesh_id][instance_id][vertex_id] 
    ::std::vector<::std::vector<::std::vector<unsigned>>> indices;

    ::std::vector<float3>     vertices;
    ::std::vector<float3>     orig_vertices;
    ::std::vector<float3>     colors;
    ::std::vector<float2>     tc;
    ::std::vector<float3>     normals;
    ::std::vector<float3>     tangents;
  } preparedData;
};
} // namespace diff_render
