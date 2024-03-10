#pragma once

#include "common_utils/body.h"
#include "common_utils/blk.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
#include "tinyEngine/texture.h"
#include <vector>
#include "graphics_utils/volumetric_occlusion.h"

namespace visualizer
{
    void body_to_model(Body *b, Mesh *m, bool fixed_tc = false, float4 tc = float4(1,0,0,0));
    void heightmap_to_model(Heightmap &h, Mesh *m, float2 detailed_size, float2 full_size, float precision,
                            int LODs);
    void box_to_model(Box *b, Mesh *m);
    void ellipsoid_to_model(Ellipsoid *b, Mesh *m, int sectors, int stacks, bool smooth = true);
    void cylinder_to_model(Cylinder *b, Mesh *m, int sectors);

    void visualize_light_voxels(LightVoxelsCube *voxels, Mesh *m, float3 shift = float3(0,0,0), float3 scale = float3(1,1,1));
    void visualize_light_voxels(LightVoxelsCube *voxels, Mesh *m, float3 pos, float3 size, float3 step, float dot_size, 
                                  float threshold = 0, float3 shift = float3(0,0,0), float3 scale = float3(1,1,1),
                                  int mip = 0);
    void visualize_aabb(AABB &box, Mesh *m, float3 &color);
    void visualize_aabb(::std::vector<AABB> &boxes, Mesh *m, ::std::vector<float3> &colors);
    void simple_mesh_to_model_332(const std::vector<float> &verts, Mesh *m);
};

namespace model_loader
{
  void load_default_blk();
  bool create_model_from_block(Block &bl, ComplexModel &mod);
  Model *create_model_by_name(std::string name, Texture &tex);
  Model *create_debug_box_model();
  Model *create_simple_grass_model();
  Model *load_model_from_obj(std::string name, Texture &tex);
  Model *load_model_from_obj_directly(std::string obj_filename);
  void normalize_model(Model *m);
  void save_model_to_obj(const Model *m, const std::string &filename);
  
  extern std::string base_path;
  extern Block obj_models_blk;
  extern bool obj_models_blk_loaded;
};
