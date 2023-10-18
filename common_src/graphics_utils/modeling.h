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
    void body_to_model(Body *b, Mesh *m, bool fixed_tc = false, glm::vec4 tc = glm::vec4(1,0,0,0));
    void heightmap_to_model(Heightmap &h, Mesh *m, glm::vec2 detailed_size, glm::vec2 full_size, float precision,
                            int LODs);
    void box_to_model(Box *b, Mesh *m);
    void ellipsoid_to_model(Ellipsoid *b, Mesh *m, int sectors, int stacks, bool smooth = true);
    void cylinder_to_model(Cylinder *b, Mesh *m, int sectors);

    void visualize_light_voxels(LightVoxelsCube *voxels, Mesh *m, glm::vec3 shift = glm::vec3(0,0,0), glm::vec3 scale = glm::vec3(1,1,1));
    void visualize_light_voxels(LightVoxelsCube *voxels, Mesh *m, glm::vec3 pos, glm::vec3 size, glm::vec3 step, float dot_size, 
                                  float threshold = 0, glm::vec3 shift = glm::vec3(0,0,0), glm::vec3 scale = glm::vec3(1,1,1),
                                  int mip = 0);
    void visualize_aabb(AABB &box, Mesh *m, glm::vec3 &color);
    void visualize_aabb(::std::vector<AABB> &boxes, Mesh *m, ::std::vector<glm::vec3> &colors);
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
  void save_model_to_obj(const Model *m, const std::string &filename);
  
  extern std::string base_path;
  extern Block obj_models_blk;
  extern bool obj_models_blk_loaded;
};
