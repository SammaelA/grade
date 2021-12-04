#pragma once
#include <glm/glm.hpp>
#include <list>
#include "core/body.h"
#include "graphics_utils/terrain.h"
#include "common_utils/distribution.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
#include <vector>
struct LightVoxelsCube : public Countable
{
public:
    LightVoxelsCube(glm::vec3 center, glm::vec3 size, float base_size, float light_precision, int mip_levels = 1, int mip_decrease = 2);
    LightVoxelsCube(glm::vec3 center, glm::ivec3 sizes, float voxel_size, int mip_levels = 1, int mip_decrease = 2);
    LightVoxelsCube(LightVoxelsCube *source);
    LightVoxelsCube(LightVoxelsCube *source, glm::vec3 pos, glm::vec3 sizes, int size_decrease = 1, 
                    glm::vec2 min_max = glm::vec2(0,1e10));
    LightVoxelsCube(LightVoxelsCube *source, glm::ivec3 vox_pos, glm::ivec3 vox_sizes, int size_decrease,
                    glm::vec2 min_max);
    ~LightVoxelsCube();
    void clear();
    void replace_occluder_voxel(glm::ivec3 voxel, float strength);
    void set_occluder_voxel(glm::ivec3 voxel, float strength);
    void set_occluder_voxel_mip(glm::ivec3 voxel, float strenght, int mip);
    void set_occluder(glm::vec3 pos, float strenght);
    void set_occluder_pyramid(glm::vec3 pos, float strenght);
    void set_occluder_pyramid2(glm::vec3 pos, float strenght, float pow_b, int max_r);
    void set_occluder_simple(glm::vec3 pos, float strenght);
    void set_occluder_simple_mip(glm::vec3 pos, float strenght, int mip);
    void set_occluder_trilinear(glm::vec3 pos, float strenght);
    void set_point_light(glm::vec3 pos, float strength);
    void set_directed_light(glm::vec3 direction, float strength);
    float get_occlusion_voxel(glm::ivec3 voxel);
    float get_occlusion_voxel_unsafe(int x, int y, int z);
    float get_occlusion(glm::vec3 pos);
    float get_occlusion_view_ray(glm::vec3 pos, glm::vec3 light);
    float get_occlusion_view_ray(glm::vec3 pos);
    float get_occlusion_simple(glm::vec3 pos);
    float get_occlusion_projection(glm::vec3 pos);
    float get_occlusion_simple_mip(glm::vec3 pos, int mip);
    float get_occlusion_trilinear(glm::vec3 pos);
    float get_occlusion_trilinear_mip(glm::vec3 pos, int mip);
    float get_occlusion_cone(glm::vec3 pos, glm::vec3 dir, glm::vec3 n, float H, float R, 
                             int num_samples, Uniform &phi_distr, Uniform &h_distr, Uniform &l_distr);
    glm::vec3 get_dir_to_bright_place(glm::vec3 pos, float *occlusion);
    glm::vec3 get_dir_to_bright_place_cone(glm::vec3 pos, float r, int cones, float *occlusion);
    glm::vec3 get_dir_to_bright_place_ext(glm::vec3 pos, int steps, float *occlusion);
    glm::vec3 get_dir_to_bright_place_ray(glm::vec3 pos, float length, int rays_sqrt, float *occlusion);
    void prepare_mip_levels();
    void print_average_occlusion();
    void add_body(Body *b, float opacity = 1e9, bool solid = true, float infl_distance = 0, float base_infl_occ = 0);
    void add_AABB(AABB &box, bool precise, float opacity = 1e9);
    void add_heightmap(Heightmap &h);
    void add_voxels_cube(LightVoxelsCube *cube);
    void calculte_precise_occlusion_from_bodies();
    void get_data(float **data, glm::ivec3 &size);
    AABB get_bbox();
    float NMSE(LightVoxelsCube *B);
    glm::vec3 get_center();
    float get_voxel_size();
    glm::ivec3 get_vox_sizes();
    int get_mip_count();
    int get_mip_decrease();
    long get_size_cnt() {return count;}
    
private:
    struct LightParams
    {
        int penumbraDepth = 20;
        float penumbraWidthInc = 1;
        float penumbraDepthDecay = 0.4;
        float penumbraWidthDecay = 0.8;
        float searchDepth = 5;
        float sunPositions = 4;
    } lightParams;
    struct Light
    {
        Light() {}
        Light(glm::vec3 vec, float intensity)
        {
            this->vec = vec;
            this->intensity = intensity;
        }
        glm::vec3 vec;
        float intensity;
    };
    glm::vec3 center;
    float voxel_size;
    int vox_x, vox_y, vox_z;
    int count;
    int block_size = 1;
    int block_x, block_y, block_z;
    int block_cnt = 1;
    int mip_levels = 1;
    int mip_decrease = 2;
    std::vector<int> mip_offsets;
    std::vector<int> mip_size_decrease;
    std::vector<glm::ivec3> mip_vox_xyz;
    float *voxels;
    std::list<Light> point_lights;
    std::list<Light> directed_lights;
    float fast_voxel_occlusion(glm::ivec3 voxel);
    float voxel_occlusion(glm::ivec3 voxel);
    glm::ivec3 pos_to_voxel(glm::vec3 pos);
    glm::ivec3 pos_to_block(glm::vec3 pos);
    glm::vec3 voxel_to_pos(glm::ivec3 voxel);
    void fill_blocks(glm::ivec3 from, glm::ivec3 to, float val);
    void fill_voxels(glm::ivec3 from, glm::ivec3 to, float val);
    bool in_voxel_cube(glm::ivec3 voxel);
    bool in_voxel_cube(int x, int y, int z);
    int v_to_i(glm::ivec3 voxel);
    int v_to_i(int x, int y, int z);
    int v_to_i_mip(glm::ivec3 voxel, int mip);
    int v_to_i_mip(int x, int y, int z, int mip);
    int v_to_i_mip_no_offset(glm::ivec3 voxel, int mip);
    int v_to_i_mip_no_offset(int x, int y, int z, int mip);
    int voxels_count() {return vox_x*vox_y*vox_z;}

    float sum_occlusion = 0.0;
    float occ_count = 0.0;
};