#pragma once
#include <glm/glm.hpp>
#include <list>
#include "body.h"
#include "terrain.h"
struct LightVoxelsCube
{
public:
    LightVoxelsCube(glm::vec3 center, glm::vec3 size, float base_size, float light_precision);
    LightVoxelsCube(glm::vec3 center, glm::ivec3 sizes, float voxel_size);
    LightVoxelsCube(LightVoxelsCube *source, glm::vec3 pos, glm::vec3 sizes);
    LightVoxelsCube(LightVoxelsCube *source, glm::ivec3 vox_pos, glm::ivec3 vox_sizes);
    ~LightVoxelsCube();
    void replace_occluder_voxel(glm::ivec3 voxel, float strength);
    void set_occluder_voxel(glm::ivec3 voxel, float strength);
    void set_occluder(glm::vec3 pos, float strenght);
    void set_occluder_pyramid(glm::vec3 pos, float strenght);
    void set_occluder_simple(glm::vec3 pos, float strenght);
    void set_occluder_trilinear(glm::vec3 pos, float strenght);
    void set_point_light(glm::vec3 pos, float strength);
    void set_directed_light(glm::vec3 direction, float strength);
    float get_occlusion_voxel(glm::ivec3 voxel);
    float get_occlusion(glm::vec3 pos);
    float get_occlusion_view_ray(glm::vec3 pos);
    float get_occlusion_simple(glm::vec3 pos);
    float get_occlusion_trilinear(glm::vec3 pos);
    glm::vec3 get_dir_to_bright_place(glm::vec3 pos, float *occlusion);
    glm::vec3 get_dir_to_bright_place_ext(glm::vec3 pos, int steps, float *occlusion);
    glm::vec3 get_dir_to_bright_place_ray(glm::vec3 pos, float length, int rays_sqrt, float *occlusion);
    void print_average_occlusion();
    void add_body(Body *b, float opacity = 1e9, bool solid = true);
    void add_heightmap(Heightmap &h);
    float NMSE(LightVoxelsCube *B);
    glm::vec3 get_center();
    float get_voxel_size();
    glm::ivec3 get_vox_sizes();
private:
    struct LightParams
    {
        int penumbraDepth = 20;
        float penumbraWidthInc = 1;
        float penumbraDepthDecay = 0.4;
        float penumbraWidthDecay = 0.8;
        float searchDepth = 5;
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
    const float voxel_size;
    int vox_x, vox_y, vox_z;
    int count;
    float *voxels;
    std::list<Light> point_lights;
    std::list<Light> directed_lights;
    float fast_voxel_occlusion(glm::ivec3 voxel);
    float voxel_occlusion(glm::ivec3 voxel);
    glm::ivec3 pos_to_voxel(glm::vec3 pos);
    glm::vec3 voxel_to_pos(glm::ivec3 voxel);
    bool in_voxel_cube(glm::ivec3 voxel);
    int v_to_i(glm::ivec3 voxel);
    int v_to_i(int x, int y, int z);
    int voxels_count() {return vox_x*vox_y*vox_z;}
    float sum_occlusion = 0.0;
    float occ_count = 0.0;
};