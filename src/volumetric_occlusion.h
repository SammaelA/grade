#pragma once
#include <glm/glm.hpp>
#include <list>

struct LightVoxelsCube
{
    public:
    LightVoxelsCube(glm::vec3 center, glm::vec3 size, float base_size, float light_precision);
    ~LightVoxelsCube();
    void set_occluder(glm::vec3 pos, float strenght);
    void set_point_light(glm::vec3 pos, float strength);
    void set_directed_light(glm::vec3 direction, float strength);
    float get_occlusion(glm::vec3 pos);
    glm::vec3 get_dir_to_bright_place(glm::vec3 pos, float *occlusion);
    void print_average_occlusion();
    private:
    struct LightParams
    {
        int penumbraDepth = 5;
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
    float voxel_size;
    int vox_x,vox_y,vox_z;
    float *voxels;
    std::list<Light> point_lights;
    std::list<Light> directed_lights;
    float fast_voxel_occlusion(glm::ivec3 voxel);
    float voxel_occlusion(glm::ivec3 voxel);
    glm::ivec3 pos_to_voxel(glm::vec3 pos);
    bool in_voxel_cube(glm::ivec3 voxel);
    int v_to_i(glm::ivec3 voxel);
    float sum_occlusion = 0.0;
    float occ_count = 0.0;
};