#include "volumetric_occlusion.h"
#include <cstdio>
#include <vector>
#include "tinyEngine/utility.h"
LightVoxelsCube::LightVoxelsCube(glm::vec3 center, glm::vec3 size, float base_size, float light_precision)
{
    float voxel_size = base_size / light_precision;
    lightParams.penumbraDepth = MAX(1, light_precision * lightParams.penumbraDepth);
    lightParams.searchDepth = MAX(1, light_precision * lightParams.searchDepth);
    this->center = center;
    this->voxel_size = voxel_size;
    set_directed_light(glm::vec3(0, 1, 0), 1);
    vox_x = (int)(size.x / voxel_size);
    vox_y = (int)(size.x / voxel_size);
    vox_z = (int)(size.x / voxel_size);

    int count = (2 * vox_x + 1) * (2 * vox_y + 1) * (2 * vox_z + 1);
    debug("trying to create light voxels cube with %dx%dx%d  %d voxels\n", vox_x, vox_y, vox_z, count);
    voxels = new float[count];
    debug("successfully created light voxels cube with %d voxels\n", count);
}
void LightVoxelsCube::print_average_occlusion()
{
    debugl(2, "average_occ %f\n", sum_occlusion / occ_count);
}
LightVoxelsCube::~LightVoxelsCube()
{
    delete (voxels);
}
void LightVoxelsCube::set_occluder(glm::vec3 pos, float strenght)
{
    float base_str = strenght;
    glm::ivec3 voxel = pos_to_voxel(pos);
    for (int i = 0; i < (int)lightParams.penumbraDepth; i++)
    {
        int wd = i * lightParams.penumbraWidthInc;
        for (int j = -wd; j <= wd; j++)
        {
            for (int k = -wd; k <= wd; k++)
            {
                glm::ivec3 vx = voxel + glm::ivec3(j, -i, k);
                if (in_voxel_cube(vx))
                {
                    voxels[v_to_i(vx)] += base_str * powf(lightParams.penumbraWidthDecay, abs(j) + abs(k));
                }
                else
                {
                    debugl(2, "missed pos = %f %f %f center = %f %f %f size = %f %f %f", pos.x, pos.y, pos.z, center.x, center.y, center.z, vox_x * voxel_size, vox_y * voxel_size, vox_z * voxel_size);
                }
            }
        }
        base_str *= lightParams.penumbraDepthDecay;
    }
}
inline void LightVoxelsCube::set_point_light(glm::vec3 pos, float strength)
{
    point_lights.push_back(Light(pos, strength));
}
inline void LightVoxelsCube::set_directed_light(glm::vec3 direction, float strength)
{
    directed_lights.push_back(Light(glm::normalize(direction), strength));
}
float LightVoxelsCube::get_occlusion(glm::vec3 pos)
{

    float sum_occ = 0.0;
    int prev_v = -1;
    for (auto &light : directed_lights)
    {
        glm::vec3 start = (pos - center) / voxel_size;
        glm::vec3 step = glm::vec3(0, 1, 0) / (2 * voxel_size);
        glm::ivec3 voxel = glm::ivec3((int)start.x, (int)start.y, (int)start.z);
        while (in_voxel_cube(voxel))
        {
            int v = v_to_i(voxel);
            if (v != prev_v)
            {
                sum_occ += voxels[v];
                prev_v = v;
            }
            start += step;
            voxel = glm::ivec3((int)start.x, (int)start.y, (int)start.z);
        }
    }
    sum_occlusion += sum_occ;
    occ_count += 1;
    return sum_occ;
}
glm::vec3 LightVoxelsCube::get_dir_to_bright_place(glm::vec3 pos, float *occlusion = nullptr)
{
    float min_occ = 1e10;
    const float BIAS = 0.001;
    glm::vec3 min_shift(1, 1, 1);
    std::vector<glm::vec3> min_shifts;
    int light_test_r = lightParams.searchDepth;
    for (int i = -light_test_r; i <= light_test_r; i++)
    {
        for (int j = -light_test_r; j <= light_test_r; j++)
        {
            for (int k = -light_test_r; k <= light_test_r; k++)
            {
                glm::vec3 shift = voxel_size * glm::vec3(i, j, k);
                float cur_occ = get_occlusion(pos + shift);
                if (abs(cur_occ - min_occ) < BIAS)
                {
                    min_shifts.push_back(shift);
                }
                if (cur_occ < min_occ)
                {
                    min_occ = cur_occ;
                    min_shifts.clear();
                    min_shifts.push_back(shift);
                }
            }
        }
    }
    int min_sz = min_shifts.size();
    min_shift = min_shifts[rand() % min_sz];
    if (occlusion)
        *occlusion = min_occ;
    return (glm::normalize(min_shift + glm::vec3(0.0, 0.0001, 0.0)));
}
inline bool LightVoxelsCube::in_voxel_cube(glm::ivec3 voxel)
{
    return (abs(voxel.x) <= vox_x) && (abs(voxel.y) <= vox_y) && (abs(voxel.z) <= vox_z);
}
inline float LightVoxelsCube::fast_voxel_occlusion(glm::ivec3 voxel)
{
    return voxels[v_to_i(voxel)];
}
inline float LightVoxelsCube::voxel_occlusion(glm::ivec3 voxel)
{
    if (in_voxel_cube(voxel))
        return voxels[v_to_i(voxel)];
}
glm::ivec3 LightVoxelsCube::pos_to_voxel(glm::vec3 pos)
{
    pos -= center;
    glm::ivec3 voxel = glm::ivec3((int)(pos.x / voxel_size), (int)(pos.y / voxel_size), (int)(pos.z / voxel_size));
    return voxel;
}
inline int LightVoxelsCube::v_to_i(glm::ivec3 voxel)
{
    return (2 * vox_x + 1) * (2 * vox_y + 1) * (vox_z + voxel.z) + (2 * vox_x + 1) * (vox_y + voxel.y) + (vox_x + voxel.x);
}