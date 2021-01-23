#include "volumetric_occlusion.h"
#include <cstdio>
#include <vector>
#include "tinyEngine/utility.h"
#include "distribution.h"

LightVoxelsCube::LightVoxelsCube(glm::vec3 center, glm::vec3 size, float base_size, float light_precision)
{
    float voxel_size = base_size / light_precision;
    lightParams.penumbraDepth = MAX(1, light_precision * lightParams.penumbraDepth);
    lightParams.searchDepth = MAX(1, light_precision * lightParams.searchDepth);
    this->center = center;
    this->voxel_size = voxel_size;
    set_directed_light(glm::vec3(0, 1, 0), 1);
    vox_x = (int)(size.x / voxel_size);
    vox_y = (int)(size.y / voxel_size);
    vox_z = (int)(size.z / voxel_size);

    int count = (2 * vox_x + 1) * (2 * vox_y + 1) * (2 * vox_z + 1);
    //debug("trying to create light voxels cube with %dx%dx%d  %d voxels\n", vox_x, vox_y, vox_z, count);
    voxels = new float[count];
    //debug("successfully created light voxels cube with %d voxels\n", count);
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
    set_occluder_pyramid(pos,strenght);
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
    return get_occlusion_trilinear(pos);
}
glm::vec3 LightVoxelsCube::get_dir_to_bright_place_ext(glm::vec3 pos, int light_test_r, float *occlusion = nullptr)
{
    float min_occ = 1e10;
    const float BIAS = 0.001;
    glm::vec3 min_shift(1, 1, 1);
    std::vector<glm::vec3> min_shifts;
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
    min_shift = min_shifts[urandi(0,min_sz)];
    if (occlusion)
        *occlusion = min_occ;
    return (glm::normalize(min_shift + glm::vec3(0.0, 0.0001, 0.0)));
}
glm::vec3 LightVoxelsCube::get_dir_to_bright_place(glm::vec3 pos, float *occlusion = nullptr)
{
    return get_dir_to_bright_place_ext(pos, lightParams.searchDepth, occlusion);
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
int LightVoxelsCube::v_to_i(int x, int y, int z)
{
    return (2 * vox_x + 1) * (2 * vox_y + 1) * (vox_z + z) + (2 * vox_x + 1) * (vox_y + y) + (vox_x + x);
}
int LightVoxelsCube::v_to_i(glm::ivec3 voxel)
{
    return (2 * vox_x + 1) * (2 * vox_y + 1) * (vox_z + voxel.z) + (2 * vox_x + 1) * (vox_y + voxel.y) + (vox_x + voxel.x);
}
void LightVoxelsCube::set_occluder_pyramid(glm::vec3 pos, float strenght)
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
void LightVoxelsCube::set_occluder_simple(glm::vec3 pos, float strenght)
{
    glm::ivec3 voxel = pos_to_voxel(pos);
    if (in_voxel_cube(voxel))
        voxels[v_to_i(voxel)] += strenght;
}
void LightVoxelsCube::set_occluder_trilinear(glm::vec3 pos, float strenght)
{
    glm::vec3 flpos = glm::vec3(floorf(pos.x),floorf(pos.y),floorf(pos.z));
    glm::vec3 dp = pos - flpos;
    glm::ivec3 voxel = pos_to_voxel(flpos);
    int x = voxel.x;
    int y = voxel.y;
    int z = voxel.z;
    #define C(i,j,k) voxels[v_to_i(x+i,y+j,z+k)]
    if (in_voxel_cube(voxel) && in_voxel_cube(voxel + glm::ivec3(1,1,1)))
    {
        C(0,0,0) += strenght*(1 - dp.x)*(1 - dp.y)*(1 - dp.z);
        C(0,0,1) += strenght*(1 - dp.x)*(1 - dp.y)*(dp.z);
        C(0,1,0) += strenght*(1 - dp.x)*(dp.y)*(1 - dp.z);
        C(0,1,1) += strenght*(1 - dp.x)*(dp.y)*(dp.z);
        C(1,0,0) += strenght*(dp.x)*(1 - dp.y)*(1 - dp.z);
        C(1,0,1) += strenght*(dp.x)*(1 - dp.y)*(dp.z);
        C(1,1,0) += strenght*(dp.x)*(dp.y)*(1 - dp.z);
        C(1,1,1) += strenght*(dp.x)*(dp.y)*(dp.z);

    }
    else if (in_voxel_cube(voxel))
    {
        C(0,0,0) += strenght;
    }
    
}
float LightVoxelsCube::get_occlusion_view_ray(glm::vec3 pos)
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
float LightVoxelsCube::get_occlusion_simple(glm::vec3 pos)
{
    glm::ivec3 voxel = pos_to_voxel(pos);
    if (in_voxel_cube(voxel))
    {
        int v = v_to_i(voxel);
        return voxels[v];
    }
    else
    {
        return 1000; 
    }
}
float LightVoxelsCube::get_occlusion_trilinear(glm::vec3 pos)
{
    glm::vec3 flpos = glm::vec3(floorf(pos.x),floorf(pos.y),floorf(pos.z));
    glm::vec3 dp = pos - flpos;
    glm::ivec3 voxel = pos_to_voxel(flpos);
    int x = voxel.x;
    int y = voxel.y;
    int z = voxel.z;
    #define C(i,j,k) voxels[v_to_i(x+i,y+j,z+k)]
    if (in_voxel_cube(voxel) && in_voxel_cube(voxel + glm::ivec3(1,1,1)))
    {
        float c00 = C(0,0,0)*(1-dp.x) + C(1,0,0)*dp.x;
        float c01 = C(0,0,1)*(1-dp.x) + C(1,0,1)*dp.x;
        float c10 = C(0,1,0)*(1-dp.x) + C(1,1,0)*dp.x;
        float c11 = C(0,1,1)*(1-dp.x) + C(1,1,1)*dp.x;

        float c0 = c00*(1-dp.y) + c10*dp.y;
        float c1 = c01*(1-dp.y) + c11*dp.y;

        float c = c0*(1-dp.z) + c1*dp.z;

        return c;
    }
    else if (in_voxel_cube(voxel))
    {
        return C(0,0,0);
    }
    else
    {
        return 1000; 
    }
}
float LightVoxelsCube::NMSE(LightVoxelsCube *B)
{
    if (voxels_count() != B->voxels_count())
    {
        logerr("Unable to find NMSE between different-sized voxel cubes sizes:(%d %d %d) and (%d %d %d)",
                vox_x,vox_y,vox_z,B->vox_x,B->vox_y,B->vox_z);
        return 1;
    }
    float MSE = 0.0;
    float SQ = 0.001;
    for (int i=0;i<voxels_count();i++)
    {
        MSE += abs((voxels[i] - B->voxels[i])*(voxels[i] - B->voxels[i]));
        SQ += voxels[i]*voxels[i] + (B->voxels[i])*(B->voxels[i]);
    }
    //debug("NMSE = %f %f %f\n", MSE, SQ, MSE/SQ);
    return MSE/SQ;
}
void LightVoxelsCube::add_body(Body *b, float opacity)
{
    glm::vec3 dv = glm::vec3(voxel_size,voxel_size,voxel_size);
    for (int i=-vox_x;i<=vox_x;i++)
    {
        for (int j=-vox_y;j<=vox_y;j++)
        {
            for (int k=-vox_z;k<=vox_z;k++)
            {
                if (b->in_body(center + voxel_size*glm::vec3(i,j,k)))   
                {
                    logerr("in body %d %d %d",i,j,k);
                    voxels[v_to_i(i,j,k)] += opacity;
                }
            }
        }
    }
}