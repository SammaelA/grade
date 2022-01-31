#include "graphics_utils/volumetric_occlusion.h"
#include <cstdio>
#include <vector>
#include "common_utils/utility.h"
#include "common_utils/distribution.h"
#include "common_utils/sun.h"
#include <chrono>
#include <atomic>

using namespace glm;

std::atomic<int> sum_memory(0);
std::atomic<int> sum_allocs(1);
glm::ivec3 vox_sizes(glm::vec3 sizes, float voxel_size)
{
    return glm::ivec3(MAX(1,sizes.x/voxel_size),MAX(1,sizes.y/voxel_size),MAX(1,sizes.z/voxel_size));
}
int LVC_divisible(int a, int b)
{
    if (b == 1)
        return a;
    int c = a;
    int cnt = 0;
    while (((2*c + 1) % b) && (cnt <= b))
    {
        c++;
        cnt++;
    }
    if (cnt <= b)
    {
        return c;
    }
    else
    {
        logerr("cannot find div %d %d",a,b);
        return a;
    }
}
LightVoxelsCube::LightVoxelsCube(glm::vec3 center, glm::vec3 size, float base_size, float light_precision,
                                 int mip_levels, int mip_decrease):
LightVoxelsCube(center, vox_sizes(size,base_size / light_precision), base_size / light_precision, mip_levels, mip_decrease)
{
    lightParams.penumbraDepth = MAX(1, lightParams.penumbraDepth/voxel_size);
    lightParams.searchDepth = MAX(1, light_precision * lightParams.searchDepth);
}
LightVoxelsCube::LightVoxelsCube(glm::vec3 cent, glm::ivec3 sizes, float vox_size, int _mip_levels, int _mip_decrease):
voxel_size(vox_size)
{
    center = cent;
    mip_levels = _mip_levels;
    mip_decrease = _mip_decrease;
    block_size = (mip_levels == 1) ? 5 : 1;
    block_cnt = block_size*block_size*block_size;
    vox_x = LVC_divisible(sizes.x, block_size);
    vox_y = LVC_divisible(sizes.y, block_size);
    vox_z = LVC_divisible(sizes.z, block_size);

    block_x = ceil((2*vox_x + 1.0)/block_size);
    block_y = ceil((2*vox_y + 1.0)/block_size);
    block_z = ceil((2*vox_z + 1.0)/block_size);

    set_directed_light(glm::vec3(0, 1, 0), 1);

    count = 0;
    for (int i=0;i<mip_levels;i++)
    {
        int decrease = powf(mip_decrease,i);
        int mip_x = ceil((float)(block_x*block_size)/decrease);
        int mip_y = ceil((float)(block_y*block_size)/decrease);
        int mip_z = ceil((float)(block_z*block_size)/decrease);
        mip_offsets.push_back(count);
        mip_vox_xyz.push_back(glm::ivec3(mip_x,mip_y,mip_z));
        mip_size_decrease.push_back(decrease);
        count+= mip_x*mip_y*mip_z;
    }
    voxels = new float[count];
    fill(0);
    sum_memory.fetch_add(count*sizeof(float));
    sum_allocs.fetch_add(1);
    //debug("[%f %f %f] created voxels cube %d x %d x %d voxels, %dx%dx%d blocks\n", center.x, center.y, center.z,
    //      2*vox_x + 1, 2*vox_y + 1, 2*vox_z + 1, block_x, block_y, block_z);
}

void LightVoxelsCube::fill(float val)
{
    std::fill(voxels,voxels+count,0);
}
LightVoxelsCube::LightVoxelsCube(LightVoxelsCube *source):
LightVoxelsCube(source,
                glm::ivec3(0,0,0),
                source->get_vox_sizes())
{

}
LightVoxelsCube::LightVoxelsCube(LightVoxelsCube *source, glm::vec3 pos, glm::vec3 sizes, int size_decrease, glm::vec2 min_max):
LightVoxelsCube(source,
                source->pos_to_voxel(pos),
                vox_sizes(sizes,source->voxel_size*size_decrease),
                size_decrease,
                min_max)

{
    glm::ivec3 vsz = vox_sizes(sizes, source->voxel_size);
    glm::ivec3 vox_pos = source->pos_to_voxel(pos);
    if (size_decrease > 1 && (vsz.x % size_decrease || vsz.y % size_decrease || vsz.z % size_decrease))
    {
        debug("warning: trying to get slice [cent: %d %d %d, size: %d %d %d] of volumetric occlusion \
               with size_decrease = %d. It will make the real size smaller %d %d %d\n",vox_pos.x,vox_pos.y,vox_pos.z,
               vsz.x,vsz.y,vsz.z, size_decrease, source->vox_x, source->vox_y, source->vox_z);
    }
    center = pos;
}
LightVoxelsCube::LightVoxelsCube(LightVoxelsCube *source, glm::ivec3 vox_pos, glm::ivec3 src_vox_sizes, int size_decrease, glm::vec2 min_max):
LightVoxelsCube(source->voxel_to_pos(vox_pos),src_vox_sizes/size_decrease,source->voxel_size*size_decrease,
                source->mip_levels, source->mip_decrease)
{
    glm::ivec3 vox_sizes = glm::ivec3(vox_x, vox_y, vox_z);
    //logerr("creating voxels %d %d %d from %d %d %d",get_vox_sizes().x,get_vox_sizes().y,get_vox_sizes().z,
    //source->get_vox_sizes().x,source->get_vox_sizes().y,source->get_vox_sizes().z);

    for (int i = -vox_sizes.x; i <= +vox_sizes.x; i++)
    {
        for (int j = -vox_sizes.y; j <= +vox_sizes.y; j++)
        {
            for (int k = -vox_sizes.z; k <= +vox_sizes.z; k++)
            {
                float occ = 0;
                int sm = 0;
                for (int x0 = 0; x0 < size_decrease; x0+=(size_decrease - 1)/2)
                {
                    for (int y0 = 0; y0 < size_decrease; y0+=(size_decrease - 1)/2)
                    {
                        for (int z0 = 0; z0 < size_decrease; z0+=(size_decrease - 1)/2)
                        {
                            int x = size_decrease*i + x0 + vox_pos.x;
                            int y = size_decrease*j + y0 + vox_pos.y;
                            int z = size_decrease*k + z0 + vox_pos.z;
                            if (source->in_voxel_cube(x,y,z))
                            {
                                float f = source->get_occlusion_voxel_unsafe(x,y,z);
                                if (f >= min_max.x && f < min_max.y)
                                {
                                    occ += f;
                                    sm ++;
                                }
                            }
                        }
                    }
                }
                voxels[v_to_i(i,j,k)] = sm > 0 ? occ/sm : 0;
            }
        }
    }
    prepare_mip_levels();
}
void LightVoxelsCube::get_data(float **data, glm::ivec3 &size)
{
    //TODO - adjust to block structure
    *data = voxels;
    for (int i=0;i<voxels_count();i++)
    {
        if (voxels[i] < 0)
        {
            logerr("VOL error %d %f ",i,voxels[i]);
        }
    }
    size.x = vox_x;
    size.y = vox_y;
    size.z = vox_z;
}
void LightVoxelsCube::set_occluder_voxel(glm::ivec3 voxel, float strength)
{
    if (in_voxel_cube(voxel))
    {
        voxels[v_to_i(voxel)] += strength;
    }
}
void LightVoxelsCube::replace_occluder_voxel(glm::ivec3 voxel, float strength)
{
    if (in_voxel_cube(voxel))
    {
        voxels[v_to_i(voxel)] = strength;
    }
}
float LightVoxelsCube::get_occlusion_voxel(glm::ivec3 voxel)
{
    if (in_voxel_cube(voxel))
    {
        return voxels[v_to_i(voxel)];
    }
    else 
        return 1e9;
}

float LightVoxelsCube::get_occlusion_voxel_unsafe(int x, int y, int z)
{
    return voxels[v_to_i(x, y, z)];
}

glm::vec3 LightVoxelsCube::get_center()
{
    return center;
}
AABB LightVoxelsCube::get_bbox()
{
    return AABB(voxel_to_pos(-get_vox_sizes()), voxel_to_pos(get_vox_sizes()));
}
float LightVoxelsCube::get_voxel_size()
{
    return voxel_size;
}
glm::vec3 LightVoxelsCube::voxel_to_pos(glm::ivec3 voxel)
{
    return center + voxel_size*glm::vec3(voxel.x,voxel.y,voxel.z);
}
glm::ivec3 LightVoxelsCube::get_vox_sizes()
{
    return glm::ivec3(vox_x,vox_y,vox_z);
}

int LightVoxelsCube::get_mip_count()
{
    return mip_levels;
}

int LightVoxelsCube::get_mip_decrease()
{
    return mip_decrease;
}

void LightVoxelsCube::print_average_occlusion()
{
    debugl(2, "average_occ %f\n", sum_occlusion / occ_count);
}
void LightVoxelsCube::clear()
{
    if (voxels)
    {
        int count = (2 * vox_x + 1) * (2 * vox_y + 1) * (2 * vox_z + 1);
        sum_memory.fetch_add(-count*sizeof(float));
        sum_allocs.fetch_add(-1);
        delete[] voxels;
        voxels = nullptr;
    }
}
LightVoxelsCube::~LightVoxelsCube()
{
    clear();
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
    if (pos.x != pos.x || pos.y != pos.y || pos.z != pos.z)//if pos in NAN
        return 1e9;
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
    min_shift = glm::normalize(min_shift + glm::vec3(0.0, 0.0001, 0.0));
    if (occlusion)
        *occlusion = get_occlusion(pos - voxel_size*light_test_r*min_shift);
    return min_shift;
}
float LightVoxelsCube::get_occlusion_cone(glm::vec3 pos, glm::vec3 dir, glm::vec3 n, float H, float R, 
                         int num_samples, Uniform &phi_distr, Uniform &h_distr, Uniform &l_distr)
{
    float occ = 0;
    float w_sum = 0;
    for (int i = 0; i<num_samples; i++)
    {
        float phi = phi_distr.get();
        float h = h_distr.get();
        float l = l_distr.get();
        float wl = l/R;
        l /= (h/H);
        float wh = (h/H)*(h/H);
        glm::vec3 nr = glm::rotate(glm::mat4(1),phi,dir)*glm::vec4(n,0);
        glm::vec3 sample = pos + h*dir + l*nr;
        occ += wl*wh*get_occlusion_trilinear(sample);
        w_sum += wl*wh;
    }
    if (w_sum < 1e-6)
        return 10000;
    else
        return occ/w_sum;
}
glm::vec3 LightVoxelsCube::get_dir_to_bright_place_cone(glm::vec3 pos, float r, int cones, float *occlusion)
{
    const int num_samples = 64;
    const float BIAS = 0.1;
    float alpha = PI/10;
    float cone_R = r*sin(alpha); 
    float cone_H = r;
    Uniform phi_distr = Uniform(0,2*PI);
    Uniform h_distr = Uniform(0,cone_H);
    Uniform l_distr = Uniform(0,cone_R);
    Uniform box_distr = Uniform(-r,r);
    float min_occ = 1e10;

    glm::vec3 min_shift = glm::normalize(glm::vec3(1,1,1));
    std::vector<glm::vec3> min_shifts;
    glm::vec3 prev_dir = glm::normalize(glm::vec3(1,1,1));
    for (int i = 0; i < cones; i++)
    {
        bool in_sph = false;
        glm::vec3 dir; 
        while (!in_sph)
        {
            dir.x = box_distr.get();
            dir.y = box_distr.get();
            dir.z = box_distr.get();
            in_sph = (glm::length(dir) < r) && (glm::dot(glm::normalize(dir),prev_dir) < 0.9);
        }

        dir = glm::normalize(dir);
        glm::vec3 n = glm::cross(dir,prev_dir);

        float occ = get_occlusion_cone(pos,dir,n,cone_H,cone_R,num_samples,phi_distr,h_distr,l_distr);
        if (occ < (1-BIAS)*min_occ)
        {
            min_shifts.clear();
            min_shifts.push_back(dir);
            min_occ = occ;
        }
        else if (occ < (1+BIAS)*min_occ)
        {
            min_shifts.push_back(dir);
            min_occ = MIN(occ,min_occ);
        }
        prev_dir = dir;
    }

    int min_sz = min_shifts.size();
    min_shift = min_shifts[urandi(0,min_sz)];
    min_shift = glm::normalize(min_shift + glm::vec3(0.0, 0.0001, 0.0));
    if (occlusion)
        *occlusion = min_occ;
    return min_shift;
}
glm::vec3 LightVoxelsCube::get_dir_to_bright_place(glm::vec3 pos, float *occlusion = nullptr)
{
    return get_dir_to_bright_place_ext(pos, lightParams.searchDepth, occlusion);
}
inline bool LightVoxelsCube::in_voxel_cube(glm::ivec3 voxel)
{
    return (abs(voxel.x) <= vox_x) && (abs(voxel.y) <= vox_y) && (abs(voxel.z) <= vox_z);
}
inline bool LightVoxelsCube::in_voxel_cube(int x, int y, int z)
{
    return (x >= -vox_x) && (x <= vox_x) && (y >= -vox_y) && (y <= vox_y) && (z >= -vox_z) && (z <= vox_z);
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
    glm::ivec3 voxel = glm::ivec3((pos.x / voxel_size), (pos.y / voxel_size), (pos.z / voxel_size));
    return voxel;
}
glm::ivec3 LightVoxelsCube::pos_to_block(glm::vec3 pos)
{
    pos -= center;
    return glm::ivec3(((int)(pos.x / voxel_size) + vox_x) / block_size,
                      ((int)(pos.y / voxel_size) + vox_y) / block_size,
                      ((int)(pos.z / voxel_size) + vox_z) / block_size);
}

#define LIN(x,y,z, mx, my, mz) ((mx)*(my)*(z) + (mx)*(y) + (x))

void LightVoxelsCube::fill_blocks(glm::ivec3 from, glm::ivec3 to, float val)
{
    //logerr("filling blocks %d %d %d -- %d %d %d", from.x, from.y, from.z, to.x, to.y, to.z);
    if (from.x >= to.x || from.y >= to.y || from.z >= to.z)
        return;
    int base_start = block_cnt*LIN(from.x, from.y, from.z, block_x, block_y, block_z);
    int size = block_cnt*(to.x - from.x + 1);
    std::fill_n(voxels + base_start, size, val);

    for (int i = from.y; i<= to.y;i++)
    {
        for (int j = from.z; j<= to.z;j++)
        {
            int start = block_cnt*LIN(from.x, i, j, block_x, block_y, block_z);
            memcpy(voxels + start, voxels + base_start, sizeof(float)*size);
        }
    }
}
void LightVoxelsCube::fill_voxels(glm::ivec3 from, glm::ivec3 to, float val)
{
    for (int i = from.y; i<= to.y;i++)
    {
        for (int j = from.z; j<= to.z;j++)
        {
            for (int k = from.x; k<= to.x;k++)
            {
                voxels[v_to_i(k, i, j)] = val;
            }
        }
    }
}

int LightVoxelsCube::v_to_i(int x, int y, int z)
{
    int bl_x = (x + vox_x) / block_size;
    int bl_y = (y + vox_y) / block_size;
    int bl_z = (z + vox_z) / block_size;

    int self_x = (x + vox_x) % block_size;
    int self_y = (y + vox_y) % block_size;
    int self_z = (z + vox_z) % block_size;

    return block_cnt*LIN(bl_x, bl_y, bl_z, block_x, block_y, block_z)+
           LIN(self_x, self_y, self_z, block_size, block_size, block_size);
}
int LightVoxelsCube::v_to_i(glm::ivec3 voxel)
{
    int bl_x = (voxel.x + vox_x) / block_size;
    int bl_y = (voxel.y + vox_y) / block_size;
    int bl_z = (voxel.z + vox_z) / block_size;

    int self_x = (voxel.x + vox_x) % block_size;
    int self_y = (voxel.y + vox_y) % block_size;
    int self_z = (voxel.z + vox_z) % block_size;
    return block_cnt*LIN(bl_x, bl_y, bl_z, block_x, block_y, block_z)+
           LIN(self_x, self_y, self_z, block_size, block_size, block_size);
}

/*
int LightVoxelsCube::v_to_i(int x, int y, int z)
{
    return (2 * vox_x + 1) * (2 * vox_y + 1) * (vox_z + z) + (2 * vox_x + 1) * (vox_y + y) + (vox_x + x);
}
int LightVoxelsCube::v_to_i(glm::ivec3 voxel)
{
    return (2 * vox_x + 1) * (2 * vox_y + 1) * (vox_z + voxel.z) + (2 * vox_x + 1) * (vox_y + voxel.y) + (vox_x + voxel.x);
}
*/
int LightVoxelsCube::v_to_i_mip(glm::ivec3 voxel, int mip)
{
    return mip == 0 ? v_to_i(voxel) :
          (mip_offsets[mip] + 
           mip_vox_xyz[mip].x * mip_vox_xyz[mip].y * (vox_z + voxel.z)/mip_size_decrease[mip] + 
           mip_vox_xyz[mip].x * (vox_y + voxel.y)/mip_size_decrease[mip] + 
           (vox_x + voxel.x)/mip_size_decrease[mip]);   
}
int LightVoxelsCube::v_to_i_mip(int x, int y, int z, int mip)
{
    return mip == 0 ? v_to_i(x,y,z) : 
          (mip_offsets[mip] + 
           mip_vox_xyz[mip].x * mip_vox_xyz[mip].y * ((vox_z + z)/mip_size_decrease[mip]) + 
           mip_vox_xyz[mip].x * ((vox_y + y)/mip_size_decrease[mip]) + 
           (vox_x + x)/mip_size_decrease[mip]);    
}
int LightVoxelsCube::v_to_i_mip_no_offset(glm::ivec3 voxel, int mip)
{
    return mip == 0 ? v_to_i(voxel) :
          (mip_offsets[mip] + 
           mip_vox_xyz[mip].x * mip_vox_xyz[mip].y * (voxel.z) + 
           mip_vox_xyz[mip].x * (voxel.y) + 
           (voxel.x));   
}
int LightVoxelsCube::v_to_i_mip_no_offset(int x, int y, int z, int mip)
{
    return mip == 0 ? v_to_i(x,y,z) : 
          (mip_offsets[mip] + 
           mip_vox_xyz[mip].x * mip_vox_xyz[mip].y * (z) + 
           mip_vox_xyz[mip].x * (y) + 
           (x));    
}
void LightVoxelsCube::set_occluder_pyramid(glm::vec3 pos, float strenght)
{
    if (pos.x != pos.x || pos.y != pos.y || pos.z != pos.z)//if pos in NAN
        return;
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
                    float dist_sq = (i*i + j*j + k*k) + 1;
                    voxels[v_to_i(vx)] += base_str / dist_sq;
                }
            }
        }
        base_str *= lightParams.penumbraDepthDecay;
    }
}

//int cnt = 0;
//double tm = 0.0;

void LightVoxelsCube::set_occluder_pyramid2(glm::vec3 pos, float strenght, float pow_b, int max_r)
{
    #define CH 0.6
    float occ = 0.0;
    glm::ivec3 voxel = pos_to_voxel(pos);
    int x0 = voxel.x;
    int y0 = voxel.y;
    int z0 = voxel.z;
    for (int i = 0; voxel.y - i > -vox_y; i++)
    {
        int wd = MIN(i, max_r);
        if (wd == max_r) 
        {
            if (urand() < CH)
                continue;
            else 
                occ = 1/(1-CH);
        }
        else
            occ = 1;

        if (abs(pow_b - 2) < 1e-4)
            occ *= strenght / SQR(i + 2);
        else
            occ *= strenght * pow(i + 2, -pow_b);
        if (abs(occ) < 0.00025)
            return;
        int y = y0 - i;
        for (int j = MAX(z0 - wd, -vox_z); j <= MIN(z0 + wd, vox_z); j++)
        {
            for (int k = MAX(x0 - wd, -vox_x); k <= MIN(x0 + wd, vox_x); k++)
            {
                voxels[v_to_i(k, y, j)] += occ;
            }
        }
    }
}

void LightVoxelsCube::set_occluder_pyramid_tail_5(glm::ivec3 voxel, float strenght, int max_r, int rnd_seed)
{
    #define BLOCK_SIZE_5 5
    int x0 = voxel.x;
    int y0 = voxel.y;
    int z0 = voxel.z;
    int wd = max_r;
    int max_depth = sqrt(abs(strenght)/0.005) - 2;

    int x1 = MAX(x0 - wd, -vox_x);
    int y1 = MAX(-vox_y, y0 - max_depth);
    int z1 = MAX(z0 - wd, -vox_z);

    int bl1_x = (x1 + vox_x) / BLOCK_SIZE_5;
    int bl1_y = (y1 + vox_y) / BLOCK_SIZE_5;
    int bl1_z = (z1 + vox_z) / BLOCK_SIZE_5;

    int self1_x = (x1 + vox_x) % BLOCK_SIZE_5;
    int self1_y = (y1 + vox_y) % BLOCK_SIZE_5;
    int self1_z = (z1 + vox_z) % BLOCK_SIZE_5;

    int bl1_pos = LIN(bl1_x, bl1_y, bl1_z, block_x, block_y, block_z);

    int x2 = MIN(x0 + wd, vox_x);
    int y2 = MAX(y0 - wd, MAX(-vox_y, y0 - max_depth));
    int z2 = MIN(z0 + wd, vox_z);

    int bl2_x = (x2 + vox_x) / BLOCK_SIZE_5;
    int bl2_y = (y2 + vox_y) / BLOCK_SIZE_5;
    int bl2_z = (z2 + vox_z) / BLOCK_SIZE_5;

    int self2_x = (x2 + vox_x) % BLOCK_SIZE_5;
    int self2_y = (y2 + vox_y) % BLOCK_SIZE_5;
    int self2_z = (z2 + vox_z) % BLOCK_SIZE_5;

    bl1_z += rnd_seed % 2;

    for (int bl_z = bl1_z; bl_z <= bl2_z; bl_z+=2)
    {
        int vst_z = bl_z == bl1_z ? self1_z : 0;
        int ven_z = bl_z == bl2_z ? self2_z : BLOCK_SIZE_5 - 1;
        for (int bl_y = bl1_y; bl_y <= bl2_y; bl_y++)
        {
            int vst_y = bl_y == bl1_y ? self1_y : 0;
            int ven_y = bl_y == bl2_y ? self2_y : BLOCK_SIZE_5 - 1;
            for (int bl_x = bl1_x; bl_x <= bl2_x; bl_x++)
            {
                int vst_x = bl_x == bl1_x ? self1_x : 0;
                int ven_x = bl_x == bl2_x ? self2_x : BLOCK_SIZE_5 - 1;
                int bl_pos = LIN(bl_x, bl_y, bl_z, block_x, block_y, block_z);
                int bl_st_y = -bl_y*BLOCK_SIZE_5 + vox_y;
                for (int self_z = vst_z; self_z <= ven_z; self_z ++)
                {
                    int self_addr_z = block_cnt*bl_pos + BLOCK_SIZE_5*BLOCK_SIZE_5*self_z;
                    for (int self_y = vst_y; self_y <= ven_y; self_y++)
                    {
                        int self_addr_y = self_addr_z + BLOCK_SIZE_5*self_y;
                        float occ = strenght / SQR(bl_st_y - self_y + y0 + 2);
                        for (int self_x = vst_x; self_x <= ven_x; self_x++)
                        {
                            voxels[self_addr_y + self_x] += occ;
                        }
                    }
                }
            }
        }
    }
}

void LightVoxelsCube::set_occluder_pyramid_head_5(glm::ivec3 voxel, float strenght, int max_r)
    {
        #define BLOCK_SIZE_5 5
    int x0 = voxel.x;
    int y0 = voxel.y;
    int z0 = voxel.z;
        int wd = max_r - 1;
        int x1 = MAX(x0 - wd, -vox_x);
        int y1 = MAX(-vox_y, y0 - wd);
        int z1 = MAX(z0 - wd, -vox_z);

        int bl1_x = (x1 + vox_x) / BLOCK_SIZE_5;
        int bl1_y = (y1 + vox_y) / BLOCK_SIZE_5;
        int bl1_z = (z1 + vox_z) / BLOCK_SIZE_5;

        int self1_x = (x1 + vox_x) % BLOCK_SIZE_5;
        int self1_y = (y1 + vox_y) % BLOCK_SIZE_5;
        int self1_z = (z1 + vox_z) % BLOCK_SIZE_5;

        int bl1_pos = LIN(bl1_x, bl1_y, bl1_z, block_x, block_y, block_z);

        int x2 = MIN(x0 + wd, vox_x);
        int y2 = y0;
        int z2 = MIN(z0 + wd, vox_z);

        int bl2_x = (x2 + vox_x) / BLOCK_SIZE_5;
        int bl2_y = (y2 + vox_y) / BLOCK_SIZE_5;
        int bl2_z = (z2 + vox_z) / BLOCK_SIZE_5;

        int self2_x = (x2 + vox_x) % BLOCK_SIZE_5;
        int self2_y = (y2 + vox_y) % BLOCK_SIZE_5;
        int self2_z = (z2 + vox_z) % BLOCK_SIZE_5;    
        for (int bl_z = bl1_z; bl_z <= bl2_z; bl_z++)
        {
            int bl_z_sh = bl_z*BLOCK_SIZE_5 - vox_z - z0;
            int vst_z = bl_z == bl1_z ? self1_z : 0;
            int ven_z = bl_z == bl2_z ? self2_z : BLOCK_SIZE_5 - 1;
            for (int bl_y = bl1_y; bl_y <= bl2_y; bl_y++)
            {
                int bl_y_sh = y0 - bl_y*BLOCK_SIZE_5 + vox_y;
                int vst_y = bl_y == bl1_y ? self1_y : 0;
                int ven_y = bl_y == bl2_y ? self2_y : BLOCK_SIZE_5 - 1;
                for (int bl_x = bl1_x; bl_x <= bl2_x; bl_x++)
                {
                    int vst_x = bl_x == bl1_x ? self1_x : 0;
                    int ven_x = bl_x == bl2_x ? self2_x : BLOCK_SIZE_5 - 1;
                    int bl_pos = LIN(bl_x, bl_y, bl_z, block_x, block_y, block_z);
                    int bl_x_sh = bl_x*BLOCK_SIZE_5 - vox_x - x0;
                    for (int self_z = vst_z; self_z <= ven_z; self_z ++)
                    {
                        int self_addr_z = block_cnt*bl_pos + BLOCK_SIZE_5*BLOCK_SIZE_5*self_z;
                        int self_sh_z = abs(bl_z_sh + self_z);
                        int en_local_y = MIN(ven_y,bl_y_sh-self_sh_z);
                        for (int self_y = vst_y; self_y <= en_local_y; self_y++)
                        {
                            int self_addr_y = self_addr_z + BLOCK_SIZE_5*self_y;
                            int self_sh_y = bl_y_sh - self_y;
                            float occ = strenght / SQR(bl_y_sh - self_y + 2);
                            int st_local_x = MAX(vst_x, -self_sh_y - bl_x_sh);
                            int en_local_x = MIN(ven_x, -bl_x_sh + self_sh_y);
                            for (int self_x = st_local_x; self_x <= en_local_x; self_x++)
                            {
                                voxels[self_addr_y + self_x] += occ;
                                /*
                                if (length(pos - glm::vec3(17-200,19-200,23-200))<0.5)
                                {
                                    logerr("[%d-%d][%d-%d][%d-%d]voxels[%d %d %d] += %f shifts %d %d %d", 
                                    st_local_x, en_local_x, vst_y, en_local_y, vst_z, ven_z,
                                    bl_x*block_size +self_x - vox_x, 
                                    bl_y*block_size + self_y - vox_y, bl_z*block_size +self_z - vox_z, occ, bl_x_sh + self_x, -self_sh_y,
                                    bl_z_sh + self_z);
                                }*/
                                
                            }
                        }
                    }
                }
            }
        }
    }

void LightVoxelsCube::set_occluder_pyramid_fast(glm::vec3 pos, float strenght, int max_r, int rnd_seed)
{
    float occ = 0.0;
    glm::ivec3 voxel = pos_to_voxel(pos);
    int x0 = voxel.x;
    int y0 = voxel.y;
    int z0 = voxel.z;

    if (block_size == 5)
    {
        set_occluder_pyramid_head_5(voxel, strenght, max_r);
        set_occluder_pyramid_tail_5(voxel, strenght, max_r, rnd_seed);
        return;
    }

    for (int i = 0; i < MIN(voxel.y+vox_y, max_r);i++)
    {
        occ = strenght / SQR(i + 2);
        int y = y0 - i;
        for (int j = MAX(z0 - i, -vox_z); j <= MIN(z0 + i, vox_z); j++)
        {
            for (int k = MAX(x0 - i, -vox_x); k <= MIN(x0 + i, vox_x); k++)
            {
                voxels[v_to_i(k, y, j)] += occ;
            }
        }
    }

    int wd = max_r;
    int max_depth = sqrt(abs(strenght)/0.00025) - 2;

    int x1 = MAX(x0 - wd, -vox_x);
    int y1 = MAX(-vox_y, y0 - max_depth);
    int z1 = MAX(z0 - wd, -vox_z);

    int bl1_x = (x1 + vox_x) / block_size;
    int bl1_y = (y1 + vox_y) / block_size;
    int bl1_z = (z1 + vox_z) / block_size;

    int self1_x = (x1 + vox_x) % block_size;
    int self1_y = (y1 + vox_y) % block_size;
    int self1_z = (z1 + vox_z) % block_size;

    int bl1_pos = LIN(bl1_x, bl1_y, bl1_z, block_x, block_y, block_z);

    int x2 = MIN(x0 + wd, vox_x);
    int y2 = MAX(y0 - wd, MAX(-vox_y, y0 - max_depth));
    int z2 = MIN(z0 + wd, vox_z);

    int bl2_x = (x2 + vox_x) / block_size;
    int bl2_y = (y2 + vox_y) / block_size;
    int bl2_z = (z2 + vox_z) / block_size;

    int self2_x = (x2 + vox_x) % block_size;
    int self2_y = (y2 + vox_y) % block_size;
    int self2_z = (z2 + vox_z) % block_size;

    bl1_z += rnd_seed % 2;

    for (int bl_z = bl1_z; bl_z <= bl2_z; bl_z+=2)
    {
        int vst_z = bl_z == bl1_z ? self1_z : 0;
        int ven_z = bl_z == bl2_z ? self2_z : block_size - 1;
        for (int bl_y = bl1_y; bl_y <= bl2_y; bl_y++)
        {
            int vst_y = bl_y == bl1_y ? self1_y : 0;
            int ven_y = bl_y == bl2_y ? self2_y : block_size - 1;
            for (int bl_x = bl1_x; bl_x <= bl2_x; bl_x++)
            {
                int vst_x = bl_x == bl1_x ? self1_x : 0;
                int ven_x = bl_x == bl2_x ? self2_x : block_size - 1;
                int bl_pos = LIN(bl_x, bl_y, bl_z, block_x, block_y, block_z);
                int bl_st_y = -bl_y*block_size + vox_y;
                for (int self_z = vst_z; self_z <= ven_z; self_z ++)
                {
                    int self_addr_z = block_cnt*bl_pos + block_size*block_size*self_z;
                    for (int self_y = vst_y; self_y <= ven_y; self_y++)
                    {
                        int self_addr_y = self_addr_z + block_size*self_y;
                        float occ = strenght / SQR(bl_st_y - self_y + y0 + 2);
                        for (int self_x = vst_x; self_x <= ven_x; self_x++)
                        {
                            voxels[self_addr_y + self_x] += occ;
                        }
                    }
                }
            }
        }
    }
}

void LightVoxelsCube::set_occluder_simple(glm::vec3 pos, float strenght)
{
    glm::ivec3 voxel = pos_to_voxel(pos);
    if (in_voxel_cube(voxel))
        voxels[v_to_i(voxel)] += strenght;
}

void LightVoxelsCube::set_occluder_voxel_mip(glm::ivec3 voxel, float strenght, int mip)
{
    mip = CLAMP(mip, 0, mip_levels);
    if (in_voxel_cube(voxel))
        voxels[v_to_i_mip(voxel, mip)] += strenght;
}

void LightVoxelsCube::set_occluder_simple_mip(glm::vec3 pos, float strenght, int mip)
{
    mip = CLAMP(mip, 0, mip_levels);
    glm::ivec3 voxel = pos_to_voxel(pos);
    if (in_voxel_cube(voxel))
        voxels[v_to_i_mip(voxel, mip)] += strenght;
}

void LightVoxelsCube::set_occluder_trilinear(glm::vec3 pos, float strenght)
{
    if (pos.x != pos.x || pos.y != pos.y || pos.z != pos.z)//if pos in NAN
        return;
    glm::vec3 flpos = glm::vec3(floorf(pos.x),floorf(pos.y),floorf(pos.z));
    glm::vec3 dp = pos - flpos;
    glm::ivec3 voxel = pos_to_voxel(flpos);
    int x = voxel.x;
    int y = voxel.y;
    int z = voxel.z;
    if (dp.x < 0 || dp.x > 1 || dp.y < 0 ||dp.y > 1|| dp.z<0 || dp.z >1)
    {
        logerr("trilinear error %f %f %f pos =  %f %f %f",dp.x,dp.y,dp.z,pos.x,pos.y,pos.z);
    }
    #define C(i,j,k) voxels[v_to_i(x+i,y+j,z+k)]
    if (false && in_voxel_cube(voxel) && in_voxel_cube(voxel + glm::ivec3(1,1,1)))
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

        if (strenght<0 || C(0,0,0) < 0)
        {
            logerr("wrong value %f %f",strenght,C(0,0,0));
        }
        else
                C(0,0,0) += strenght;
    }
    
}
float LightVoxelsCube::get_occlusion_view_ray(glm::vec3 pos, glm::vec3 light)
{
    int prev_v = -1;
    float sum_occ = 0;
    glm::vec3 start = (pos - center) / voxel_size;
    glm::vec3 step = light / (2 * voxel_size);
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
    return sum_occ;
}
float LightVoxelsCube::get_occlusion_view_ray(glm::vec3 pos)
{
    float sum_occ = 0.0;
    for (auto &light : directed_lights)
    {
        sum_occ += get_occlusion_view_ray(pos, light.vec);
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
float LightVoxelsCube::get_occlusion_projection(glm::vec3 pos)
{
    glm::ivec3 voxel = pos_to_voxel(pos);
    voxel.y = 0;
    float res = 0;
    if (in_voxel_cube(voxel))
    {
        for (int j=vox_y;j>=-vox_y;j--)
        {
            voxel.y = j;
            float o = voxels[v_to_i(voxel)];
            if (o < 5*1e9)
            {
                res += o;
            }
            else //reached heightmap
            {
                break;
            }
        }
        //logerr("occ projection %f %f %f - %f",pos.x,pos.y,pos.z,res);
        return res;
    }
    else
    {
        return 1e9; 
    }
}
float LightVoxelsCube::get_occlusion_simple_mip(glm::vec3 pos, int mip)
{
    glm::ivec3 voxel = pos_to_voxel(pos);
    mip = CLAMP(mip, 0, mip_levels);
    if (in_voxel_cube(voxel))
    {
        int v = v_to_i_mip(voxel, mip);
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
    //logerr("pos voxel pos to voxel %f %f %f %d %d %d %d %d %d",pos.x, pos.y, pos.z, x,y,z, vox_x, vox_y, vox_z);
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

float LightVoxelsCube::get_occlusion_trilinear_mip(glm::vec3 pos, int mip)
{
    glm::vec3 flpos = glm::vec3(floorf(pos.x),floorf(pos.y),floorf(pos.z));
    glm::vec3 dp = pos - flpos;
    glm::ivec3 voxel = pos_to_voxel(flpos);
    int x = voxel.x;
    int y = voxel.y;
    int z = voxel.z;
    mip = CLAMP(mip, 0, mip_levels);
    #define CM(i,j,k) voxels[v_to_i_mip(x+i,y+j,z+k,mip)]
    if (in_voxel_cube(voxel) && in_voxel_cube(voxel + glm::ivec3(1,1,1)))
    {
        float c00 = CM(0,0,0)*(1-dp.x) + CM(1,0,0)*dp.x;
        float c01 = CM(0,0,1)*(1-dp.x) + CM(1,0,1)*dp.x;
        float c10 = CM(0,1,0)*(1-dp.x) + CM(1,1,0)*dp.x;
        float c11 = CM(0,1,1)*(1-dp.x) + CM(1,1,1)*dp.x;

        float c0 = c00*(1-dp.y) + c10*dp.y;
        float c1 = c01*(1-dp.y) + c11*dp.y;

        float c = c0*(1-dp.z) + c1*dp.z;

        return c;
    }
    else if (in_voxel_cube(voxel))
    {
        return CM(0,0,0);
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
    return MSE/SQ;
}
void LightVoxelsCube::add_heightmap(Heightmap &h)
{
    float thr = voxel_size + 0.01;
    for (int i=-vox_x;i<=vox_x;i++)
    {
        for (int k=-vox_z;k<=vox_z;k++)
        {
            glm::vec3 terr_pos = center + voxel_size*glm::vec3(i,0,k);
            terr_pos.y = h.get_height_simple(terr_pos);
            for (int j=-vox_y;j<=vox_y;j++)
            {
                glm::vec3 pos = center + voxel_size*glm::vec3(i,j,k);
                if (pos.y < terr_pos.y - thr)
                {
                    voxels[v_to_i(i,j,k)] = MAX(1e10, voxels[v_to_i(i,j,k)]);
                }
                else
                {
                    break;
                }
            }
        }
    }
}
void LightVoxelsCube::add_voxels_cube(LightVoxelsCube *cube)
{
    if (!cube)
        return;
    AABB box1 = get_bbox();
    AABB box2 = cube->get_bbox();
    glm::vec3 mn = max(box1.min_pos,box2.min_pos);
    glm::vec3 mx = min(box1.max_pos,box2.max_pos);
    if (mn.x >= mx.x || mn.y >= mx.y || mn.z >= mx.z)
        return;
    glm::ivec3 vox_mn = pos_to_voxel(mn);
    glm::ivec3 vox_mx = pos_to_voxel(mx);
    for (int i=vox_mn.x;i<=vox_mx.x;i++)
    {
        for (int j=vox_mn.y;j<=vox_mx.y;j++)
        {
            for (int k=vox_mn.z;k<=vox_mx.z;k++)
            {
                glm::vec3 pos = center + voxel_size*glm::vec3(i,j,k);
                voxels[v_to_i(i,j,k)] = MAX(cube->get_occlusion_trilinear(pos), voxels[v_to_i(i,j,k)]);
            }
        }
    }
}
void LightVoxelsCube::add_body(Body *b, float opacity, bool solid, float infl_distance, float base_infl_occ)
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    if (!b)
        return;
    Box *box = dynamic_cast<Box *>(b);
    if (box)
    {
        vec3 max_pos = box->get_Bbox().position + box->get_Bbox().sizes;
        AABB bbox = AABB(box->get_Bbox().position, max_pos);
        add_AABB(bbox, false, opacity);
    }
    else
    {
        int infl_dist = infl_distance/voxel_size;

        for (int i=-vox_x;i<=vox_x;i++)
        {
            for (int j=-vox_y;j<=vox_y;j++)
            {
                for (int k=-vox_z;k<=vox_z;k++)
                {
                    glm::vec3 pos = center + voxel_size*glm::vec3(i,j,k);
                    if (b->in_body(pos))   
                    {
                        if (solid)
                            voxels[v_to_i(i,j,k)] += 1e7;
                        if (opacity > 0.1)
                        {
                            set_occluder_simple(pos,opacity);

                            if (infl_dist > 0)
                            {
                                if (!b->in_body(pos + glm::vec3(voxel_size, 0, 0)))
                                {
                                    for (int l = 1; l<=infl_dist;l++)
                                        set_occluder_simple(pos + glm::vec3(l*voxel_size, 0, 0),
                                                            (base_infl_occ*(infl_dist-l))/infl_dist);
                                }
                                if (!b->in_body(pos + glm::vec3(-voxel_size, 0, 0)))
                                {
                                    for (int l = 1; l<=infl_dist;l++)
                                        set_occluder_simple(pos + glm::vec3(-l*voxel_size, 0, 0),
                                                            (base_infl_occ*(infl_dist-l))/infl_dist);
                                }
                                if (!b->in_body(pos + glm::vec3(0, 0, voxel_size)))
                                {
                                    for (int l = 1; l<=infl_dist;l++)
                                        set_occluder_simple(pos + glm::vec3(0, 0, l*voxel_size),
                                                            (base_infl_occ*(infl_dist-l))/infl_dist);
                                }
                                if (!b->in_body(pos + glm::vec3(0, 0, -voxel_size)))
                                {
                                    for (int l = 1; l<=infl_dist;l++)
                                        set_occluder_simple(pos + glm::vec3(0, 0, -l*voxel_size),
                                                            (base_infl_occ*(infl_dist-l))/infl_dist);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logerr("body added took %.2f ms", 0.001*time);
}

void LightVoxelsCube::add_AABB(const AABB &box, bool precise, float opacity)
{
    if (!precise)
    {
        ivec3 st_b = clamp(pos_to_block(box.min_pos), ivec3(0,0,0), ivec3(block_x-1, block_y-1, block_z-1));
        ivec3 end_b = clamp(pos_to_block(box.max_pos), ivec3(0,0,0), ivec3(block_x-1, block_y-1, block_z-1));

        fill_blocks(st_b, end_b, opacity);
    }
    else
    {
        ivec3 st_v = clamp(pos_to_voxel(box.min_pos), ivec3(-vox_x,-vox_y,-vox_z), ivec3(vox_x, vox_x, vox_z));
        ivec3 end_v = clamp(pos_to_voxel(box.max_pos), ivec3(-vox_x,-vox_y,-vox_z), ivec3(vox_x, vox_x, vox_z));

        fill_voxels(st_v, end_v, opacity);        
    }
}

void LightVoxelsCube::calculte_precise_occlusion_from_bodies()
{
    glm::vec3 dv = glm::vec3(voxel_size,voxel_size,voxel_size);
    std::vector<glm::vec3> sun_dirs;
    int steps = lightParams.sunPositions;
    steps = 6;
    logerr("body occlusion %d steps", steps);
    if (steps < 3)
        return;
    for (int i=0;i<steps;i++)
    {
        EnvironmentParameters params;
        params.hours = 6 + (12.0*i)/steps;
        glm::vec3 sun = Sun::sun_direction(params);
        if (sun.y < 0) //sun above horizon
            sun_dirs.push_back(-sun);   
    }
    if (sun_dirs.size() < 3)
        return;
    std::vector<glm::vec3> main_dirs = {sun_dirs.front(),sun_dirs[sun_dirs.size()/2],sun_dirs.back()};
    logerr("body occlusion %d dirs", sun_dirs.size());
    for (int i=-vox_x;i<=vox_x;i++)
    {
        for (int j=-vox_y;j<=vox_y;j++)
        {
            for (int k=-vox_z;k<=vox_z;k++)
            {
                glm::vec3 pos = center + voxel_size*glm::vec3(i,j,k);
                float steps = 0;
                float occs = 0;
                float max_occ = 50;
                for (glm::vec3 &dir : main_dirs)
                {
                    steps+=dir.y;
                    float occ = get_occlusion_view_ray(pos,dir);
                    if (occ > 1e6)//body intersected
                        occs+=dir.y;
                }
                if (occs > 0)
                {
                    for (glm::vec3 &dir : sun_dirs)
                    {
                        steps+=dir.y;
                        float occ = get_occlusion_view_ray(pos,dir);
                        if (occ > 1e6)//body intersected
                            occs+=dir.y;
                    }
                    set_occluder_simple(pos,max_occ*occs/(steps));
                }
            }
        }
    }
}
glm::vec3 LightVoxelsCube::get_dir_to_bright_place_ray(glm::vec3 pos, float length, int rays_sqrt, float *occlusion)
{
    const float BIAS = 0.01;
    float step = 2*PI/rays_sqrt;
    float r_step = 0.75*voxel_size;
    float *occlusions = new float[rays_sqrt*rays_sqrt];
    for (int i=0;i<rays_sqrt*rays_sqrt;i++)
    {
        occlusions[i] = 0;
    }
    for (float r = r_step; r < length; r+= r_step)
    {
        for (float phi=0;phi<2*PI - step;phi += step)
        {
            for (float psi=0;psi<2*PI;psi += step)
            {
                float x = pos.x + r*cosf(phi)*cosf(psi);
                float y = pos.y + r*sinf(phi)*cosf(psi);
                float z = pos.z + r*sinf(psi);
                occlusions[((int)(phi/step))*rays_sqrt + (int)(psi/step)] += get_occlusion_trilinear(glm::vec3(x,y,z))/(r*r);
            }
        }
    }
    std::vector<int> mps;
    float mo = occlusions[0];
    for (int i=0;i<rays_sqrt*rays_sqrt;i++)
    {
        if (occlusions[i] < mo - BIAS)
        {
            mo = occlusions[i];
            mps.clear();
            mps.push_back(i);
        }
        else if (abs(occlusions[i] - mo) < BIAS)
            mps.push_back(i);
    }
    int min_sz = mps.size();
    int mp = mps[urandi(0,min_sz)];
    float phi = step * (mp / rays_sqrt);
    float psi = step * (mp % rays_sqrt);
    float x = cosf(phi)*cosf(psi);
    float y = sinf(phi)*cosf(psi);
    float z = sinf(psi);

    *occlusion = mo;
    if (occlusions)
        delete[] occlusions;
    return glm::vec3(x,y,z);
}

void LightVoxelsCube::prepare_mip_levels()
{
    int p = mip_decrease;
    for (int i=1;i<mip_levels;i++)
    {
        logerr("prepare mip level %d %d %d %d %d %d",
        mip_vox_xyz[i].x,mip_vox_xyz[i].y, mip_vox_xyz[i].z,
        mip_vox_xyz[i-1].x,mip_vox_xyz[i-1].y,mip_vox_xyz[i-1].z);
        for (int x = 0;x<mip_vox_xyz[i].x;x++)
        {
            for (int y = 0;y<mip_vox_xyz[i].y;y++)
            {
                for (int z = 0;z<mip_vox_xyz[i].z;z++)
                {
                    if (i == 1)
                    {
                        int a = v_to_i_mip_no_offset(x,y,z,i-1);
                        int b = v_to_i(x-vox_x,y-vox_y,z-vox_z);
                        if (a!=b)
                        logerr("%d %d",a,b);
                    }
                    int x1 = MIN(p*x+1,mip_vox_xyz[i-1].x - 1);
                    int y1 = MIN(p*y+1,mip_vox_xyz[i-1].y - 1);
                    int z1 = MIN(p*z+1,mip_vox_xyz[i-1].z - 1);

                    voxels[v_to_i_mip_no_offset(x,y,z,i)] = 
                    0.125*voxels[v_to_i_mip_no_offset(p*x  ,p*y  ,p*z  ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(p*x  ,p*y  ,z1   ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(p*x  ,y1   ,p*z  ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(p*x  ,y1   ,z1   ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(x1   ,p*y  ,p*z  ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(x1   ,p*y  ,z1   ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(x1   ,y1   ,p*z  ,i-1)] +
                    0.125*voxels[v_to_i_mip_no_offset(x1   ,y1   ,z1   ,i-1)];
                    debug("%.1f ",voxels[v_to_i_mip_no_offset(x,y,z,i)]);
                }
                debug("\n");
            }
            debug("\n");
            debug("\n");
        }
    }
}

void LightVoxelsCube::read_func(const std::function<void(glm::vec3 &, float )> reader)
{
    for (int i=-vox_x;i<=vox_x;i++)
    {
        for (int j=-vox_y;j<=vox_y;j++)
        {
            for (int k=-vox_z;k<=vox_z;k++)
            {
                glm::vec3 pos = center + voxel_size*glm::vec3(i,j,k);
                reader(pos, voxels[v_to_i(i,j,k)]);
            }
        }
    }
}

void LightVoxelsCube::read_func_simple(const std::function<void(float )> reader)
{
    for (int i=0;i<count;i++)
        reader(voxels[i]);
}