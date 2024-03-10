#pragma once
#include "common_utils/LiteMath_ext.h"
#include <list>
#include "common_utils/body.h"
#include "graphics_utils/terrain.h"
#include "common_utils/distribution.h"
#include "tinyEngine/model.h"
#include "tinyEngine/shader.h"
#include <functional>
#include <vector>
struct LightVoxelsCube
{
public:
    LightVoxelsCube(float3 center, float3 size, float base_size, int mip_levels = 1, int mip_decrease = 2, int preferred_block_size=-1);
    LightVoxelsCube(float3 center, int3 sizes, float voxel_size, int mip_levels = 1, int mip_decrease = 2, int preferred_block_size=-1);
    LightVoxelsCube(LightVoxelsCube *source);
    LightVoxelsCube(LightVoxelsCube *source, float3 pos, float3 sizes, int size_decrease = 1, 
                    float2 min_max = float2(0,1e10));
    LightVoxelsCube(LightVoxelsCube *source, int3 vox_pos, int3 vox_sizes, int size_decrease,
                    float2 min_max);
    LightVoxelsCube() {};
    ~LightVoxelsCube();
    void clear();
    void fill(float val);
    void clamp_values(float min = -1, float max = -1);
    void set_occluder(float3 pos, float strenght);
    void set_occluder_pyramid_fast(float3 pos, float strenght, int max_r, int rnd_seed);
    void set_occluder_simple(float3 pos, float strenght);
    void set_occluder_simple_mip(float3 pos, float strenght, int mip);
    void set_occluder_trilinear(float3 pos, float strenght);
    float get_occlusion_voxel(int3 voxel);
    float get_occlusion_voxel_unsafe(int x, int y, int z);
    float get_occlusion(float3 pos);
    float get_occlusion_simple(float3 pos);
    float get_occlusion_projection(float3 pos);
    float get_occlusion_simple_mip(float3 pos, int mip);
    float get_occlusion_trilinear(float3 pos);
    float get_occlusion_trilinear_mip(float3 pos, int mip);
    void prepare_mip_levels();
    void add_body(Body *b, float opacity = 1e9, bool solid = true, float infl_distance = 0, float base_infl_occ = 0);
    void add_AABB(const AABB &box, bool precise, float opacity = 1e9);
    void add_heightmap(Heightmap &h);
    void add_voxels_cube(LightVoxelsCube *cube, bool fast_fill_expected = false);
    void get_data(float **data, int3 &size);
    AABB get_bbox();
    float NMSE(LightVoxelsCube *B);
    float3 get_center();
    float get_voxel_size();
    int3 get_vox_sizes();
    int get_mip_count();
    int get_mip_decrease();
    long get_size_cnt() {return count;}
    void read_func(const std::function<void(float3 &, float )> reader);
    void read_func_simple(const std::function<void(float )> reader);
    void relocate(float3 new_pos) { center = new_pos; }
    int get_block_size() { return block_size; }
    bool empty() {return (voxels == nullptr);}
    static int get_default_block_size() { return 5; }

    float *voxels = nullptr;

    float3 center;
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
    std::vector<int3> mip_vox_xyz;

private:

    float fast_voxel_occlusion(int3 voxel);
    float voxel_occlusion(int3 voxel);
    int3 pos_to_voxel(float3 pos);
    int3 pos_to_block(float3 pos);
    float3 voxel_to_pos(int3 voxel);
    void fill_blocks(int3 from, int3 to, float val);
    void fill_voxels(int3 from, int3 to, float val);
    bool in_voxel_cube(int3 voxel);
    bool in_voxel_cube(int x, int y, int z);
    int v_to_i(int3 voxel);
    int v_to_i(int x, int y, int z);
    int v_to_i_mip(int3 voxel, int mip);
    int v_to_i_mip(int x, int y, int z, int mip);
    int v_to_i_mip_no_offset(int3 voxel, int mip);
    int v_to_i_mip_no_offset(int x, int y, int z, int mip);
    int voxels_count() {return vox_x*vox_y*vox_z;}
    void set_occluder_pyramid_tail_5(int3 voxel, float strenght, int max_r, int rnd_seed);
    void set_occluder_pyramid_head_5(int3 voxel, float strenght, int max_r);
};