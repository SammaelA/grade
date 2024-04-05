#pragma once
#include <vector>
#include <functional>
#include <fstream>
#include "common_utils/LiteMath_ext.h"
#include "LiteRT/sdfScene/sdf_scene.h"
/*
A class that is able to represent and arbitrary function f : R^3 -> T 
as a sparse octree, where every leaf contains value of function in it's
center. 
T - is a type that can be meaningfully interpolated (i.e. float, double, float2, float3)
It should have T+T and float*T operator and be POD of course. 
Octree always represents unit cube [-1,1]^3
*/

struct SparseOctreeSettings
{
  unsigned max_nodes = 1000000;
  unsigned depth = 1;
  unsigned min_remove_level = 4;
  float remove_thr = 0.0001;
};

template <typename T>
struct BlockSparseOctree
{
  //Page size in Vulkan sparse texture
  static constexpr unsigned BLOCK_SIZE_X = 32;
  static constexpr unsigned BLOCK_SIZE_Y = 32;
  static constexpr unsigned BLOCK_SIZE_Z = 16;
  struct BlockInfo
  {
    uint3 coords; //offset in blocks inside mip
    unsigned mip;
    unsigned data_offset; //offset in T array (not bytes!). All block data stored together starting with this offset
  };

  uint3 block_size = uint3(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
  uint3 top_mip_size = uint3(1,1,2); //to create cube grid from 32*32*16 blocks
  unsigned top_mip;

  std::vector<BlockInfo> blocks;
  std::vector<T> data; //for each block values is stored in z y x order to comply with Vulkan sparse texture
};

template <typename T>
void save_block_sparse_octree(const std::string &path, const BlockSparseOctree<T> &bso)
{
  std::ofstream fs(path, std::ios::binary);
  std::vector<unsigned> metadata(9);
  metadata[0] = bso.block_size.x; 
  metadata[1] = bso.block_size.y;
  metadata[2] = bso.block_size.z;

  metadata[3] = bso.top_mip_size.x; 
  metadata[4] = bso.top_mip_size.y;
  metadata[5] = bso.top_mip_size.z;  

  metadata[6] = bso.top_mip; 
  metadata[7] = bso.blocks.size();
  metadata[8] = bso.data.size();

  fs.write((const char *)metadata.data(), metadata.size()*sizeof(unsigned));
  fs.write((const char *)bso.blocks.data(), bso.blocks.size() * sizeof(typename BlockSparseOctree<T>::BlockInfo));
  fs.write((const char *)bso.data.data(), bso.data.size() * sizeof(T));
  fs.flush();
  fs.close();
}

template <typename T>
void load_block_sparse_octree(const std::string &path, /*out*/BlockSparseOctree<T> &bso)
{
  std::ifstream fs(path, std::ios::binary);
  std::vector<unsigned> metadata(9);

  fs.read((char *)metadata.data(), metadata.size()*sizeof(unsigned));

  bso.block_size.x = metadata[0];
  bso.block_size.y = metadata[1];
  bso.block_size.z = metadata[2];

  bso.top_mip_size.x = metadata[3];
  bso.top_mip_size.y = metadata[4];
  bso.top_mip_size.z = metadata[5];

  bso.top_mip = metadata[6];
  bso.blocks.resize(metadata[7]);
  bso.data.resize(metadata[8]);

  fs.read((char *)bso.blocks.data(), bso.blocks.size() * sizeof(typename BlockSparseOctree<T>::BlockInfo));
  fs.read((char *)bso.data.data(), bso.data.size() * sizeof(T));

  fs.close();
}

class SparseOctreeBuilder
{
public:
  using Node = SdfOctreeNode;
  using T = float;

  static bool is_border(float distance, unsigned level);

  SparseOctreeBuilder();
  void construct_bottom_up(std::function<T(const float3 &)> f, SparseOctreeSettings settings);
  void construct_bottom_up_blocks(std::function<T(const float3 &)> f, SparseOctreeSettings settings, 
                                  BlockSparseOctree<T> &out_bso);
  T sample(const float3 &pos, unsigned max_level = 1000) const;
  T sample_closest(const float3 &pos) const;

  void print_stat() const;
  std::pair<float,float> estimate_quality(std::function<T(const float3 &)> reference_f, float dist_thr = 0.01, unsigned samples = 10000) const;
  std::vector<Node> &get_nodes() { return octree_f->get_nodes(); }
  const std::vector<Node> &get_nodes() const { return octree_f->get_nodes(); }

  Node &get_node(unsigned idx) { return octree_f->get_nodes()[idx]; }
  const Node &get_node(unsigned idx) const { return octree_f->get_nodes()[idx]; }
protected:
  void add_node_rec(std::function<T(const float3 &)> f, unsigned node_idx, unsigned depth,
                    unsigned max_depth, float3 p, float d);
  void split_children(std::function<T(const float3 &)> f, unsigned node_idx,
                      float threshold, float3 p, float d, unsigned level);

  void construct_bottom_up_base(std::function<T(const float3 &)> f, SparseOctreeSettings settings);
  void construct_bottom_up_finish(std::function<T(const float3 &)> f, SparseOctreeSettings settings);

  std::shared_ptr<ISdfOctreeFunction> octree_f; //0 node is root
};