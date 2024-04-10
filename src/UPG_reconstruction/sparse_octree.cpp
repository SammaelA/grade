#include "sparse_octree.h"
#include "common_utils/distribution.h"
#include "common_utils/stat_box.h"
#include <cassert>
#include <map>

SparseOctreeBuilder::SparseOctreeBuilder()
{
  std::vector<Node> nodes;
  octree_f = get_SdfOctreeFunction({0, nodes.data()});
}

bool SparseOctreeBuilder::is_border(float distance, int level)
{
  return level < 2  ? true : abs(distance) < sqrt(2)*pow(2, -level);
}

constexpr unsigned INVALID_IDX = 1u<<31u;
bool is_leaf(unsigned offset)
{
  return (offset == 0) || (offset & INVALID_IDX);
}

void SparseOctreeBuilder::add_node_rec(unsigned node_idx, unsigned depth, unsigned max_depth, float3 p, float d)
{
  auto &nodes = get_nodes();
  nodes[node_idx].value = sdf(2.0f * ((p + float3(0.5, 0.5, 0.5)) * d) - float3(1, 1, 1));

  if (depth < max_depth && is_border(nodes[node_idx].value, depth))
  {
    nodes[node_idx].offset = nodes.size();
    
    nodes.resize(nodes.size() + 8);
    unsigned idx = nodes[node_idx].offset;
    for (unsigned cid = 0; cid < 8; cid++)
    {
      add_node_rec(idx + cid, depth + 1, max_depth, 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1), d / 2);
    }
  }
}

bool DBG = false;

void SparseOctreeBuilder::split_children(unsigned node_idx, float threshold, float3 p, float d, unsigned level)
{
  auto &nodes = get_nodes();
  assert(!is_leaf(nodes[node_idx].offset));
  unsigned idx = nodes[node_idx].offset;
  T n_distances[8];
  for (unsigned cid = 0; cid < 8; cid++)
    n_distances[cid] = nodes[idx + cid].value;
  for (unsigned cid = 0; cid < 8; cid++)
  {
    if (!is_leaf(nodes[idx + cid].offset)) // go deeper
    {
      split_children(idx + cid, threshold, 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1), d/2, level+1);
    }
    else if (is_border(abs(nodes[idx + cid].value), level+1)) // child is leaf, check if we should split it
    {
      //printf("LOL %u %u %u\n",(cid & 4) >> 2, (cid & 2) >> 1, cid & 1);
      float3 p1 = 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1);
      float3 p2 = (p1+float3(0.5,0.5,0.5)) * (d/2);
      float d1 = sdf(2.0f*p2-1.0f);
      float d2 = sample_closest(2.0f*p2-1.0f);
      if (abs(d1-d2)>1e-4)
        printf("ERRRR %f %f %f --- %f %f\n",p2.x, p2.y, p2.z, d1,d2);

      bool need_split = false;
      unsigned samples = 64;
      float av_diff = 0;
      for (unsigned s=0;s<samples;s++)
      {
        float3 n_pos = 0.5f*float3(0.5f + ((s & 4) >> 2), 0.5f + ((s & 2) >> 1), 0.5f + (s & 1));
        float3 pos = (p1 + n_pos) * (d/2);
        pos = (p1 + float3(urand(), urand(), urand())) * (d/2);
        float d_ref = sdf(2.0f*pos-1.0f);
        float d_sample = sample(2.0f*pos-1.0f);
        float diff = abs(d_ref - d_sample);
        av_diff += diff;
        //if (diff > threshold)
        //  printf("%f %f %f %f %f %f diff = %f (%f %f)\n", p1.x, p1.y, p1.z, pos.x, pos.y, pos.z, diff, d_ref, d_sample);
      }
      av_diff /= samples;
      if (av_diff > threshold)
      {
        printf("Performing split. Level %u size %d\n", level, (int)nodes.size());
        need_split = true;
      }

      if (need_split && level<10)
      {
        nodes[idx + cid].offset = nodes.size();
        
        nodes.resize(nodes.size() + 8);
        unsigned gc_idx = nodes[idx + cid].offset;
        for (unsigned gcid = 0; gcid < 8; gcid++)
        {
          float3 pos = (p1 + 0.5f*float3(0.5f + ((gcid & 4) >> 2), 0.5f + ((gcid & 2) >> 1), 0.5f + (gcid & 1))) * (d/2);
          nodes[gc_idx + gcid].value = sdf(2.0f*pos-1.0f);
        }        
      }
    }
  }
}

SparseOctreeBuilder::T SparseOctreeBuilder::sample(const float3 &position, unsigned max_level) const
{
  return octree_f->eval_distance_level(position, max_level);
}

SparseOctreeBuilder::T SparseOctreeBuilder::sample_closest(const float3 &position) const
{
  auto &nodes = get_nodes();
  //if (abs(position).x > 0.9 || abs(position).y > 0.9 || abs(position).z > 0.9)
  //  return 0.1;
  float3 pos = LiteMath::clamp(0.5f*(position + 1.0f), 0.0f, 1.0f);
  unsigned idx = 0;
  float d = 1;
  float3 p = float3(0,0,0);
  while (nodes[idx].offset != 0)
  {
    float3 pindf = pos/d - p;
    unsigned ch_index = 4*(pindf.x >= 0.5) + 2*(pindf.y >= 0.5) + (pindf.z >= 0.5);
    //printf("%u pindf %f %f %f %u\n",idx, pindf.x, pindf.y, pindf.z, ch_index);
    idx = nodes[idx].offset + ch_index;
    d = d/2;
    p = 2*p + float3((ch_index & 4) >> 2, (ch_index & 2) >> 1, ch_index & 1);
  }
  //printf("\n");
  //printf("%u last pindf \n",idx);

  return nodes[idx].value;
}

struct SparseOctreeCounts
{
  std::vector<unsigned> count_all;
  std::vector<unsigned> count_border;
  std::vector<unsigned> count_leaf;
  std::vector<unsigned> count_border_leaf;
};

void check_border_reappearance_rec(const std::vector<SparseOctreeBuilder::Node> &nodes, unsigned idx, unsigned level)
{
  bool b = SparseOctreeBuilder::is_border(nodes[idx].value, level); 
  if (!is_leaf(nodes[idx].offset))
  {
    for (unsigned ch_idx=0; ch_idx<8; ch_idx++)
    {
      bool ch_b = SparseOctreeBuilder::is_border(nodes[nodes[idx].offset + ch_idx].value, level+1);
      if (!b && ch_b)
        printf("reappeared border (level %d-%d): %f %f\n", level,level+1, nodes[idx].value, nodes[nodes[idx].offset + ch_idx].value);
      check_border_reappearance_rec(nodes, nodes[idx].offset + ch_idx, level+1);
    }
  }
}

void print_stat_rec(const std::vector<SparseOctreeBuilder::Node> &nodes, SparseOctreeCounts &counts, unsigned idx, unsigned level)
{
  if (level == counts.count_all.size())
  {
    counts.count_all.push_back(0);
    counts.count_border.push_back(0);
    counts.count_leaf.push_back(0);
    counts.count_border_leaf.push_back(0);
  }
  counts.count_all[level] += 1;
  counts.count_border[level] += SparseOctreeBuilder::is_border(nodes[idx].value, level);
  counts.count_leaf[level] += is_leaf(nodes[idx].offset);
  counts.count_border_leaf[level] += SparseOctreeBuilder::is_border(nodes[idx].value, level) && is_leaf(nodes[idx].offset);
  if (!is_leaf(nodes[idx].offset))
  {
    for (unsigned ch_idx=0; ch_idx<8; ch_idx++)
      print_stat_rec(nodes, counts, nodes[idx].offset + ch_idx, level+1);
  }
}

void SparseOctreeBuilder::print_stat() const
{
  auto &nodes = get_nodes();
  check_border_reappearance_rec(nodes, 0, 0);

  SparseOctreeCounts counts;
  print_stat_rec(nodes, counts, 0, 0);
  printf("SparseOctreeBuilder::print_stat\n");
  printf("Levels: %d\n", (int)counts.count_all.size());
  float tmax = 1;
  for (int i=0;i<counts.count_all.size();i++)
  {
    float m = 100.0f/counts.count_all[i];
    printf("Level %2d: cnt %6u, b %5u (%5.1f%%), l %5u (%5.1f%%), bl %5u (%5.1f%%)  occ %5.1f%%\n", i, counts.count_all[i], 
           counts.count_border[i], m*counts.count_border[i], 
           counts.count_leaf[i], m*counts.count_leaf[i],
           counts.count_border_leaf[i], m*counts.count_border_leaf[i],
           100.0f*counts.count_all[i]/tmax);
    tmax *= 8;
  }
}

std::pair<float,float> SparseOctreeBuilder::estimate_quality(float dist_thr, unsigned samples) const
{
  stat::Bins<float> stat_box(0.0f, 0.1f);
  unsigned i = 0;
  std::vector<float> differences(samples, 0);
  while (i < samples)
  {
    float3 p = float3(urand(-1,1),urand(-1,1),urand(-1,1));
    float ref_d = sdf(p);
    if (abs(ref_d) < dist_thr)
    {
      differences[i] = sample(p) - ref_d;
      stat_box.add(abs(differences[i]));
      //if (abs(differences[i]) > 0.005)
      //{
      //  printf("AAAA %f %f %f, d = %f %f\n", p.x, p.y, p.z, sample(p), ref_d);
      //  DBG = true;
      //  sample(p);
      //  DBG = false;
      //}
      i++;
    }
  }

  double av_diff = 0;
  double max_diff = 0;
  for (auto &d : differences)
  {
    av_diff += abs(d);
    max_diff = std::max(max_diff, (double)d);
  }
  av_diff /= samples;
  //stat_box.print_bins();
  return {av_diff, max_diff};
}

float2 invalidate_node_rec(std::vector<SparseOctreeBuilder::Node> &nodes, unsigned idx)
{
  unsigned ofs = nodes[idx].offset;
  if (!is_leaf(ofs)) 
  {
    float min_val = FLT_MAX;
    float avg_val = 0;
    for (int i=0;i<8;i++)
    {
      float2 min_avg = invalidate_node_rec(nodes, ofs + i);
      min_val = std::min(min_val, min_avg.x);
      avg_val += min_avg.y;
      nodes[ofs + i].offset |= INVALID_IDX;
    }
    return float2(min_val, avg_val/8);
  }
  else
    return float2(nodes[idx].value);
}

void remove_linear_rec(SparseOctreeBuilder &octree, float thr, unsigned min_level_to_remove, unsigned idx, unsigned level, float3 p, float d)
{
  unsigned ofs = octree.get_nodes()[idx].offset;
  if (is_leaf(ofs)) 
    return;
  
  bool diff_less = true;
  if (level >= min_level_to_remove)
  {
    for (int i=0;i<8;i++)
    {
      float ch_d = d/2;
      float3 ch_p = 2*p + float3((i & 4) >> 2, (i & 2) >> 1, i & 1);
      float dist = octree.get_nodes()[ofs + i].value;
      float sampled_dist = octree.sample(2.0f*(ch_d*(ch_p + float3(0.5,0.5,0.5))) - 1.0f, level-1);
      //printf("d sd %f %f\n",dist, sampled_dist);
      if (abs(dist - sampled_dist) > thr)
      {
        diff_less = false;
        break;
      }
    }  
  }
  else
    diff_less = false;

  if (!diff_less)
  {
    for (int i=0;i<8;i++)
    {
      float ch_d = d/2;
      float3 ch_p = 2*p + float3((i & 4) >> 2, (i & 2) >> 1, i & 1);
      remove_linear_rec(octree, thr, min_level_to_remove, ofs + i, level+1, ch_p, ch_d);
    }
  }
  else
  {
    octree.get_nodes()[idx].value = invalidate_node_rec(octree.get_nodes(), idx).y;
    octree.get_nodes()[idx].offset |= INVALID_IDX;
  }
}

void all_to_valid_nodes_remap_rec(const std::vector<SparseOctreeBuilder::Node> &all_nodes, 
                                  std::vector<int> &all_to_valid_remap,
                                  unsigned *valid_counter,
                                  unsigned idx)
{
  unsigned ch_idx = all_nodes[idx].offset;
  if (is_leaf(ch_idx))
    return;
  else
  {
    for (int i=0;i<8;i++)
      all_to_valid_remap[ch_idx+i] = (*valid_counter) + i;
    
    (*valid_counter) += 8;

    for (int i=0;i<8;i++)
      all_to_valid_nodes_remap_rec(all_nodes, all_to_valid_remap, valid_counter, ch_idx+i);
  }
}

void SparseOctreeBuilder::construct_bottom_up_finish()
{
  auto &nodes = get_nodes();
  //auto p = estimate_quality(10, 100000);
  //printf("estimate_quality: %f %f\n", (float)p.first, (float)p.second);

  std::vector<int> valid_remap(nodes.size(), -1);
  unsigned valid_count = 1;

  valid_remap[0] = 0; //root is always valid
  all_to_valid_nodes_remap_rec(nodes, valid_remap, &valid_count, 0);

  std::vector<Node> valid_nodes(valid_count);
  for (int i=0;i<nodes.size();i++)
  {
    if (valid_remap[i] >= 0)
    {
      valid_nodes[valid_remap[i]] = nodes[i];
      if (is_leaf(valid_nodes[valid_remap[i]].offset))
        valid_nodes[valid_remap[i]].offset = 0;
      else
        valid_nodes[valid_remap[i]].offset = valid_remap[valid_nodes[valid_remap[i]].offset];
    }
  }
  
  nodes = valid_nodes;

  //printf("SDF octree created with %u (%5.1f%%) nodes (dense one would have %u)\n", 
  //       valid_count, 100.0f*valid_count/(unsigned)valid_remap.size(), (unsigned)valid_remap.size());
}

void SparseOctreeBuilder::construct_bottom_up_base(unsigned start_depth, float3 start_p, float start_d)
{
  auto &nodes = get_nodes();
  nodes.clear();

  // create octree on border nodes
  nodes.emplace_back();
  add_node_rec(0, start_depth, settings.depth, start_p, start_d);

  //remove some nodes with estimated error is lower than given threshold
  if (settings.remove_thr > 0)
    remove_linear_rec(*this, settings.remove_thr, settings.min_remove_level, 0, start_depth, start_p, start_d);  
}

void SparseOctreeBuilder::construct_bottom_up(std::function<T(const float3 &)> _sdf, SparseOctreeSettings _settings)
{
  sdf = _sdf;
  settings = _settings;
  octree_f->get_nodes().clear();

  construct_bottom_up_base(0u, float3(0,0,0), 1.0f);
  construct_bottom_up_finish();
}

struct cmpUint3 {
    bool operator()(const uint3& a, const uint3& b) const 
    {
      if (a.x < b.x)
        return true;
      else if (a.x > b.x)
        return false;
      if (a.y < b.y)
        return true;
      else if (a.y > b.y)
        return false;
      return a.z < b.z;
      
    }
};

void add_all_blocks_rec(std::vector<std::vector<BlockSparseOctree<float>::BlockInfo>> &all_blocks, 
                        std::vector<std::vector<bool>> &block_active, 
                        std::vector<std::map<uint3, unsigned, cmpUint3>> &block_idx_to_block_n, 
                        unsigned mip, unsigned idx)
{
  if (all_blocks[mip][idx].mip > 0)
  {
    for (int i=0;i<8;i++)
    {
      uint3 off = uint3((i & 4) >> 2, (i & 2) >> 1, i & 1); 
      uint3 pidx = 2*all_blocks[mip][idx].coords + off;
      all_blocks[mip-1].push_back({pidx, mip-1, 0});
      block_active[mip-1].push_back(false);
      block_idx_to_block_n[mip-1][pidx] = all_blocks[mip-1].size()-1;
      add_all_blocks_rec(all_blocks, block_active, block_idx_to_block_n, mip-1, all_blocks[mip-1].size()-1);
    }
  }
}

void find_active_blocks_rec(const SparseOctreeBuilder &octree, std::vector<std::vector<bool>> &block_active, 
                            std::vector<std::map<uint3, unsigned, cmpUint3>> &block_idx_to_block_n,
                            unsigned min_depth, unsigned max_block_mip, 
                            unsigned idx, unsigned level, uint3 pixel_idx)
{
  if ((octree.get_nodes()[idx].offset & INVALID_IDX) == 0)
  {
    if (level >= min_depth)
    {
      unsigned block_mip = max_block_mip - (level - min_depth);
      uint3 block_idx = pixel_idx/uint3(BlockSparseOctree<float>::BLOCK_SIZE_X, BlockSparseOctree<float>::BLOCK_SIZE_Y, BlockSparseOctree<float>::BLOCK_SIZE_Z);
      block_active[block_mip][block_idx_to_block_n[block_mip][block_idx]] = true;
    }
    if (!is_leaf(octree.get_nodes()[idx].offset))
    {
      for (int i=0;i<8;i++)
      {
        uint3 off = uint3((i & 4) >> 2, (i & 2) >> 1, i & 1); 
        find_active_blocks_rec(octree, block_active, block_idx_to_block_n, min_depth, max_block_mip, octree.get_nodes()[idx].offset+i, level+1, 2*pixel_idx+off);
      }
    }
  }
}

void restore_invalid_nodes_in_active_blocks_rec(SparseOctreeBuilder &octree, std::vector<std::vector<bool>> &block_active, 
                                                std::vector<std::map<uint3, unsigned, cmpUint3>> &block_idx_to_block_n,
                                                unsigned min_depth, unsigned max_block_mip, 
                                                unsigned idx, unsigned level, uint3 pixel_idx)
{
  unsigned real_offset = octree.get_nodes()[idx].offset & (~INVALID_IDX);
  //printf("offset %u %u \n", octree.get_nodes()[idx].offset, real_offset);

  if (level >= min_depth)
  {
    unsigned block_mip = max_block_mip - (level - min_depth);
    uint3 block_idx = pixel_idx/uint3(BlockSparseOctree<float>::BLOCK_SIZE_X, BlockSparseOctree<float>::BLOCK_SIZE_Y, BlockSparseOctree<float>::BLOCK_SIZE_Z);

    //invalid node in active block
    if ((octree.get_nodes()[idx].offset & INVALID_IDX) != 0 && block_active[block_mip][block_idx_to_block_n[block_mip][block_idx]])
    {
      //printf("reactivated node %u in block %u %u %u\n", idx, block_idx.x, block_idx.y, block_idx.z);
      octree.get_nodes()[idx].offset = real_offset;
    }
  }

  if (real_offset != 0)
  {
    for (int i=0;i<8;i++)
    {
      uint3 off = uint3((i & 4) >> 2, (i & 2) >> 1, i & 1); 
      restore_invalid_nodes_in_active_blocks_rec(octree, block_active, block_idx_to_block_n, min_depth, max_block_mip, real_offset+i, level+1, 2*pixel_idx+off);
    }    
  }
}

void node_values_to_blocks_rec(SparseOctreeBuilder &octree, std::vector<std::vector<bool>> &block_active, 
                               std::vector<std::vector<BlockSparseOctree<float>::BlockInfo>> &all_blocks,
                               std::vector<std::map<uint3, unsigned, cmpUint3>> &block_idx_to_block_n,
                               std::vector<float> &blocks_data,
                               unsigned min_depth, unsigned max_block_mip, 
                               unsigned idx, unsigned level, uint3 pixel_idx)
{
  unsigned real_offset = octree.get_nodes()[idx].offset & (~INVALID_IDX);

  if (level >= min_depth)
  {
    unsigned block_mip = max_block_mip - (level - min_depth);
    uint3 block_idx = pixel_idx/uint3(BlockSparseOctree<float>::BLOCK_SIZE_X, BlockSparseOctree<float>::BLOCK_SIZE_Y, BlockSparseOctree<float>::BLOCK_SIZE_Z);
    uint3 local_idx = pixel_idx % uint3(BlockSparseOctree<float>::BLOCK_SIZE_X, BlockSparseOctree<float>::BLOCK_SIZE_Y, BlockSparseOctree<float>::BLOCK_SIZE_Z);
    unsigned local_offset = local_idx.z * BlockSparseOctree<float>::BLOCK_SIZE_X * BlockSparseOctree<float>::BLOCK_SIZE_Y +
                            local_idx.y * BlockSparseOctree<float>::BLOCK_SIZE_X + 
                            local_idx.x;

    unsigned block_pos = block_idx_to_block_n[block_mip][block_idx];
    if (block_active[block_mip][block_pos])
      blocks_data[all_blocks[block_mip][block_pos].data_offset + local_offset] = octree.get_nodes()[idx].value;
  }

  if (real_offset != 0)
  {
    for (int i=0;i<8;i++)
    {
      uint3 off = uint3((i & 4) >> 2, (i & 2) >> 1, i & 1); 
      node_values_to_blocks_rec(octree, block_active, all_blocks, block_idx_to_block_n, blocks_data, min_depth, max_block_mip, 
                                real_offset+i, level+1, 2*pixel_idx+off);
    }    
  }
}

void SparseOctreeBuilder::construct_bottom_up_blocks(std::function<T(const float3 &)> _sdf, SparseOctreeSettings _settings, 
                                              BlockSparseOctree<T> &out_bso)
{
  sdf = _sdf;
  settings = _settings;
  octree_f->get_nodes().clear();

  unsigned max_block_size = std::max(BlockSparseOctree<float>::BLOCK_SIZE_X, std::max(BlockSparseOctree<float>::BLOCK_SIZE_Y, BlockSparseOctree<float>::BLOCK_SIZE_Z));
  unsigned max_depth_blocks = settings.depth - log2(max_block_size);

  //I - construct octree, leaving removed modes intact but with INVALID_IDX flad for their offsets
  construct_bottom_up_base(0u, float3(0,0,0), 1.0f);

  //II - restore invalid nodes in active blocks
  std::vector<std::vector<bool>> block_active;
  std::vector<std::vector<BlockSparseOctree<float>::BlockInfo>> all_blocks;
  std::vector<std::map<uint3, unsigned, cmpUint3>> block_idx_to_block_n;
  all_blocks.resize(max_depth_blocks + 1);
  block_active.resize(max_depth_blocks + 1);
  block_idx_to_block_n.resize(max_depth_blocks + 1);

  //convert SparseOctree to BlockSparseOctree
  out_bso.top_mip = max_depth_blocks;

  //create two root blocks and fill all_blocks with all possible blocks
  unsigned z_off = 0;
  while (z_off*BlockSparseOctree<float>::BLOCK_SIZE_Z != BlockSparseOctree<float>::BLOCK_SIZE_X)
  {
    all_blocks[max_depth_blocks].push_back({uint3(0,0,z_off), max_depth_blocks, 0}); 
    block_active[max_depth_blocks].push_back(false);
    block_idx_to_block_n[max_depth_blocks][uint3(0,0,z_off)] = z_off;
    add_all_blocks_rec(all_blocks, block_active, block_idx_to_block_n, max_depth_blocks, z_off);
    z_off++;
  }

  //check all nodes by block to find active blocks
  find_active_blocks_rec(*this, block_active, block_idx_to_block_n, log2(BlockSparseOctree<float>::BLOCK_SIZE_X), max_depth_blocks, 0, 0, uint3(0,0,0));

  //restore invalid nodes in active blocks
  restore_invalid_nodes_in_active_blocks_rec(*this, block_active, block_idx_to_block_n, 
                                             log2(BlockSparseOctree<float>::BLOCK_SIZE_X), max_depth_blocks, 
                                             0, 0, uint3(0,0,0));

  //calculate offsets for active blocks and allocate memory for it
  unsigned offset = 0;
  for (int i=0;i<max_depth_blocks + 1; i++)
  {
    for (int j=0;j<all_blocks[i].size();j++)
    {
      if (block_active[i][j])
      {
        all_blocks[i][j].data_offset = offset;
        offset += out_bso.block_size.x*out_bso.block_size.y*out_bso.block_size.z;
      }
    }
  }
  out_bso.data.resize(offset, T(FLT_MAX));

  //fill buffer for data per block
  node_values_to_blocks_rec(*this, block_active, all_blocks, block_idx_to_block_n, out_bso.data, 
                            log2(BlockSparseOctree<float>::BLOCK_SIZE_X), max_depth_blocks, 
                            0, 0, uint3(0,0,0));

  //collect all active 
  for (int i=0;i<max_depth_blocks + 1; i++)
  {
    for (int j=0;j<all_blocks[i].size();j++)
    {
      if (block_active[i][j])
        out_bso.blocks.push_back(all_blocks[i][j]);
    }
  }
  
  //DEBUG:  check is all the data for blocks is collected
  /*
  for (int i=0;i<max_depth_blocks + 1; i++)
  {
    for (int j=0;j<all_blocks[i].size();j++)
    {
      if (block_active[i][j])
      {
        for (int z=0;z<16;z++)
        {
          for (int y=0;y<32;y++)
          {
            for (int x=0;x<32;x++)
            {
              unsigned off = all_blocks[i][j].data_offset + z*32*32 + y*32 + x;
              float val = out_bso.data[off];
              if (val > 1e9)
              {
                printf("block [%u %u %u] in lod %u missed value in (%d %d %d)\n",
                       all_blocks[i][j].coords.x, all_blocks[i][j].coords.y, all_blocks[i][j].coords.z,
                       all_blocks[i][j].mip, x,y,z);
              }
            }
          }
        }
      }
    }
  }
  */

  //III - finish building octree
  construct_bottom_up_finish();
}

void fill_octree_frame_rec(std::function<SparseOctreeBuilder::T(const float3 &)> f,   
                           const std::vector<SdfOctreeNode> &nodes,
                           std::vector<SdfFrameOctreeNode> &frame,
                           unsigned idx, float3 p, float d)
{
  unsigned ofs = nodes[idx].offset;
  frame[idx].offset = ofs;
  if (is_leaf(ofs)) 
  {
    float3 pos = 2.0f*(d*p) - 1.0f;
    for (int i = 0; i < 8; i++)
    {
      float3 ch_pos = pos + 2*d*float3((i & 4) >> 2, (i & 2) >> 1, i & 1);
      frame[idx].values[i] = f(ch_pos);
      //frame[idx].values[i] = nodes[idx].value;
    }
  }
  else
  {
    for (int i = 0; i < 8; i++)
    {
      float ch_d = d / 2;
      float3 ch_p = 2 * p + float3((i & 4) >> 2, (i & 2) >> 1, i & 1);
      fill_octree_frame_rec(f, nodes, frame, ofs + i, ch_p, ch_d);
    }
  }
}

void SparseOctreeBuilder::construct_large_cell_rec(std::vector<Node> &final_nodes, unsigned root_idx, unsigned level, float3 p, float d)
{
  //construct octree in this large cell
  if (level == settings.min_remove_level)
  {

  }
}

//large node is a smallest unit for which octree can be built independently from other units
struct LargeNode
{
  float3 p;
  float d;
  unsigned level;
  unsigned root_idx;
  unsigned children_idx;
};

void SparseOctreeBuilder::construct(std::function<T(const float3 &)> _sdf, SparseOctreeSettings _settings)
{
  sdf = _sdf;
  settings = _settings;

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

  std::vector<LargeNode> large_nodes;
  large_nodes.push_back({float3(0,0,0), 1.0f, 0u, 0u, 0u});

  unsigned i = 0;
  while (i < large_nodes.size())
  {
    if (large_nodes[i].level < settings.min_remove_level)
    {
      large_nodes[i].children_idx = large_nodes.size();
      for (int j=0;j<8;j++)
      {
        float ch_d = large_nodes[i].d / 2;
        float3 ch_p = 2 * large_nodes[i].p + float3((j & 4) >> 2, (j & 2) >> 1, j & 1);
        large_nodes.push_back({ch_p, ch_d, large_nodes[i].level+1, i, 0u});
      }
    }
    i++;
  }

  std::vector<Node> all_nodes;
  all_nodes.resize(large_nodes.size());
  for (int i=0;i<large_nodes.size();i++)
  {
    all_nodes[i].offset = large_nodes[i].children_idx;
    all_nodes[i].value = sdf(2.0f * ((large_nodes[i].p + float3(0.5, 0.5, 0.5)) * large_nodes[i].d) - float3(1, 1, 1));
  }

  for (int i=0;i<large_nodes.size();i++)
  {
    //if it is a leaf (as a large node) and it can contain borders (with some additional confidence mult 1.05)
    if (large_nodes[i].children_idx == 0 && is_border(all_nodes[i].value, large_nodes[i].level))
    {
      construct_bottom_up_base(large_nodes[i].level, large_nodes[i].p, large_nodes[i].d);
      construct_bottom_up_finish();

      if (get_nodes().size() > 1)
      {
        //printf("create octree for node (%f %f %f) - %d nodes\n", large_nodes[i].p.x, large_nodes[i].p.y, large_nodes[i].p.z,
        //       (int)(get_nodes().size()));
        unsigned off = all_nodes.size();
        all_nodes[i].offset = off;
        all_nodes.reserve(all_nodes.size() + get_nodes().size() - 1);
        all_nodes.insert(all_nodes.end(), get_nodes().begin()+1, get_nodes().end()); //skip root node of get_nodes() as it is already in all_nodes
        for (int j=off; j<all_nodes.size(); j++) 
        {
          if (!is_leaf(all_nodes[j].offset))
            all_nodes[j].offset += off-1;
        }
      }
    }
  }

  get_nodes() = all_nodes;
  
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  
  float time_1 = 1e-3f*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  unsigned long n_real = get_nodes().size();
  unsigned long n_max = pow(8, settings.depth);
  
  printf("SDF octree created with %lu (%5.2lf%%) nodes (dense one would have %lu)\n", n_real, 100.0*n_real/n_max, n_max);
  printf("time spent (ms) %.1f\n", time_1);
}

void SparseOctreeBuilder::convert_to_frame_octree(std::vector<SdfFrameOctreeNode> &out_frame)
{
  auto &nodes = get_nodes();
  out_frame.resize(nodes.size());
  fill_octree_frame_rec(sdf, nodes, out_frame, 0, float3(0,0,0), 1);
}

void SparseOctreeBuilder::construct_bottom_up_frame(std::function<T(const float3 &)> _sdf, SparseOctreeSettings _settings, 
                                                    std::vector<SdfFrameOctreeNode> &out_frame)
{
  sdf = _sdf;
  settings = _settings;
  octree_f->get_nodes().clear();

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  construct_bottom_up_base(0u, float3(0,0,0), 1.0f);
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  construct_bottom_up_finish();
std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

  convert_to_frame_octree(out_frame);
  
std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
  
  float time_1 = 1e-3f*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  float time_2 = 1e-3f*std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
  float time_3 = 1e-3f*std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
  unsigned long n_real = get_nodes().size();
  unsigned long n_max = pow(8, settings.depth);
  
  printf("SDF octree created with %lu (%5.2lf%%) nodes (dense one would have %lu)\n", n_real, 100.0*n_real/n_max, n_max);
  printf("time spent (ms) %.1f %.1f %.1f\n", time_1, time_2, time_3);
}