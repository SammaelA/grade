#include "sparse_octree.h"
#include "distribution.h"
#include <cassert>

bool SparseOctree::is_border(float distance, unsigned level)
{
  return level < 2  ? true : abs(distance) < sqrt(2)*(1/(pow(2, level-2)));
}

void SparseOctree::add_node_rec(std::function<T(const float3 &)> f,
                                unsigned node_idx,
                                unsigned depth,
                                unsigned max_depth,
                                float3 p,
                                float d)
{
  nodes[node_idx].value = f(2.0f * ((p + float3(0.5, 0.5, 0.5)) * d) - float3(1, 1, 1));

  if (depth < max_depth)
  {
    nodes[node_idx].offset = nodes.size();
    
    nodes.resize(nodes.size() + 8);
    unsigned idx = nodes[node_idx].offset;
    for (unsigned cid = 0; cid < 8; cid++)
    {
      add_node_rec(f, idx + cid, depth + 1, max_depth, 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1), d / 2);
    }
  }
}

bool DBG = false;

void SparseOctree::split_children(std::function<T(const float3 &)> f,
                                  unsigned node_idx,
                                  float threshold,
                                  float3 p,
                                  float d,
                                  unsigned level)
{
  assert(nodes[node_idx].offset > 0);
  unsigned idx = nodes[node_idx].offset;
  T n_distances[8];
  for (unsigned cid = 0; cid < 8; cid++)
    n_distances[cid] = nodes[idx + cid].value;
  for (unsigned cid = 0; cid < 8; cid++)
  {
    if (nodes[idx + cid].offset > 0) // go deeper
    {
      split_children(f, idx + cid, threshold, 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1), d/2, level+1);
    }
    else if (is_border(abs(nodes[idx + cid].value), level)) // child is leaf, check if we should split it
    {
      //printf("LOL %u %u %u\n",(cid & 4) >> 2, (cid & 2) >> 1, cid & 1);
      float3 p1 = 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1);
      float3 p2 = (p1+float3(0.5,0.5,0.5)) * (d/2);
      float d1 = f(2.0f*p2-1.0f);
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
        float d_ref = f(2.0f*pos-1.0f);
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
          nodes[gc_idx + gcid].value = f(2.0f*pos-1.0f);
          //printf("%u putf %f %f %f -- %f\n",gc_idx + gcid, pos.x,pos.y,pos.z, nodes[gc_idx + gcid].value);
          //add_node_rec(f, idx + cid, depth + 1, max_depth, 2 * p + float3((cid & 4) >> 2, (cid & 2) >> 1, cid & 1), d / 2);
        }        
      }
    }
  }
  //    add_node_rec(f, octree, idx+cid, depth+1, max_depth, 2*p + float3((cid&4)<<2,(cid&2)<<1,cid&1), d/2);
}

void SparseOctree::construct(std::function<T(const float3 &)> f,
                             SparseOctreeSettings settings)
{
  //DBG = true;
  //DBG = false;
  nodes.clear();
  nodes.reserve(std::min(10000u, settings.max_nodes));

  // adding initial nodes
  nodes.emplace_back();
  add_node_rec(f, 0, 0, settings.min_depth, float3(0, 0, 0), 1.0f);

  // splitting octree threshold is reached
  //if (settings.threshold > 0)
  unsigned prev_size = 0;
  while (prev_size < nodes.size())
  {
    prev_size = nodes.size();
    split_children(f, 0, 0.005, float3(0,0,0), 1, 1);
  }
  
  print_stat();
  for (int i=1;i<=10;i++)
  {
    auto p = estimate_quality(f, 0.002f*i);
    printf("estimate_quality(%f): %f %f\n",0.002f*i, (float)p.first, (float)p.second);
  }
  //split_children(f, 0, 0.01, float3(0,0,0), 1, 1);
  //split_children(f, 0, 0.01, float3(0,0,0), 1, 1);
  //DBG = true;
  //split_children(f, 0, 0.01, float3(0,0,0), 1, 1);
  //DBG = false;
  
  //DBG = true;
}

SparseOctree::T SparseOctree::sample(const float3 &position, unsigned max_level) const
{
  if (DBG) printf("start pos %f %f %f\n", position.x, position.y, position.z);
  float3 n_pos = LiteMath::clamp(0.5f*(position + 1.0f), 0.0f, 1.0f);//position in current neighborhood
  float d = 1;//size of current neighborhood
  T n_distances[8];
  unsigned n_indices[8]; //0 index means that it is "virtual" node
  unsigned non_leaf_nodes = 0; //how many nodes in neighborhood have children

  T prev_n_distances[8];
  unsigned prev_n_indices[8]; //0 index means that it is "virtual" node

  //start with root's childer as a neighborhood
  unsigned r_idx = nodes[0].offset;
  for (int i=0;i<8;i++)
  {
    n_distances[i] = nodes[r_idx+i].value;
    n_indices[i] = r_idx+i;
    non_leaf_nodes += (nodes[r_idx+i].offset > 0);
  }

  int level = 1;
  while (non_leaf_nodes > 0 && level <= max_level)
  {
    if (DBG) printf("level %d\n", level);
    for (int i=0;i<8;i++)
    {
      prev_n_distances[i] = n_distances[i];
      prev_n_indices[i] = n_indices[i];
    }
    //go 1 level deeper every iteration
    non_leaf_nodes = 0;

    uint3 idx8 = uint3(8.0f*fract(n_pos));
    //if (idx8.x + idx8.y + idx8.z > 0)
    //  printf("%u %u %u\n", idx8.x, idx8.y, idx8.z);
    //assert(LiteMath::all_of(idx8 >= uint3(2u)) && LiteMath::all_of(idx8 < uint3(6u)));
    if (DBG) printf("idx8 %u %u %u\n", idx8.x, idx8.y, idx8.z);
    uint3 pidx[2], chidx[2];
    float3 n_pos_sh;
    for (int i=0;i<3;i++)
    {
      if (idx8[i] == 0)
        { pidx[0][i] = 0; chidx[0][i] = 0; pidx[1][i] = 0; chidx[1][i] = 0; n_pos_sh[i] = 0;}
      else if (idx8[i] <= 2)
        { pidx[0][i] = 0; chidx[0][i] = 0; pidx[1][i] = 0; chidx[1][i] = 1; n_pos_sh[i] = 0;}
      else if (idx8[i] <= 4)
        { pidx[0][i] = 0; chidx[0][i] = 1; pidx[1][i] = 1; chidx[1][i] = 0; n_pos_sh[i] = 0.25;}
      else if (idx8[i] <= 6)
        { pidx[0][i] = 1; chidx[0][i] = 0; pidx[1][i] = 1; chidx[1][i] = 1; n_pos_sh[i] = 0.5;}
      else //if (idx8[i] == 7)
        { pidx[0][i] = 1; chidx[0][i] = 1; pidx[1][i] = 1; chidx[1][i] = 1; n_pos_sh[i] = 0.5;}
    }

    if (DBG) printf("pidx %u %u %u -- %u %u %u\n", pidx[0].x, pidx[0].y, pidx[0].z, pidx[1].x, pidx[1].y, pidx[1].z);
    if (DBG) printf("chidx %u %u %u -- %u %u %u\n", chidx[0].x, chidx[0].y, chidx[0].z, chidx[1].x, chidx[1].y, chidx[1].z);
    if (DBG) printf("npos %f %f %f\n",n_pos.x, n_pos.y, n_pos.z);
    
    //create new neighborhood
    for (unsigned i=0;i<8;i++)
    {
      uint3 n_idx((i & 4) >> 2, (i & 2) >> 1, i & 1);
      uint3 cur_pidx = uint3(pidx[n_idx[0]][0], pidx[n_idx[1]][1], pidx[n_idx[2]][2]);
      uint3 cur_chidx = uint3(chidx[n_idx[0]][0], chidx[n_idx[1]][1], chidx[n_idx[2]][2]);

      unsigned p_index = 4*cur_pidx.x + 2*cur_pidx.y + cur_pidx.z;
      if (prev_n_indices[p_index] > 0 &&                //p_index is a real node
          nodes[prev_n_indices[p_index]].offset > 0)    //p_index has children
      {
        if (DBG) printf("child %u real\n", i);
        unsigned ch_index = nodes[prev_n_indices[p_index]].offset + 4*cur_chidx.x + 2*cur_chidx.y + cur_chidx.z;
        n_distances[i] = nodes[ch_index].value;
        n_indices[i] = ch_index;      
        non_leaf_nodes += nodes[ch_index].offset > 0;
      }
      else                                              //p_index is a leaf node
      {
        //assert(prev_n_indices[p_index] > 0);
        if (DBG) printf("child %u fake\n", i);
        float3 smp_pos = 0.5f*float3(cur_pidx) + 0.25f*float3(cur_chidx);
        n_distances[i] = prev_n_distances[p_index];
        n_indices[i] = 0;   
        //assert(false);
      }
    }

    n_pos = fract(2.0f*(n_pos - n_pos_sh));
    d /= 2;

    level++;
  }

  //float d1 = sample_closest(position);
  //unsigned ch_index = 4*(n_pos.x >= 0.5) + 2*(n_pos.y >= 0.5) + (n_pos.z >= 0.5);
  //float d2 = n_distances[ch_index];
  //if (std::abs(d1 - d2) > 1e-5)
  //  printf("AAAA %f %f %f\n", position.x, position.y, position.z);
  if (DBG)
  {
    printf("%f %f  %f %f\n%f %f  %f %f\n", n_distances[2], n_distances[6], n_distances[3], n_distances[7],
    n_distances[0],n_distances[4],n_distances[1],n_distances[5]);
    printf("npos %f %f %f\n", n_pos.x, n_pos.y, n_pos.z);
  }
  //if (DBG) printf("D = %f\n",sample_neighborhood_bilinear(n_pos, n_distances));
  return sample_neighborhood_bilinear(clamp(2.0f*n_pos - float3(0.5, 0.5, 0.5), 0.0f, 1.0f), n_distances);
}


SparseOctree::T SparseOctree::sample_2(const float3 &position, unsigned max_level) const
{
  if (DBG) printf("start pos %f %f %f\n", position.x, position.y, position.z);
  float3 n_pos = LiteMath::clamp(0.5f*(position + 1.0f), 0.0f, 1.0f);//position in current neighborhood
  float d = 1;//size of current neighborhood
  T n_distances[8];
  unsigned n_indices[8]; //0 index means that it is "virtual" node
  unsigned non_leaf_nodes = 0; //how many nodes in neighborhood have children

  T prev_n_distances[8];
  unsigned prev_n_indices[8]; //0 index means that it is "virtual" node

  //start with root's childer as a neighborhood
  unsigned r_idx = nodes[0].offset;
  for (int i=0;i<8;i++)
  {
    n_distances[i] = nodes[r_idx+i].value;
    n_indices[i] = r_idx+i;
    non_leaf_nodes += (nodes[r_idx+i].offset > 0);
  }

  int level = 1;
  while (non_leaf_nodes > 0 && level <= max_level)
  {
    if (DBG) printf("level %d\n", level);
    for (int i=0;i<8;i++)
    {
      prev_n_distances[i] = n_distances[i];
      prev_n_indices[i] = n_indices[i];
    }
    //go 1 level deeper every iteration
    non_leaf_nodes = 0;

    uint3 idx8 = uint3(8.0f*fract(n_pos));
    //if (idx8.x + idx8.y + idx8.z > 0)
    //  printf("%u %u %u\n", idx8.x, idx8.y, idx8.z);
    //assert(LiteMath::all_of(idx8 >= uint3(2u)) && LiteMath::all_of(idx8 < uint3(6u)));
    if (DBG) printf("idx8 %u %u %u\n", idx8.x, idx8.y, idx8.z);
    uint3 pidx[2], chidx[2];
    float3 n_pos_sh;
    for (int i=0;i<3;i++)
    {
      if (idx8[i] == 0)
        { pidx[0][i] = 0; chidx[0][i] = 0; pidx[1][i] = 0; chidx[1][i] = 0; n_pos_sh[i] = 0;}
      else if (idx8[i] <= 2)
        { pidx[0][i] = 0; chidx[0][i] = 0; pidx[1][i] = 0; chidx[1][i] = 1; n_pos_sh[i] = 0;}
      else if (idx8[i] <= 4)
        { pidx[0][i] = 0; chidx[0][i] = 1; pidx[1][i] = 1; chidx[1][i] = 0; n_pos_sh[i] = 0.25;}
      else if (idx8[i] <= 6)
        { pidx[0][i] = 1; chidx[0][i] = 0; pidx[1][i] = 1; chidx[1][i] = 1; n_pos_sh[i] = 0.5;}
      else //if (idx8[i] == 7)
        { pidx[0][i] = 1; chidx[0][i] = 1; pidx[1][i] = 1; chidx[1][i] = 1; n_pos_sh[i] = 0.5;}
    }

    if (DBG) printf("pidx %u %u %u -- %u %u %u\n", pidx[0].x, pidx[0].y, pidx[0].z, pidx[1].x, pidx[1].y, pidx[1].z);
    if (DBG) printf("chidx %u %u %u -- %u %u %u\n", chidx[0].x, chidx[0].y, chidx[0].z, chidx[1].x, chidx[1].y, chidx[1].z);
    if (DBG) printf("npos %f %f %f\n",n_pos.x, n_pos.y, n_pos.z);
    
    //create new neighborhood
    for (unsigned i=0;i<8;i++)
    {
      uint3 n_idx((i & 4) >> 2, (i & 2) >> 1, i & 1);
      uint3 cur_pidx = uint3(pidx[n_idx[0]][0], pidx[n_idx[1]][1], pidx[n_idx[2]][2]);
      uint3 cur_chidx = uint3(chidx[n_idx[0]][0], chidx[n_idx[1]][1], chidx[n_idx[2]][2]);

      unsigned p_index = 4*cur_pidx.x + 2*cur_pidx.y + cur_pidx.z;
      if (prev_n_indices[p_index] > 0 &&                //p_index is a real node
          nodes[prev_n_indices[p_index]].offset > 0)    //p_index has children
      {
        if (DBG) printf("child %u real\n", i);
        unsigned ch_index = nodes[prev_n_indices[p_index]].offset + 4*cur_chidx.x + 2*cur_chidx.y + cur_chidx.z;
        n_distances[i] = nodes[ch_index].value;
        n_indices[i] = ch_index;      
        non_leaf_nodes += nodes[ch_index].offset > 0;
      }
      else                                              //p_index is a leaf node
      {
        assert(prev_n_indices[p_index] > 0);
        if (DBG) printf("child %u fake\n", i);
        float3 smp_pos = 0.5f*float3(cur_pidx) + 0.25f*float3(cur_chidx);
        n_distances[i] = sample_neighborhood_bilinear(smp_pos, prev_n_distances);
        n_indices[i] = 0;   
        //assert(false);
      }
    }

    n_pos = fract(2.0f*(n_pos - n_pos_sh));
    d /= 2;

    level++;
  }

  //float d1 = sample_closest(position);
  //unsigned ch_index = 4*(n_pos.x >= 0.5) + 2*(n_pos.y >= 0.5) + (n_pos.z >= 0.5);
  //float d2 = n_distances[ch_index];
  //if (std::abs(d1 - d2) > 1e-5)
  //  printf("AAAA %f %f %f\n", position.x, position.y, position.z);
  if (DBG)
  {
    printf("%f %f  %f %f\n%f %f  %f %f\n", n_distances[2], n_distances[6], n_distances[3], n_distances[7],
    n_distances[0],n_distances[4],n_distances[1],n_distances[5]);
    printf("npos %f %f %f\n", n_pos.x, n_pos.y, n_pos.z);
  }
  //if (DBG) printf("D = %f\n",sample_neighborhood_bilinear(n_pos, n_distances));
  return sample_neighborhood_bilinear(clamp(2.0f*n_pos - float3(0.5, 0.5, 0.5), 0.0f, 1.0f), n_distances);
}

SparseOctree::T SparseOctree::sample_neighborhood_bilinear(const float3 &dp, T n_distances[8]) const
{
  //unsigned ch_index = 4*(dp.x >= 0.5) + 2*(dp.y >= 0.5) + (dp.z >= 0.5);
  //return n_distances[ch_index];
  //float t = 0.620838*0.620838 + 0.323323*0.323323 + 0.627405*0.627405;
    if (DBG) printf("sample %f %f %f\n", dp.x, dp.y, dp.z);
  return (1-dp.x)*(1-dp.y)*(1-dp.z)*n_distances[0] + 
         (1-dp.x)*(1-dp.y)*(  dp.z)*n_distances[1] + 
         (1-dp.x)*(  dp.y)*(1-dp.z)*n_distances[2] + 
         (1-dp.x)*(  dp.y)*(  dp.z)*n_distances[3] + 
         (  dp.x)*(1-dp.y)*(1-dp.z)*n_distances[4] + 
         (  dp.x)*(1-dp.y)*(  dp.z)*n_distances[5] + 
         (  dp.x)*(  dp.y)*(1-dp.z)*n_distances[6] + 
         (  dp.x)*(  dp.y)*(  dp.z)*n_distances[7];
}

SparseOctree::T SparseOctree::sample_closest(const float3 &position) const
{
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

void check_border_reappearance_rec(const std::vector<SparseOctree::Node> &nodes, unsigned idx, unsigned level)
{
  bool b = SparseOctree::is_border(nodes[idx].value, level); 
  if (nodes[idx].offset > 0)
  {
    for (unsigned ch_idx=0; ch_idx<8; ch_idx++)
    {
      bool ch_b = SparseOctree::is_border(nodes[nodes[idx].offset + ch_idx].value, level+1);
      if (!b && ch_b)
        printf("reappeared border (level %d-%d): %f %f\n", level,level+1, nodes[idx].value, nodes[nodes[idx].offset + ch_idx].value);
      check_border_reappearance_rec(nodes, nodes[idx].offset + ch_idx, level+1);
    }
  }
}

void print_stat_rec(const std::vector<SparseOctree::Node> &nodes, SparseOctreeCounts &counts, unsigned idx, unsigned level)
{
  if (level == counts.count_all.size())
  {
    counts.count_all.push_back(0);
    counts.count_border.push_back(0);
    counts.count_leaf.push_back(0);
    counts.count_border_leaf.push_back(0);
  }
  counts.count_all[level] += 1;
  counts.count_border[level] += SparseOctree::is_border(nodes[idx].value, level);
  counts.count_leaf[level] += nodes[idx].offset==0;
  counts.count_border_leaf[level] += SparseOctree::is_border(nodes[idx].value, level) && nodes[idx].offset==0;
  if (nodes[idx].offset > 0)
  {
    for (unsigned ch_idx=0; ch_idx<8; ch_idx++)
      print_stat_rec(nodes, counts, nodes[idx].offset + ch_idx, level+1);
  }
}

void SparseOctree::print_stat() const
{
  check_border_reappearance_rec(nodes, 0, 0);

  SparseOctreeCounts counts;
  print_stat_rec(nodes, counts, 0, 0);
  printf("SparseOctree::print_stat\n");
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

std::pair<float,float> SparseOctree::estimate_quality(std::function<T(const float3 &)> reference_f, float dist_thr, unsigned samples) const
{
  unsigned i = 0;
  std::vector<float> differences(samples, 0);
  while (i < samples)
  {
    float3 p = float3(urand(-1,1),urand(-1,1),urand(-1,1));
    float ref_d = reference_f(p);
    if (abs(ref_d) < dist_thr)
    {
      differences[i] = sample(p) - ref_d;
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
  return {av_diff, max_diff};
}