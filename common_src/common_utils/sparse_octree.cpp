#include "sparse_octree.h"
#include "distribution.h"
#include "stat_box.h"
#include <cassert>
#include <map>

bool SparseOctree::is_border(float distance, unsigned level)
{
  return level < 2  ? true : abs(distance) < sqrt(2)*(1/(pow(2, level-2)));
}

constexpr unsigned INVALID_IDX = 1u<<31u;
bool is_leaf(unsigned offset)
{
  return (offset == 0) || (offset & INVALID_IDX);
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
  assert(!is_leaf(nodes[node_idx].offset));
  unsigned idx = nodes[node_idx].offset;
  T n_distances[8];
  for (unsigned cid = 0; cid < 8; cid++)
    n_distances[cid] = nodes[idx + cid].value;
  for (unsigned cid = 0; cid < 8; cid++)
  {
    if (!is_leaf(nodes[idx + cid].offset)) // go deeper
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

void SparseOctree::construct_top_down(std::function<T(const float3 &)> f,
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
  return sample_mip_skip_3x3(position, max_level);
}
SparseOctree::T SparseOctree::sample_mip_skip_closest(const float3 &position, unsigned max_level) const
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
    non_leaf_nodes += (!is_leaf(nodes[r_idx+i].offset));
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
          !is_leaf(nodes[prev_n_indices[p_index]].offset))    //p_index has children
      {
        if (DBG) printf("child %u real\n", i);
        unsigned ch_index = nodes[prev_n_indices[p_index]].offset + 4*cur_chidx.x + 2*cur_chidx.y + cur_chidx.z;
        n_distances[i] = nodes[ch_index].value;
        n_indices[i] = ch_index;      
        non_leaf_nodes += !is_leaf(nodes[ch_index].offset);
      }
      else                                              //p_index is a leaf node
      {
        if (DBG) printf("child %u fake\n", i);
        n_distances[i] = prev_n_distances[p_index];
        n_indices[i] = 0;   
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

static inline unsigned idx3x3(int3 c)
{
  return 9*c.x + 3*c.y + c.z;
}

enum Overshoot
{
  X_L = 1<<0,
  X_H = 1<<1,
  Y_L = 1<<2,
  Y_H = 1<<3,
  Z_L = 1<<4,
  Z_H = 1<<5
};
struct Neighbor
{
  SparseOctree::Node node;
  unsigned char overshoot;
};
static std::vector<unsigned char> overshoot_array = {
  X_L|Y_L|Z_L, X_L|Y_L| 0 , X_L|Y_L|Z_H,
  X_L| 0 |Z_L, X_L| 0 | 0 , X_L| 0 |Z_H,
  X_L|Y_H|Z_L, X_L|Y_H| 0 , X_L|Y_H|Z_H,

   0 |Y_L|Z_L,  0 |Y_L| 0 ,  0 |Y_L|Z_H,
   0 | 0 |Z_L,  0 | 0 | 0 ,  0 | 0 |Z_H,
   0 |Y_H|Z_L,  0 |Y_H| 0 ,  0 |Y_H|Z_H,

  X_H|Y_L|Z_L, X_H|Y_L| 0 , X_H|Y_L|Z_H,
  X_H| 0 |Z_L, X_H| 0 | 0 , X_H| 0 |Z_H,
  X_H|Y_H|Z_L, X_H|Y_H| 0 , X_H|Y_H|Z_H,
};

float sample_neighborhood(Neighbor *neighbors, float3 n_pos)
{
  float3 qx = clamp(float3(0.5-n_pos.x,std::min(0.5f + n_pos.x, 1.5f - n_pos.x),-0.5+n_pos.x),0.0f,1.0f);
  float3 qy = clamp(float3(0.5-n_pos.y,std::min(0.5f + n_pos.y, 1.5f - n_pos.y),-0.5+n_pos.y),0.0f,1.0f);
  float3 qz = clamp(float3(0.5-n_pos.z,std::min(0.5f + n_pos.z, 1.5f - n_pos.z),-0.5+n_pos.z),0.0f,1.0f);

  if (DBG) 
  {
    printf("sample %f %f %f\n",n_pos.x, n_pos.y, n_pos.z);
    printf("qx %f %f %f\nqy %f %f %f\nqz %f %f %f\n", qx.x,qx.y,qx.z, qy.x,qy.y,qy.z, qz.x,qz.y,qz.z);

    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3;j++)
      {
        for (int k=0;k<3;k++)
        {
          float v = neighbors[9*i + 3*j + k].node.value;
          printf("%f ", v);
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  float res = 0.0;
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
        res += qx[i]*qy[j]*qz[k]*neighbors[9*i + 3*j + k].node.value;
  return res;
}

SparseOctree::T SparseOctree::sample_mip_skip_3x3(const float3 &position, unsigned max_level) const
{
  constexpr unsigned CENTER = 9 + 3 + 1;
  constexpr float EPS = 1e-6;
  Neighbor neighbors[27];
  Neighbor new_neighbors[27];


  float3 n_pos = LiteMath::clamp(0.5f*(position + 1.0f), EPS, 1.0f-EPS);//position in current neighborhood
  float d = 1;//size of current neighborhood
  unsigned level = 0;
  for (int i=0;i<27;i++)
  {
    neighbors[i].node = nodes[0];
    neighbors[i].overshoot = overshoot_array[i];
  }
  neighbors[CENTER].overshoot = 0;

  while (!is_leaf(neighbors[CENTER].node.offset) && level < max_level)
  {
    int3 ch_shift = int3(n_pos.x >= 0.5, n_pos.y >= 0.5, n_pos.z >= 0.5);

    for (int i=0;i<27;i++)
    {
      int3 n_offset = int3(i/9, i/3%3, i%3); //[0,2]^3
      int3 p_idx = (n_offset + ch_shift + 1) / 2;
      int3 ch_idx = (n_offset + ch_shift + 1) - 2*p_idx;
      unsigned p_offset = idx3x3(p_idx);
      if (DBG) printf("level %u, npos (%f %f %f), ch_sift (%d %d %d), i %d p_idx (%d %d %d) ch_idx (%d %d %d)\n",
            level, n_pos.x, n_pos.y, n_pos.z,  ch_shift.x, ch_shift.y, ch_shift.z,
             i,  p_idx.x, p_idx.y, p_idx.z,  ch_idx.x, ch_idx.y, ch_idx.z);
    
      if (is_leaf(neighbors[p_offset].node.offset)) //resample
      {
        float3 rs_pos = 0.5f*float3(2*p_idx + ch_idx) - 1.0f + 0.25f;//in [-1,2]^3
        if (DBG) printf("resample, rs_pos %f %f %f\n", rs_pos.x, rs_pos.y, rs_pos.z);
        new_neighbors[i].node.value = sample_neighborhood(neighbors, rs_pos);
        new_neighbors[i].node.offset = 0;
        new_neighbors[i].overshoot = 0;
      }
      else if (is_leaf(neighbors[p_offset].overshoot)) //pick child node
      {
        if (DBG) printf("base\n");
        unsigned ch_offset = 4*ch_idx.x + 2*ch_idx.y + ch_idx.z;
        unsigned off = neighbors[p_offset].node.offset;
        new_neighbors[i].node = nodes[off + ch_offset];
        new_neighbors[i].overshoot = 0;
        if (DBG) printf("child[%d %d %d] = parent(%u)[%d %d %d][%d %d %d]\n", i/9,i/3%3,i%3, neighbors[p_offset].node.offset, 
                        p_idx.x, p_idx.y, p_idx.z,
                        ch_idx.x, ch_idx.y, ch_idx.z);
      }
      else //pick child node, but mind the overshoot
      {
        /**/
        if (DBG) printf("overshoot\n");
        assert(p_offset != CENTER);
        int3 ch_idx_overshoot = ch_idx;
        unsigned char osh = neighbors[p_offset].overshoot;
        unsigned char new_osh = 0;
        if ((osh&X_L) && p_idx.x == 0) 
          {ch_idx_overshoot.x = 0; new_osh |= X_L; }
        else if ((osh&X_H) && p_idx.x == 2) 
          {ch_idx_overshoot.x = 1; new_osh |= X_H; }
        if ((osh&Y_L) && p_idx.y == 0) 
          {ch_idx_overshoot.y = 0; new_osh |= Y_L; }
        else if ((osh&Y_H) && p_idx.y == 2) 
          {ch_idx_overshoot.y = 1; new_osh |= Y_H; }
        if ((osh&Z_L) && p_idx.z == 0) 
          {ch_idx_overshoot.z = 0; new_osh |= Z_L; }
        else if ((osh&Z_H) && p_idx.z == 2) 
          {ch_idx_overshoot.z = 1; new_osh |= Z_H; }

        unsigned ch_offset = 4*ch_idx_overshoot.x + 2*ch_idx_overshoot.y + ch_idx_overshoot.z;
        unsigned off = neighbors[p_offset].node.offset;
        new_neighbors[i].node = nodes[off + ch_offset];
        new_neighbors[i].overshoot = new_osh;
        if (DBG) printf("child[%d %d %d] = parent(%u)[%d %d %d][%d %d %d]\n", i/9,i/3%3,i%3, neighbors[p_offset].node.offset, 
                        p_idx.x, p_idx.y, p_idx.z,
                        ch_idx_overshoot.x, ch_idx_overshoot.y, ch_idx_overshoot.z);
      }
    }

    for (int i=0;i<27;i++)
      neighbors[i] = new_neighbors[i];

    n_pos = fract(2.0f*(n_pos - 0.5f*float3(ch_shift)));
    d /= 2;
    level++;
  }

  return sample_neighborhood(neighbors, n_pos);
}

SparseOctree::T SparseOctree::sample_mip_skip_2x2(const float3 &position, unsigned max_level) const
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
    non_leaf_nodes += (!is_leaf(nodes[r_idx+i].offset));
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
      if (prev_n_indices[p_index] > 0 &&                    //p_index is a real node
          !is_leaf(nodes[prev_n_indices[p_index]].offset))  //p_index has children
      {
        if (DBG) printf("child %u real\n", i);
        unsigned ch_index = nodes[prev_n_indices[p_index]].offset + 4*cur_chidx.x + 2*cur_chidx.y + cur_chidx.z;
        n_distances[i] = nodes[ch_index].value;
        n_indices[i] = ch_index;      
        non_leaf_nodes += !is_leaf(nodes[ch_index].offset);
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
  if (!is_leaf(nodes[idx].offset))
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
  counts.count_leaf[level] += is_leaf(nodes[idx].offset);
  counts.count_border_leaf[level] += SparseOctree::is_border(nodes[idx].value, level) && is_leaf(nodes[idx].offset);
  if (!is_leaf(nodes[idx].offset))
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
  stat::Bins<float> stat_box(0.0f, 0.1f);
  unsigned i = 0;
  std::vector<float> differences(samples, 0);
  while (i < samples)
  {
    float3 p = float3(urand(-1,1),urand(-1,1),urand(-1,1));
    float ref_d = reference_f(p);
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

float2 invalidate_node_rec(std::vector<SparseOctree::Node> &nodes, unsigned idx)
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

void remove_non_border_rec(std::vector<SparseOctree::Node> &nodes, unsigned min_level_to_remove, unsigned idx, unsigned level)
{
  if (is_leaf(nodes[idx].offset)) 
    return;
  
  bool is_border = false;
  for (int i=0;i<8;i++)
  {
    float dist = nodes[nodes[idx].offset + i].value;
    if (abs(dist) < sqrt(2)*(1/(pow(2, level-1))))
    {
      is_border = true;
      break;
    }
  }
  
  if (is_border || level < min_level_to_remove)
  {
    for (int i=0;i<8;i++)
      remove_non_border_rec(nodes, min_level_to_remove, nodes[idx].offset + i, level+1);
  }
  else
  {
    nodes[idx].value = invalidate_node_rec(nodes, idx).y;
    nodes[idx].offset |= INVALID_IDX;
  }
}

void remove_linear_rec(SparseOctree &octree, float thr, unsigned min_level_to_remove, unsigned idx, unsigned level, float3 p, float d)
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

void check_validity_rec(std::function<SparseOctree::T(const float3 &)> f, const SparseOctree &octree, unsigned idx, float3 p, float d)
{
  unsigned offset = octree.get_nodes()[idx].offset;
  if (is_leaf(offset))
  {
    float3 pos = 2.0f*((p + float3(0.5,0.5,0.5))*d) - 1.0f;
    float d1 = f(pos);
    float d2 = octree.sample_mip_skip_closest(pos);
    float d3 = octree.sample(pos);
    if (abs(d1-d2) > 1e-6 || abs(d1-d3) > 1e-6)
      printf("invalid sample in %f %f %f, d = %f %f %f\n",pos.x, pos.y, pos.z, d1,d2,d3);
  }
  else
  {
    for (int i=0;i<8;i++)
      check_validity_rec(f, octree, offset+i, 2*p + float3((i & 4) >> 2, (i & 2) >> 1, i & 1), d/2);
  }
}

void SparseOctree::construct_bottom_up(std::function<T(const float3 &)> f, SparseOctreeSettings settings)
{
  nodes.clear();
  nodes.reserve((8.0/7)*std::pow(8, settings.min_depth));

  // create regular grid
  nodes.emplace_back();
  add_node_rec(f, 0, 0, settings.min_depth, float3(0, 0, 0), 1.0f);

  //remove non-borders nodes
  remove_non_border_rec(nodes, 4, 0, 0);

  //remove when a little is changed
  remove_linear_rec(*this, 0.0001, 3, 0, 0, float3(0, 0, 0), 1.0f);

  //check_validity_rec(f, *this, 0, float3(0,0,0), 1);
  auto p = estimate_quality(f, 10, 100000);
  printf("estimate_quality: %f %f\n", (float)p.first, (float)p.second);

  for (auto &n : nodes)
    if (is_leaf(n.offset))
      n.offset = 0;
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

void find_active_blocks_rec(const SparseOctree &octree, std::vector<std::vector<bool>> &block_active, 
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

void restore_invalid_nodes_in_active_blocks_rec(SparseOctree &octree, std::vector<std::vector<bool>> &block_active, 
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

void node_values_to_blocks_rec(SparseOctree &octree, std::vector<std::vector<bool>> &block_active, 
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

void SparseOctree::construct_bottom_up_blocks(std::function<T(const float3 &)> f, BlockedSparseOctreeSettings settings, 
                                              BlockSparseOctree<T> &out_bso)
{
  unsigned size = pow(2, settings.max_depth_blocks)*std::max(BlockSparseOctree<float>::BLOCK_SIZE_X, std::max(BlockSparseOctree<float>::BLOCK_SIZE_Y, BlockSparseOctree<float>::BLOCK_SIZE_Z));
  unsigned level = log2(size);
  assert(size == (unsigned)(pow(2, level)));

  nodes.clear();
  nodes.reserve((8.0/7)*std::pow(8, level)); 

  // create regular grid
  nodes.emplace_back();
  add_node_rec(f, 0, 0, level, float3(0, 0, 0), 1.0f);

  //remove non-borders nodes
  remove_non_border_rec(nodes, settings.min_remove_level, 0, 0);

  //remove if difference is lower that threshold
  remove_linear_rec(*this, settings.remove_thr, settings.min_remove_level, 0, 0, float3(0, 0, 0), 1.0f);

  std::vector<std::vector<bool>> block_active;
  std::vector<std::vector<BlockSparseOctree<float>::BlockInfo>> all_blocks;
  std::vector<std::map<uint3, unsigned, cmpUint3>> block_idx_to_block_n;
  all_blocks.resize(settings.max_depth_blocks + 1);
  block_active.resize(settings.max_depth_blocks + 1);
  block_idx_to_block_n.resize(settings.max_depth_blocks + 1);

  //convert SparseOctree to BlockSparseOctree
  out_bso.top_mip = settings.max_depth_blocks;

  //create two root blocks and fill all_blocks with all possible blocks
  unsigned z_off = 0;
  while (z_off*BlockSparseOctree<float>::BLOCK_SIZE_Z != BlockSparseOctree<float>::BLOCK_SIZE_X)
  {
    all_blocks[settings.max_depth_blocks].push_back({uint3(0,0,z_off), settings.max_depth_blocks, 0}); 
    block_active[settings.max_depth_blocks].push_back(false);
    block_idx_to_block_n[settings.max_depth_blocks][uint3(0,0,z_off)] = z_off;
    add_all_blocks_rec(all_blocks, block_active, block_idx_to_block_n, settings.max_depth_blocks, z_off);
    z_off++;
  }

  //check all nodes by block to find active blocks
  find_active_blocks_rec(*this, block_active, block_idx_to_block_n, log2(BlockSparseOctree<float>::BLOCK_SIZE_X), settings.max_depth_blocks, 0, 0, uint3(0,0,0));

  //restore invalid nodes in active blocks
  restore_invalid_nodes_in_active_blocks_rec(*this, block_active, block_idx_to_block_n, 
                                             log2(BlockSparseOctree<float>::BLOCK_SIZE_X), settings.max_depth_blocks, 
                                             0, 0, uint3(0,0,0));

  //calculate offsets for active blocks and allocate memory for it
  unsigned offset = 0;
  for (int i=0;i<settings.max_depth_blocks + 1; i++)
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
                            log2(BlockSparseOctree<float>::BLOCK_SIZE_X), settings.max_depth_blocks, 
                            0, 0, uint3(0,0,0));

  //collect all active 
  for (int i=0;i<settings.max_depth_blocks + 1; i++)
  {
    for (int j=0;j<all_blocks[i].size();j++)
    {
      if (block_active[i][j])
        out_bso.blocks.push_back(all_blocks[i][j]);
    }
  }
  
  //DEBUG:  check is all the data for blocks is collected
  /*
  for (int i=0;i<settings.max_depth_blocks + 1; i++)
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

  //check_validity_rec(f, *this, 0, float3(0,0,0), 1);
  auto p = estimate_quality(f, 10, 100000);
  printf("estimate_quality: %f %f\n", (float)p.first, (float)p.second);

  for (auto &n : nodes)
    if (is_leaf(n.offset))
      n.offset = 0;

}