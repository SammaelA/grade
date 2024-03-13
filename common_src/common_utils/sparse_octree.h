#pragma once
#include <vector>
#include <functional>
#include "LiteMath_ext.h"
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
  unsigned max_depth = 100;
  unsigned min_depth = 1;
  float threshold = 0;
};

class SparseOctree
{
public:
  using T = float;
  using index_t = unsigned;
  struct Node
  {
    T value;
    index_t offset = 0; //offset for children (they are stored together). 0 offset means it's a leaf
    index_t parent = 0; //offset for parent, root node has itself as a parent
    float3 pos;
    float d;
  };

  void construct(std::function<T(const float3 &)> f, SparseOctreeSettings settings);
  T sample(const float3 &pos, unsigned max_level = 1000) const;
  T sample_closest(const float3 &pos) const;

protected:
  void add_node_rec(std::function<T(const float3 &)> f, unsigned node_idx, unsigned depth,
                    unsigned max_depth, float3 p, float d);
  void split_children(std::function<T(const float3 &)> f, unsigned node_idx,
                      float threshold, float3 p, float d, unsigned level);

  T sample_neighborhood_bilinear(const float3 &n_pos, T n_distances[8]) const;

  std::vector<Node> nodes; //0 node is root
};