#pragma once
#include "tree_node.h"

namespace upg
{
  class Tree
  {
  public:
    UniversalGenMesh generate()
    {
      return root->apply();
    }
    void create(const UPGStructure &structure);
    void take_params(std::span<const float> param);
    ~Tree()
    {
      tree_del();
    }
  private:
    GenNode *root;
    std::vector<my_float> all_params;
    void tree_del();
  };
}