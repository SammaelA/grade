#pragma once
#include "tree_node.h"

namespace upg
{
  class Tree
  {
    GenNode *root;
    std::vector<my_float> all_params;

    GenNode *node_by_number(uint16_t num, unsigned id);
    void tree_del();

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
  };
}