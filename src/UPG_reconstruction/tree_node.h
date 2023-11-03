#pragma once
#include "reconstruction.h"
#include "models.h"
namespace upg
{
  class GenNode
  {
  protected:
    unsigned ID;
    std::string name;
    uint16_t node_num;
    std::span<const float> p;
  public:
    GenNode(unsigned id) { ID = id; }
    unsigned get_ID()
    {
      return ID;
    }
    uint16_t get_node_num()
    {
      return node_num;
    }
    void set_param_span(std::span<my_float> s)
    {
      p = s;
    }
    virtual UniversalGenMesh apply() = 0;
    virtual bool add_child(GenNode *node) = 0;// returns the availability of free space
    inline virtual unsigned param_cnt() = 0;
    inline virtual unsigned child_cnt() = 0;
    virtual std::vector<GenNode *> childs() = 0;
  };
  GenNode *node_by_node_type_id(uint16_t num, unsigned id);
}