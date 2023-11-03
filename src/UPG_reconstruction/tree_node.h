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
    std::string get_node_name()
    {
      return name;
    }
    void set_param_span(std::span<my_float> s)
    {
      p = s;
    }
    virtual UniversalGenMesh apply() = 0;
    virtual bool add_child(GenNode *node) = 0;// returns the availability of free space
    virtual unsigned param_cnt() = 0;
    virtual unsigned child_cnt() = 0;
    virtual std::vector<GenNode *> childs() = 0;
    virtual std::vector<ParametersDescription::Param> get_parameters_block()
    {
      logerr("AAAAA");
      return {};
    }
  };
  GenNode *node_by_node_type_id(uint16_t num, unsigned id);
}