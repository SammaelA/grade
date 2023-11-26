#pragma once
#include "upg.h"
#include "generation_common.h"
namespace upg
{
  struct UniversalGenJacobian;
  struct UniversalGenMesh;
  class GenNode
  {
  protected:
    unsigned ID;
    std::string name;
    std::span<const float> p;
  public:
    GenNode(unsigned id) { ID = id; }
    virtual ~GenNode() = default;
    unsigned get_ID() const
    {
      return ID;
    }
    std::string get_node_name() const
    {
      return name;
    }
    void set_param_span(std::span<my_float> s)
    {
      p = s;
    }
    virtual UniversalGenMesh apply(UniversalGenJacobian *out_jac) = 0;
    virtual bool add_child(GenNode *node) = 0;// returns the availability of free space
    virtual unsigned param_cnt() const = 0;
    virtual unsigned child_cnt() const = 0;
    virtual std::vector<const GenNode *> get_children() const = 0;
    virtual std::vector<ParametersDescription::Param> get_parameters_block() const
    {
      logerr("AAAAA");
      return {};
    }
  };
  GenNode *node_by_node_type_id(uint16_t num, unsigned id);
}