#pragma once
#include "upg.h"
#include "generation_common.h"
#include <memory>

namespace upg
{
  class ParamNode
  {
  public:
    enum ParamNodeType
    {
      UNDEFINED,
      PRIMITIVE,
      ADD,
      SUB,
      MULT,
      DIV,
      INC,
      CMULT,
      CONST,
      NEG,
      PARAM_TYPES_COUNT
    };

    ParamNode(unsigned id, std::vector<float> p) { ID = id; param = p; }
    virtual ~ParamNode() = default;
    unsigned get_ID() const { return ID; }
    std::string get_node_name() const { return name; }

    virtual float get_res(std::vector<float> *jac) const = 0;
    bool add_source(ParamNode *node)// returns the availability of free space
    {
        if (source_cnt() > sources.size())
        {
            sources.push_back(node);
        }
        return source_cnt() > sources.size();
    }
    bool check_source()// returns the availability of free space
    {
        return source_cnt() > sources.size();
    }
    virtual unsigned param_cnt() const = 0;
    virtual unsigned source_cnt() const = 0;
    virtual std::vector<const ParamNode *> get_sources() { return sources; }
    void set_param_span(std::span<my_float> s) { in_p = s; }
  protected:
    unsigned ID;
    std::string name;
    std::vector<float> param;
    std::vector<const ParamNode *> sources;
    std::span <const float> in_p;
  };

  class ParamsGraph
  {
  public:
    ParamsGraph() {}
    ParamsGraph(const ParamsGraph& p): all_nodes(p.all_nodes)
    {
      roots = p.roots;
      structure = p.structure;
    }
    ParamsGraph(std::vector<ParamNode *> _roots, const std::vector<std::vector<float>> &_structure,
                std::span<const std::unique_ptr<ParamNode>> _all_nodes):
      all_nodes(_all_nodes)
    {
      for (int i = 0; i < _roots.size(); ++i)
      {
        roots.push_back(_roots[i]);
      }
      structure = _structure;
    }
    std::vector<float> get_params(std::vector<float> *jac = nullptr) const
    {
      std::vector<float> res_params;
      if (jac)
      {
        jac->clear();
      }
      for (int i = 0; i < roots.size(); ++i)
      {
        std::vector <float> *addr = NULL;
        std::vector <float> jac_part;
        if (jac)
        {
          addr = &jac_part;
        }
        res_params.push_back(roots[i]->get_res(addr));
        if (jac)
        {
          jac->insert(jac->end(), jac_part.begin(), jac_part.end());
        }
      }
      return res_params;
    }
    std::vector<ParamNode *> roots;
    std::vector<std::vector<float>> structure;
  private:
    std::span<const std::unique_ptr<ParamNode>> all_nodes;
  };

  class ParamGenInstance
  {
  public:
    ParamGenInstance(std::vector<std::vector<float>> structure, unsigned in, unsigned out);
    void recreate(std::vector<std::vector<float>> &structure, unsigned in, unsigned out);
    ParamsGraph get_graph(std::span<const float> parameters);

  private:
    std::vector<std::unique_ptr<ParamNode>> all_nodes;
    std::vector<ParamNode *> roots;
    
    std::vector<my_float> all_params;
    std::vector<std::vector<float>> structure;//{ {ID, type, params}, {ID, type, params}, {ID//if this id already in structure}, ...}
  };

  ParamNode *param_node_by_node_type_id(unsigned id, uint16_t num, std::vector <float> p);
}