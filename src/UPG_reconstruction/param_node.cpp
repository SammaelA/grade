#include "param_node.h"
#include <stdio.h>
#include "generation.h"
#include <set>
#include "autodiff/autodiff.h"

namespace upg
{
  ParamsGraph::ParamsGraph(std::vector<std::vector<float>> structure, unsigned in, unsigned out)
  {
    recreate(structure, in, out);
  }
  void ParamsGraph::recreate(std::vector<std::vector<float>> &structure, unsigned in, unsigned out)
  {
    this->structure = structure;
    all_params.clear();
    all_nodes.clear();
    int i, idx = 0;
    std::vector <int> stack;
    std::map <int, ParamNode *> nodes_by_id;
    all_params.resize(in);
    ParamNode *node;
    for (i = 0; i < structure.size(); ++i)
    {
      if (nodes_by_id.find(structure[i][0]) != nodes_by_id.end())
      {
        node = nodes_by_id[structure[i][0]];
      }
      else
      {
        std::vector <float> params(0);
        for (int j = 2; j < structure[i].size(); ++j)
        {
          params.push_back(structure[i][j]);
        }
        node = param_node_by_node_type_id(structure[i][0], structure[i][1], params);
        node->set_param_span(std::span<float>(all_params.data(), in));
        nodes_by_id[structure[i][0]] = node;
        stack.push_back(structure[i][0]);
        all_nodes.push_back(std::unique_ptr<ParamNode>(node));
      }
      if (i < out)
      {
        roots.push_back(node);
      }
      else
      {
        while (idx < stack.size() && !nodes_by_id[stack[idx]]->check_source())
        {
          ++idx;
        }
        if (idx < stack.size() && structure[i][0] != stack[idx])
        {
          nodes_by_id[stack[idx]]->add_source(node);
        }
      }
    }
    while (idx < stack.size() && !nodes_by_id[stack[idx]]->check_source())
    {
      ++idx;
    }
    if (idx < stack.size())
    {
      logerr("ERROR in param structure\n");
      return;
    }
  }

  void ParamsGraph::get_graph(std::span<const float> parameters)
  {
    for (int i = 0; i < all_params.size(); ++i)
    {
      if (i < parameters.size())
      {
        all_params[i] = parameters[i];
      }
      else
      {
        all_params[i] = 0;
      }
    }
  }

  class PrimitiveParamNode : public ParamNode
  {
    static constexpr int PARAM_IDX = 0;
  public:
    PrimitiveParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "input param"; }
    virtual unsigned source_cnt() const override { return 0; }
    virtual unsigned param_cnt() const override { return 1; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        for (unsigned i = 0; i < in_p.size(); ++i)
        {
          jac[i] = (i == (unsigned)param[PARAM_IDX]) ? 1 : 0;
        }
      }
      return in_p[(unsigned)param[PARAM_IDX]];
    }
  };

  class ConstParamNode : public ParamNode
  {
    static constexpr int CONST_COEFF = 0;
  public:
    ConstParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "const"; }
    virtual unsigned source_cnt() const override { return 0; }
    virtual unsigned param_cnt() const override { return 1; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
        }
      }
      return param[CONST_COEFF];
    }
  };

  void AddParamNode_apply(const float *in, float *out)
  {
    out[0] = in[0] + in[1];
  }

  class AddParamNode : public ParamNode
  {
  public:
    AddParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "add"; }
    virtual unsigned source_cnt() const override { return 2; }
    virtual unsigned param_cnt() const override { return 0; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, 2);
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        ENZYME_EVALUATE_WITH_DIFF(AddParamNode_apply, 2, 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      AddParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  void SubParamNode_apply(const float *in, float *out)
  {
    out[0] = in[0] - in[1];
  }

  class SubParamNode : public ParamNode
  {
  public:
    SubParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "sub"; }
    virtual unsigned source_cnt() const override { return 2; }
    virtual unsigned param_cnt() const override { return 0; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, source_cnt());
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        ENZYME_EVALUATE_WITH_DIFF(SubParamNode_apply, source_cnt(), 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      SubParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  void MultParamNode_apply(const float *in, float *out)
  {
    out[0] = in[0] * in[1];
  }

  class MultParamNode : public ParamNode
  {
  public:
    MultParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "mult"; }
    virtual unsigned source_cnt() const override { return 2; }
    virtual unsigned param_cnt() const override { return 0; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, source_cnt());
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        ENZYME_EVALUATE_WITH_DIFF(MultParamNode_apply, source_cnt(), 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      MultParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  void DivParamNode_apply(const float *in, float *out)
  {
    out[0] = in[0] / in[1];
  }

  class DivParamNode : public ParamNode
  {
  public:
    DivParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "div"; }
    virtual unsigned source_cnt() const override { return 2; }
    virtual unsigned param_cnt() const override { return 0; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, source_cnt());
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        ENZYME_EVALUATE_WITH_DIFF(DivParamNode_apply, source_cnt(), 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      DivParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  class IncParamNode : public ParamNode
  {
    static constexpr int INC_COEFF = 0;
  public:
    IncParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "inc"; }
    virtual unsigned source_cnt() const override { return 1; }
    virtual unsigned param_cnt() const override { return 1; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, source_cnt());
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        x.push_back(param[INC_COEFF]);
        ENZYME_EVALUATE_WITH_DIFF(AddParamNode_apply, source_cnt(), 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      x.push_back(param[INC_COEFF]);
      AddParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  class CMultParamNode : public ParamNode
  {
    static constexpr int MULT_COEFF = 0;
  public:
    CMultParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "mult"; }
    virtual unsigned source_cnt() const override { return 1; }
    virtual unsigned param_cnt() const override { return 1; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, source_cnt());
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        x.push_back(param[MULT_COEFF]);
        ENZYME_EVALUATE_WITH_DIFF(MultParamNode_apply, source_cnt(), 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      x.push_back(param[MULT_COEFF]);
      MultParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  void NegParamNode_apply(const float *in, float *out)
  {
    out[0] = -in[0];
  }

  class NegParamNode : public ParamNode
  {
  public:
    NegParamNode(unsigned id, std::vector<float> p) : ParamNode(id, p) { param = p; sources = {}; name = "div"; }
    virtual unsigned source_cnt() const override { return 1; }
    virtual unsigned param_cnt() const override { return 0; }
    virtual float get_res(float *jac) const override
    {
      if (jac != NULL)
      {
        upg::UniversalGenJacobian tmp;
        tmp.resize(1, source_cnt());
        float ans;
        std::vector<float> x(0);
        std::vector<std::vector<float>> jacs(source_cnt());
        for (int i = 0; i < source_cnt(); ++i)
        {
          jacs[i].resize(in_p.size());
          x.push_back(sources[i]->get_res(jacs[i].data()));
        }
        ENZYME_EVALUATE_WITH_DIFF(NegParamNode_apply, source_cnt(), 1, x.data(), &ans, tmp.data());
        for (int i = 0; i < in_p.size(); ++i)
        {
          jac[i] = 0;
          for (int j = 0; j < source_cnt(); ++j)
          {
            jac[i] += tmp.at(j, 0) * jacs[j][i];
          }
        }
        return ans;
      }
      float ans;
      std::vector<float> x(0);
      for (int i = 0; i < source_cnt(); ++i)
      {
        x.push_back(sources[i]->get_res(NULL));
      }
      NegParamNode_apply(x.data(), &ans);
      return ans;
    }
  };

  ParamNode *param_node_by_node_type_id(unsigned id, uint16_t num, std::vector <float> p)
  {
    ParamNode *node = NULL;
    switch(num)
    {
      case ParamNode::PRIMITIVE: 
        node = new PrimitiveParamNode(id, p);
        break;
      case ParamNode::ADD: 
        node = new AddParamNode(id, p);
        break;
      case ParamNode::SUB: 
        node = new SubParamNode(id, p);
        break;
      case ParamNode::MULT: 
        node = new MultParamNode(id, p);
        break;
      case ParamNode::DIV: 
        node = new DivParamNode(id, p);
        break;
      case ParamNode::INC: 
        node = new IncParamNode(id, p);
        break;
      case ParamNode::CMULT: 
        node = new CMultParamNode(id, p);
        break;
      case ParamNode::CONST: 
        node = new ConstParamNode(id, p);
        break;
      case ParamNode::NEG: 
        node = new NegParamNode(id, p);
        break;
      default:
        logerr("invalid node_type %u\n",id);
        node = nullptr;
        break;
    }
    return node;
  }
}