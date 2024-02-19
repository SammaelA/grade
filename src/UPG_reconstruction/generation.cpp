#include "generation.h"
#include "tinyEngine/engine.h"

namespace upg
{
  UniversalGenMesh MeshGenInstance::generate(std::span<const float> parameters, UniversalGenJacobian *jac)
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
    return root->apply(jac);
    // generator.take_params(parameters);
    // return generator.generate();
  }

  MeshGenInstance::MeshGenInstance(const UPGStructure &structure)
  {
    recreate(structure);
  }
  
  void MeshGenInstance::recreate(const UPGStructure &structure)
  {
    all_params.clear();
    all_nodes.clear();

    std::vector<GenNode *> nodes;
    std::vector<std::pair<GenNode *, unsigned>> param_startings;
    int i = 0;
    do
    {
      unsigned n = 0;
      if (i < structure.s.size())
      {
        n = structure.s[i];
      }
      GenNode *node = node_by_node_type_id(n, i);
      all_nodes.push_back(std::unique_ptr<GenNode>(node));
      desc.add_parameters(node->get_ID(), node->get_node_name(), node->get_parameters_block());
      param_startings.push_back({node, all_params.size()});
      all_params.resize(all_params.size() + node->param_cnt());
      
      if (i == 0)
      {
        root = node;
        if (node->child_cnt() > 0)
        {
          nodes.push_back(node);
        }
      }
      else
      {
        GenNode *last = nodes[nodes.size() - 1];
        
        if (!last->add_child(node))
        {
          nodes.pop_back();
        }
        if (node->child_cnt() > 0)
        {
          nodes.push_back(node);
        }
      }
      ++i;
    } while (nodes.size() > 0);

    int offset = 0;
    for (auto &nptr : all_nodes)
    {
      nptr->set_param_span(std::span<float>(all_params.data() + offset, nptr->param_cnt()));
      offset += nptr->param_cnt();
    }
  }
}