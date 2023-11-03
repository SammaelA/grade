#include "generation_impl.h"

namespace upg
{
  UniversalGenMesh UniversalGenInstance::generate(std::span<const float> parameters)
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
    return root->apply();
    // generator.take_params(parameters);
    // return generator.generate();
  }

  UniversalGenInstance::UniversalGenInstance(const UPGStructure &structure)
  {
        all_params.clear();
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
  }

  UniversalGenJacobian UniversalGenInstance::generate_jacobian(std::span<const float> parameters) // maybe it will create mesh too?
  {
    UniversalGenJacobian jac;
    jac.x_n = 9;
    jac.y_n = 9;
    jac.jacobian = std::vector<float>(9 * 9, 0);
    for (int i = 0; i < jac.y_n; i++)
      jac.jacobian[i * jac.x_n + i] = 1;

    return jac;
  }

  void UniversalGenInstance::tree_del()
  {
    std::vector<GenNode *> nodes;
    nodes.push_back(root);
    while (nodes.size() > 0)
    {
      GenNode *last = nodes[nodes.size() - 1];
      nodes.pop_back();
      if (last != NULL)
      {
        nodes.insert(nodes.end(), last->childs().begin(), last->childs().end());
        delete last;
      }
    }
  }
}