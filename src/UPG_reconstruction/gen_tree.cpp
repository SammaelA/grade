#include "gen_tree.h"

namespace upg
{
  GenNode * Tree::node_by_number(uint16_t num, unsigned id)
  {
    GenNode *node = NULL;
    switch(num)
    {
      case 2:
        node = new ScaleNode(id);
        break;
      case 3:
        node = new MoveNode(id);
        break;
      case 4:
        node = new RotateNode(id);
        break;
      case 5:
        node = new AndNode(id);
        break;
      case 6:
        node = new OrNode(id);
        break;
      case 7:
        node = new SubtrNode(id);
        break;
      default:
        node = new FigureNode(id);
        break;

    }
    return node;
  }

  void Tree::create(const UPGStructure &structure)
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
      GenNode *node = node_by_number(n, i);
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

    
    for (auto i : param_startings)
    {
      i.first->set_param_span({all_params.data() + i.second, i.first->param_cnt()});
    }
  }

  void Tree::take_params(std::span<const float> param)
  {
    for (int i = 0; i < all_params.size(); ++i)
    {
      if (i < param.size())
      {
        all_params[i] = param[i];
      }
      else
      {
        all_params[i] = 0;
      }
    }
  }

  void Tree::tree_del()
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