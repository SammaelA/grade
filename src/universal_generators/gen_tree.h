#include "tree_node.h"
template <class Model>
class Tree
{
  GenNode<Model> *root;
public:
  Model generate(const Params &p)
  {
    Params param = p;
    return root->apply(param);
  }
  void read_blk(/*blk file*/)
  {

  }
  void write_blk(/*blk file*/)
  {

  }
};