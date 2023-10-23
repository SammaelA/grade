#include "tree_node.h"

namespace u_g
{
  class Tree
  {
    GenNode *root;
  public:
    SimpleMeshData generate(const Params &p)
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
}