#include "gen_params.h"
#include "models.h"

void add_rect(u_g::vec3 point, u_g::vec3 v1, u_g::vec3 v2, SimpleMeshData &mesh);
void add_tri(u_g::vec3 point, u_g::vec3 v1, u_g::vec3 v2, SimpleMeshData &mesh);

class GenNode
{
public:
  virtual SimpleMeshData apply(Params &p) = 0;
};

class PrimitiveNode : protected GenNode
{
  
};

class FigureNode : protected PrimitiveNode
{
  SimpleMeshData apply(Params &p) override
  {
    my_float size = p.get();
    SimpleMeshData mesh;
    //creating for example cube
    add_rect({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, mesh);
    add_rect({0, 0, 0}, {0, 0, 1}, {0, 1, 0}, mesh);
    add_rect({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, mesh);

    add_rect({1, 1, 1}, {-1, 0, 0}, {0, -1, 0}, mesh);
    add_rect({1, 1, 1}, {0, -1, 0}, {0, 0, -1}, mesh);
    add_rect({1, 1, 1}, {0, 0, -1}, {-1, 0, 0}, mesh);
    return mesh;
  }
};

template <int N>
class SpinNode : protected PrimitiveNode
{
  SimpleMeshData apply(Params &p) override
  {
    my_float data[N];
    for (int i = 0; i < N; ++i)
    {
      data[i] = p.get();
    }
    SimpleMeshData mesh;
    //creating spin mesh

    return mesh;
  }
};

class TransformNode : protected GenNode
{
  
};

class OneChildNode : protected TransformNode
{
protected:
  GenNode *child;
};

class ScaleNode : protected OneChildNode
{
  SimpleMeshData apply(Params &p) override
  {
    my_float scale_x = p.get();
    my_float scale_y = p.get();
    my_float scale_z = p.get();

    SimpleMeshData mesh = child->apply(p);
    //applying scaling
    for (int i = 0; i + mesh.elem_size <= mesh.data.size(); i += mesh.elem_size)
    {
      mesh.data[i] *= scale_x;
      mesh.data[i + 1] *= scale_y;
      mesh.data[i + 2] *= scale_z;
    }
    return mesh;
  }
};

class MoveNode : protected OneChildNode
{
  SimpleMeshData apply(Params &p) override
  {
    my_float move_x = p.get();
    my_float move_y = p.get();
    my_float move_z = p.get();
    SimpleMeshData mesh = child->apply(p);
    //applying moving
    for (int i = 0; i + mesh.elem_size <= mesh.data.size(); i += mesh.elem_size)
    {
      mesh.data[i] += move_x;
      mesh.data[i + 1] += move_y;
      mesh.data[i + 2] += move_z;
    }
    return mesh;
  }
};

class RotateNode : protected OneChildNode
{
  SimpleMeshData apply(Params &p) override
  {
    my_float axis_x = p.get();
    my_float axis_y = p.get();
    my_float axis_z = p.get();
    my_float angle = p.get();

    SimpleMeshData mesh = child->apply(p);
    //applying rotating
    return mesh;
  }
};

class TwoChildNode : protected TransformNode
{
protected:
  GenNode *left;
  GenNode *right;
};

class AndNode : protected TwoChildNode
{
  SimpleMeshData apply(Params &p) override
  {
    SimpleMeshData mesh1 = left->apply(p);
    SimpleMeshData mesh2 = right->apply(p);
    //applying and
    return mesh1;
  }
};

class OrNode : protected TwoChildNode
{
  SimpleMeshData apply(Params &p) override
  {
    SimpleMeshData mesh1 = left->apply(p);
    SimpleMeshData mesh2 = right->apply(p);
    //applying or
    SimpleMeshData mesh;
    mesh.data = mesh1.data;
    mesh.data.insert(mesh.data.end(), mesh2.data.begin(), mesh2.data.end());
    return mesh1;
  }
};

class SubtrNode : protected TwoChildNode
{
  SimpleMeshData apply(Params &p) override
  {
    SimpleMeshData mesh1 = left->apply(p);
    SimpleMeshData mesh2 = right->apply(p);
    //applying subtract
    return mesh1;
  }
};