#include "gen_params.h"
#include "models.h"

template <class Model>
class GenNode
{
public:
  virtual Model apply(Params &p) = 0;
};

template <class Model>
class PrimitiveNode : protected GenNode
{
  
};

template <class Model>
class FigureNode : protected PrimitiveNode
{
  Model apply(Params &p) override
  {
    my_float size = p.get();
    Model mesh;
    //creating for example cube
    return mesh;
  }
};

template <class Model, int N>
class SpinNode : protected PrimitiveNode
{
  Model apply(Params &p) override
  {
    my_float data[N];
    for (int i = 0; i < N; ++i)
    {
      data[i] = p.get();
    }
    Model mesh;
    //creating spin mesh
    return mesh;
  }
};

template <class Model>
class TransformNode : protected GenNode
{
  
};

template <class Model>
class OneChildNode : protected TransformNode
{
protected:
  GenNode *child<Model>;
};

template <class Model>
class ScaleNode : protected OneChildNode
{
  Model apply(Params &p) override
  {
    my_float scale_x = p.get();
    my_float scale_y = p.get();
    my_float scale_z = p.get();

    Model mesh = child->apply(p);
    //applying scaling
    return mesh;
  }
};

template <class Model>
class MoveNode : protected OneChildNode
{
  Model apply(Params &p) override
  {
    my_float move_x = p.get();
    my_float move_y = p.get();
    my_float move_z = p.get();
    Model mesh = child->apply(p);
    //applying moving
    return mesh;
  }
};

template <class Model>
class RotateNode : protected OneChildNode
{
  Model apply(Params &p) override
  {
    my_float axis_x = p.get();
    my_float axis_y = p.get();
    my_float axis_z = p.get();
    my_float angle = p.get();

    Model mesh = child->apply(p);
    //applying rotating
    return mesh;
  }
};

template <class Model>
class TwoChildNode : protected TransformNode
{
protected:
  GenNode *left<Model>;
  GenNode *right<Model>;
};

template <class Model>
class AndNode : protected TwoChildNode
{
  Model apply(Params &p) override
  {
    Model mesh1 = left->apply(p);
    Model mesh2 = right->apply(p);
    //applying and
    Model mesh = mesh1 & mesh2;
    return mesh;
  }
};

template <class Model>
class OrNode : protected TwoChildNode
{
  Model apply(Params &p) override
  {
    Model mesh1 = left->apply(p);
    Model mesh2 = right->apply(p);
    //applying or
    Model mesh = mesh1 | mesh2;
    return mesh;
  }
};

template <class Model>
class SubtrNode : protected TwoChildNode
{
  Model apply(Params &p) override
  {
    Model mesh1 = left->apply(p);
    Model mesh2 = right->apply(p);
    //applying subtract
    Model mesh = mesh1 - mesh2;
    return mesh;
  }
};