#pragma once
#include "upg.h"
#include "tree_node.h"
#include "generation_common.h"
#include <memory>
namespace upg
{  
  //Structure that represents jacobian dpos/dP
  //where P is parameters list for specific
  //UniversalGenInstance. Jacobian is likely
  //to be represented as sprase or block matrix
  class UniversalGenJacobian
  {
  public:
    void resize(int xn, int yn)
    {
      x_n = xn;
      y_n = yn;
      jacobian.resize(xn*yn);
    }
    int get_xn() const 
    {
      return x_n;
    }
    int get_yn() const 
    {
      return y_n;
    }
    float &at(int y, int x)
    {
      return jacobian[y*x_n + x];
    }
    const float &at(int y, int x) const
    {
      return jacobian[y*x_n + x];
    }
    float *data()
    {
      return jacobian.data();
    }
    void clear()
    {
      std::fill(jacobian.begin(), jacobian.end(), 0);
    }
  private:
    int x_n, y_n;
    std::vector<float> jacobian;
  };

  struct UniversalGenMesh
  {
    //triangle mesh pos.size()%9 == 0
    //norm and tc can be empty
    std::vector<float> pos; //vec3
    std::vector<float> norm; //vec3
    std::vector<float> tc; //vec2
  };


  //Generator with fixed structure.
  //It contains nodes tree as well as description for parameters that
  //should be passed to generate() function
  class MeshGenInstance : public UniversalGenInstance
  {
  public:
    MeshGenInstance(const UPGStructure &structure);
    virtual void recreate(const UPGStructure &structure) override;
    UniversalGenMesh generate(std::span<const float> parameters, UniversalGenJacobian *jac = nullptr);
    ParametersDescription desc;

  private:
    std::vector<std::unique_ptr<GenNode>> all_nodes;
    GenNode *root;
    //spans from inputParams points to this container
    //put raw parameters list here to generate
    //DO NOT change size of this vector
    std::vector<my_float> all_params;
  };
}