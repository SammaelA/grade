#pragma once
#include "reconstruction.h"
#include "tree_node.h"
#include <memory>
namespace upg
{  
  //model that is produced by universal generator
  //it's should be optimal to use for reconstruction
  //
  /*struct UniversalGenMesh
  {
    //triangle mesh pos.size()%9 == 0
    //norm and tc can be empty
    std::vector<float> pos; //vec3
    std::vector<float> norm;//vec3
    std::vector<float> tc; //vec2
  };*/

  //Structure that represents jacobian dpos/dP
  //where P is parameters list for specific
  //UniversalGenInstance. Jacobian is likely
  //to be represented as sprase or block matrix
  struct UniversalGenJacobian
  {
    int x_n, y_n;
    std::vector<float> jacobian;
  };

  //Generator with fixed structure.
  //It contains nodes tree as well as description for parameters that
  //should be passed to generate() function
  class UniversalGenInstance
  {
  public:
    UniversalGenInstance(const UPGStructure &structure);
    //~UniversalGenInstance();
    UniversalGenMesh generate(std::span<const float> parameters);
    //It was used for testing
    //{
    //  UniversalGenMesh mesh;
    //  for (int i=0;i<9;i++)
    //    mesh.pos.push_back(parameters[i]);
    //  return mesh;
    //}
    UniversalGenJacobian generate_jacobian(std::span<const float> parameters);
    ParametersDescription desc;

  private:
    std::vector<std::unique_ptr<GenNode>> all_nodes;
    GenNode *root;
    //spans from inputParams points to this container
    //put raw parameters list here to generate
    //DO NOT change size of this vector
    std::vector<my_float> all_params;
    void tree_del();
  };

  //Interface for Universal generator in general
  //Use if you just want to create a specific model
  //and know both structure and parameters
  class UniversalGen
  {
  public:
    UniversalGenMesh generate(const UPGStructure &structure, std::span<const float> parameters)
    {
      UniversalGenInstance gen(structure);
      return gen.generate(parameters);
    }
  };
}