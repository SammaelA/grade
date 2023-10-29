#include "generation.h"
#include "reconstruction_impl.h"
#include "tinyEngine/engine.h"

namespace upg
{

  void mesh_to_complex_model(const UniversalGenMesh &mesh, ComplexModel &mod)
  {
    assert(mesh.pos.size()%9 == 0);
    assert(mesh.pos.size()/3 == mesh.norm.size()/3);
    assert(mesh.pos.size()/3 == mesh.tc.size()/2);

    mod.models.push_back(new Model());
    Model *m = mod.models.back();
    m->positions = mesh.pos;
    m->normals = mesh.norm;
    int sz = mesh.pos.size()/3;//number of vertices
    m->colors.resize(4*sz);
    for (int i=0;i<sz;i++)
    {
      m->colors[4*i]   = mesh.tc[2*i];
      m->colors[4*i+1] = mesh.tc[2*i+1];
      m->colors[4*i+2] = 0;
      m->colors[4*i+3] = 1;
    }

    m->indices.resize(3*sz);
    for (int i=0;i<sz;i++)
      m->indices[i] = i;

    //Some generic texture. You can choose another one (see resources.blk for available options)
    mod.materials.push_back(Material(engine::textureManager->get("porcelain")));
  }

  bool create_model_from_block(Block &bl, ComplexModel &mod)
  {
    UPGStructure structure;
    UPGParametersRaw params;
    bl.get_arr("structure", structure.s);
    bl.get_arr("params", params.p);

    //create mesh here
    //UniversalGenInstance gen(structure);
    //auto mesh = gen.generate(params.p);
    UniversalGenMesh mesh;
    mesh.pos = {0,0,0, -1,0,0, 0,1,-1};
    mesh.norm = {0,1/sqrtf(2),1/sqrtf(2), 0,1/sqrtf(2),1/sqrtf(2), 0,1/sqrtf(2),1/sqrtf(2)};
    mesh.tc = {0,0, 1,0, 0,1};
    mesh_to_complex_model(mesh, mod);

    return true;
  }
}