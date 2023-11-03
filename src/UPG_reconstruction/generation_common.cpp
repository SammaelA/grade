#include "generation_common.h"
#include "tinyEngine/engine.h"
#include "upg.h"
#include "generation.h"

namespace upg
{
  void mesh_to_complex_model(const UniversalGenMesh &mesh, ComplexModel &mod);
  void add_tri_data(upg::vec3 point, upg::vec3 n, upg::vec2 tex, UniversalGenMesh &mesh)
  {
    mesh.pos.push_back(point.x);
    mesh.pos.push_back(point.y);
    mesh.pos.push_back(point.z);
    mesh.norm.push_back(n.x);
    mesh.norm.push_back(n.y);
    mesh.norm.push_back(n.z);
    mesh.tc.push_back(tex.x);
    mesh.tc.push_back(tex.y);
  }
  void add_point_data(upg::vec3 point, UniversalGenMesh &mesh)
  {
    mesh.pos.push_back(point.x);
    mesh.pos.push_back(point.y);
    mesh.pos.push_back(point.z);
  }


  void mesh_to_complex_model(const UniversalGenMesh &mesh, ComplexModel &mod)
  {
    assert(mesh.pos.size()%9 == 0);
    //assert(mesh.pos.size()/3 == mesh.norm.size()/3);
    //assert(mesh.pos.size()/3 == mesh.tc.size()/2);

    mod.models.push_back(new Model());
    Model *m = mod.models.back();
    m->positions = mesh.pos;
    if (mesh.norm.empty())
    {
      m->normals.resize(mesh.pos.size());
      for (int i=0;i<mesh.pos.size(); i+=9)
      {
        glm::vec3 p1(mesh.pos[i], mesh.pos[i+1], mesh.pos[i+2]);
        glm::vec3 p2(mesh.pos[i+3], mesh.pos[i+4], mesh.pos[i+5]);
        glm::vec3 p3(mesh.pos[i+6], mesh.pos[i+7], mesh.pos[i+8]);
        glm::vec3 n = glm::normalize(glm::cross(p2 - p1, p3 - p1));
        for (int j=0;j<3;j++)
        {
          m->normals[i+3*j] = n.x;
          m->normals[i+3*j+1] = n.y;
          m->normals[i+3*j+2] = n.z;
        }
      }
    }
    else
    {
      assert(mesh.pos.size() == mesh.norm.size());
      m->normals = mesh.norm;
    }
    int sz = mesh.pos.size()/3;//number of vertices
    m->colors.resize(4*sz);
    if (mesh.tc.empty())
    {
      float tc[12] = {0,0,0,1, 1,0,0,1, 0,1,0,1};
      for (int i=0;i<m->colors.size();i++)
        m->colors[i] = tc[i % 12];
    }
    else
    {
      assert(mesh.pos.size()/3 == mesh.tc.size()/2);
      for (int i=0;i<sz;i++)
      {
        m->colors[4*i]   = mesh.tc[2*i];
        m->colors[4*i+1] = mesh.tc[2*i+1];
        m->colors[4*i+2] = 0;
        m->colors[4*i+3] = 1;
      }
    }

    m->indices.resize(3*sz);
    for (int i=0;i<sz;i++)
      m->indices[i] = i;

    //Some generic texture. You can choose another one (see resources.blk for available options)
    mod.materials.push_back(Material(engine::textureManager->get("porcelain")));
  }

  bool create_model_from_block(const Block &bl, ComplexModel &mod)
  {
    UPGStructure structure;
    UPGParametersRaw params;
    bl.get_arr("structure", structure.s);
    bl.get_arr("params", params.p);

    return create_model(structure, params, mod);
  }

  bool create_model(const UPGStructure &structure, const UPGParametersRaw &params,
                    ComplexModel &mod)
  {
    logerr("paramss len %d", params.p.size());
    //create mesh here
    UniversalGenInstance gen(structure);
    auto mesh = gen.generate(params.p);
    //mesh.pos = {0,0,0, -1,0,0, 0,1,-1};
    //mesh.norm = {0,1/sqrtf(2),1/sqrtf(2), 0,1/sqrtf(2),1/sqrtf(2), 0,1/sqrtf(2),1/sqrtf(2)};
    //mesh.tc = {0,0, 1,0, 0,1};
    mesh_to_complex_model(mesh, mod);

    return true;
  }
}