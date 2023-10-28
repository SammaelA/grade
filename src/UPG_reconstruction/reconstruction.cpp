#include "reconstruction.h"
#include "preprocessing.h"
#include "reconstruction_impl.h"
#include "preprocessing.h"
#include "graphics_utils/simple_model_utils.h"
#include "tinyEngine/engine.h"

namespace upg
{
  std::vector<ParametersDescription::Param> get_camera_params(const CameraSettings &cam)
  {
    std::vector<ParametersDescription::Param> params;

    params.push_back({cam.origin.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.x"});
    params.push_back({cam.origin.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.y"});
    params.push_back({cam.origin.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "origin.z"});

    params.push_back({cam.target.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.x"});
    params.push_back({cam.target.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.y"});
    params.push_back({cam.target.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "target.z"});

    params.push_back({cam.up.x, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.x"});
    params.push_back({cam.up.y, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.y"});
    params.push_back({cam.up.z, 0.0f, 1000.0f, ParameterType::DIFFERENTIABLE, "up.z"});

    params.push_back({cam.z_near, 0.001f, 1.0f, ParameterType::MUTABLE_FLOAT, "z_near"});
    params.push_back({cam.z_far, 10.0f, 1000.0f, ParameterType::MUTABLE_FLOAT, "z_far"});
    params.push_back({cam.fov_rad, 0.1f, 1.5f, ParameterType::MUTABLE_FLOAT, "fov_rad"});

    return params;
  }

  ParametersDescription get_cameras_parameter_description(const Block &cameras_blk)
  {
    ParametersDescription pd;
    static constexpr unsigned camera_block_id = 1 << 16;
    assert(cameras_blk.size() >= 1);
    for (int i = 0; i < cameras_blk.size(); i++)
    {
      Block *cam_blk = cameras_blk.get_block(i);
      if (!cam_blk)
      {
        logerr("UPG Reconstruction: views block should contain only blocks with data for each view");
        continue;
      }
      CameraSettings cam = visualizer::load_camera_settings(*cam_blk);
      auto pv = get_camera_params(cam);
      if (cam_blk->get_bool("camera.fixed", true))
      {
        for (auto &p : pv)
          p.type = ParameterType::CONST;
      }
      pd.add_parameters(camera_block_id + (unsigned)i, "camera_" + std::to_string(i), pv);
    }
    return pd;
  }

  std::vector<ReferenceView> get_reference(const Block &input_blk)
  {
    std::vector<ReferenceView> reference;
    for (int i = 0; i < input_blk.size(); i++)
    {
      Block *cam_blk = input_blk.get_block(i);
      if (!cam_blk)
      {
        logerr("UPG Reconstruction: views block should contain only blocks with data for each view");
        continue;
      }
      reference.push_back(preprocess_get_reference_view(*cam_blk));
    }

    return reference;
  }

  //all data created and used during the optimization process
  //such as population, best_structures and parameter sets
  struct UPGOptimizationContext
  {

  };

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk)
  {
    //we store settings in separate block for each stage of the algorithm
    Block *input_blk = blk.get_block("input");
    Block *gen_blk = blk.get_block("generator");
    Block *opt_blk = blk.get_block("optimization");
    if (!input_blk || !gen_blk || !opt_blk)
    {
      logerr("UPG Reconstruction: input, generator, optimization blocks should exist in configuration");
      return {};
    }

    ParametersDescription cameras_pd = get_cameras_parameter_description(*input_blk);
    std::vector<ReferenceView> reference = get_reference(*input_blk);

    std::vector<UPGReconstructionResult> res;

    cameras_pd.print_info();
    return res;
  }

  //TODO: move it somewhere from here
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