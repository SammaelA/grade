#include "reconstruction.h"
#include "preprocessing.h"
#include "reconstruction_impl.h"
#include "preprocessing.h"
#include "graphics_utils/simple_model_utils.h"

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
    static constexpr unsigned camera_block_id = 1<<16;
    assert(cameras_blk.size()>=1);
    for (int i=0; i<cameras_blk.size(); i++)
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
      pd.add_parameters(camera_block_id + (unsigned)i, "camera_"+std::to_string(i), pv);
    }
    return pd;
  }

  std::vector<ReferenceView> get_reference(const Block &input_blk)
  {
    std::vector<ReferenceView> reference;
    for (int i=0; i<input_blk.size(); i++)
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

  std::vector<UPGReconstructionResult> reconstruct(const Block &blk)
  {
    Block *input_blk = blk.get_block("input");
    Block *gen_blk = blk.get_block("generator");
    Block *opt_blk = blk.get_block("optimization");
    if (!input_blk || !gen_blk || !opt_blk)
    {
      logerr("UPG Reconstruction: input, generator, optimization blocks should exist in configuration");
      return {};
    }

    ParametersDescription cameras_pd = get_cameras_parameter_description(*input_blk);
    cameras_pd.print_info();

    std::vector<ReferenceView> reference = get_reference(*input_blk);

    std::vector<UPGReconstructionResult> res;
    return res;
  }
}