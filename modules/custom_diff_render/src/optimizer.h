#pragma once
#include "dmesh.h"
#include "drender.h"
#include <memory> // for shared pointers

namespace diff_render
{
struct OptimizerParameters
{ 
  enum OPT_ALGORITHM {GD_Naive=1, GD_AdaGrad=2, GD_RMSProp=3, GD_Adam=4};
  OptimizerParameters() = default;
  OptimizerParameters(OPT_ALGORITHM _alg)
  {
    alg = _alg;
    set_default();
  }
  OptimizerParameters(OPT_ALGORITHM _alg, float pos_lr, float tex_lr, float _transforms_lr = 0.1, 
                      ::std::vector<int> diff_mesh_ids = {0},
                      bool _verbose = false)
  {
    alg = _alg;
    set_default();
    position_lr = pos_lr;
    textures_lr = tex_lr;
    transforms_lr = _transforms_lr;
    verbose = _verbose;
    differentiable_mesh_ids = diff_mesh_ids;
  }
  OPT_ALGORITHM alg = GD_Naive;
  int decayPeriod   = 25;
  float decay_mult  = 0.75;
  float base_lr     = 0.1;
  float position_lr = 0.1;
  float textures_lr = 0.1;
  float transforms_lr = 0.1;
  ::std::vector<int> differentiable_mesh_ids = {0};
  bool verbose = true;
private:
  void set_default();
};

struct IOptimizer
{
  virtual ~IOptimizer() {};
  virtual void Init(const Scene& a_scene, ::std::shared_ptr<IDiffRender> a_pDRImpl, 
                    const CamInfo* a_cams, const Img* a_images, int a_numViews, OptimizerParameters a_params) = 0;

  virtual Scene Run (size_t a_numIters, float &final_error, ::std::vector<Scene> *iterations_dump = nullptr) = 0;
};

extern IOptimizer* CreateSimpleOptimizer();
}