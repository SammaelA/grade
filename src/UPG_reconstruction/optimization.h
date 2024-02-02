#pragma once
#include "upg.h"
#include "generation_common.h"
#include <vector>
#include <memory>

namespace upg
{
  // all parameters that can be changed by optimizer structured
  // in a convenient (for optimizer) way
  // This class is basically a vector that is deliberately put into
  // a separate structure to clarify it's meaning (i.e. - it should not be put 
  // into a generator directly!)
  struct OptParams
  {
    int size() const
    {
      return data.size();
    }
    void resize(int size)
    {
      data.resize(size);
    }
    float &operator[](int index)
    {
      return data[index];
    }
    const float &operator[](int index) const
    {
      return data[index];
    }
    void set(const std::span<const float> &p)
    {
      data = std::vector<float>(p.begin(), p.end());
    }
  private:
    std::vector<float> data;
  };

  /*
    UPGOptimizableFunction class encapsulates a bunch of things that optimizer shouldn't be aware of:
    1. What data structure is used to represent reference
    2. What method is used to actually calculate the function
    3. what parameters are used to represent cameras (if we even care about cameras)
    However, it cannot be treated as an arbitrary function, because optimizer still knows about UPG,
    structure and other stuff.
  */
  class UPGOptimizableFunction
  {
  public:
    virtual ~UPGOptimizableFunction() {};
    static std::vector<float> opt_params_to_gen_params(const OptParams &params, const ParametersDescription &pd,
                                                       std::vector<float> *out_scene_params = nullptr)
    {
      //iterate all parameters' groups from description
      //map orders them by block_id, so generator's params are first, and scene's are after it
      int i = 0;
      std::vector<float> gen_params;
      std::vector<float> scene_params;
      gen_params.reserve(pd.get_total_params_count());
      for (const auto &p : pd.get_block_params())
      {
        bool is_scene_block = p.first >= get_scene_params_block_offset();
        std::vector<float> &p_v = is_scene_block ? scene_params : gen_params;
        for (auto &par : p.second.p)
        {
          if (par.type == ParameterType::CONST)
          {
            p_v.push_back(par.value);
          }
          else
          {
            float v = CLAMP(params[i], par.min_val, par.max_val);
            p_v.push_back(v);
            i++;
          }
        }
      }
      if (out_scene_params)
        *out_scene_params = scene_params;

      return gen_params;
    }
    static OptParams gen_params_to_opt_params(const std::vector<float> &params, const ParametersDescription &pd)
    {
      OptParams p;
      std::vector<float> non_const_params;
      non_const_params.reserve(params.size());

      int p_id = 0;
      for (const auto &p : pd.get_block_params())
      {
        bool is_scene_block = p.first >= get_scene_params_block_offset();
        for (auto &par : p.second.p)
        {
          if (par.type != ParameterType::CONST)
            non_const_params.push_back(params[p_id]);
          p_id++;
        }
      }

      p.set(non_const_params);
      return p;
    }
    //calculate function that we optimize and it's gradient (put into given span)
    //requires already created UniversalGenInstance
    //size of out_grad - is a number of differentiable parameters in params
    virtual float f_grad_f(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params, 
                           std::span<float> out_grad) = 0;
    virtual float f_no_grad(UniversalGenInstance *gen, const ParametersDescription &pd, const OptParams &params) = 0;

    //based on GenInstance create full parameters description, that includes both generator's and 
    //scene's parameter blocks
    virtual ParametersDescription get_full_parameters_description(const UniversalGenInstance *gen) const = 0;
    virtual std::shared_ptr<UniversalGenInstance> get_generator(const UPGStructure &structure) const = 0;
    virtual float estimate_positioning_quality(const UPGStructure &structure,
                                               const UPGPart &part, std::span<const float> parameters,
                                               float border_sigma = 0.01f,
                                               float inner_point_penalty = 10) const 
    {
      return 0;
    }
  protected:
    //we prepare parameters descriptions for generator's and scene (i.e. camera) parameters
    //independently, so we must assure that their block ids do not overlap. To do so,
    //add offset to all parametrs block related to scene
    static unsigned get_scene_params_block_offset() { return 1<<16; }
  };

  class UPGOptimizer
  {
  public:
    virtual ~UPGOptimizer() = default;
    virtual void optimize(int iters = -1) = 0; //iters = -1 mean that we take it from settings
    virtual std::vector<UPGReconstructionResult> get_best_results() = 0;
    UPGOptimizer(UPGOptimizableFunction *_func):
    func(_func)
    {

    }
  protected:
    UPGOptimizableFunction *func;
  };

  std::shared_ptr<UPGOptimizer> get_optimizer_adam(UPGOptimizableFunction *_func, 
                                                   const Block &settings, const UPGReconstructionResult &start_params);
  std::shared_ptr<UPGOptimizer> get_optimizer_memetic(UPGOptimizableFunction *_func, 
                                                      const Block &settings, const UPGStructure &_structure);
  std::shared_ptr<UPGOptimizer> get_optimizer_CHC(UPGOptimizableFunction *_func, 
                                                  const Block &settings, const UPGStructure &_structure);
  std::shared_ptr<UPGOptimizer> get_optimizer_particle_swarm(UPGOptimizableFunction *_func, 
                                                             const Block &settings, const UPGStructure &_structure);
  std::shared_ptr<UPGOptimizer> get_optimizer_CC(UPGOptimizableFunction *_func, 
                                                 const Block &settings, const UPGStructure &_structure);
  std::shared_ptr<UPGOptimizer> get_optimizer_differentiable_evolution(UPGOptimizableFunction *_func, 
                                                                       const Block &settings, const UPGStructure &_structure);
  std::shared_ptr<UPGOptimizer> get_optimizer_iterative_fitting(UPGOptimizableFunction *_func, 
                                                                const Block &settings, const UPGStructure &_structure);
  std::shared_ptr<UPGOptimizer> get_optimizer_multistep_adam(UPGOptimizableFunction *_func, 
                                                             const Block &settings, const UPGReconstructionResult &start_params);

  std::shared_ptr<UPGOptimizer> optimization_test_stand(UPGOptimizableFunction *_func, const Block &settings, 
                                                        const UPGStructure &_structure, const UPGParametersRaw &_params);
}