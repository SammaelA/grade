#include "upg.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/image_metrics.h"
#include "graphics_utils/modeling.h"
#include "generation.h"
#include "sdf_node.h"
#include <unistd.h>
#include <functional>
#include <chrono>

namespace upg
{
  using SceneDesc = std::pair<UPGStructure, UPGParametersRaw>;
  UPGReconstructionResult benchmark_single(std::string optimizer_name, int total_iterations,
                                           const SceneDesc &ref_scene, bool fixed_structure)
  {
    srand(time(NULL));
    debug("\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 10000
        } 
    }
    generator {

    }
    optimization {
        start {
        }
        step_0 {
            optimizer_name:s = "memetic"
            iterations:i = 100
            verbose:b = false
        }
        step_1 {
            learning_rate:r = 0.003
            iterations:i = 500
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = false
        check_model_quality:b = false
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);

    settings_blk.get_block("optimization")->get_block("step_0")->set_string("optimizer_name", optimizer_name);
    settings_blk.get_block("optimization")->get_block("step_0")->set_int("iterations", total_iterations);
    settings_blk.get_block("optimization")->get_block("step_0")->set_double("finish_threshold", 1e-6);
    settings_blk.get_block("input")->get_block("synthetic_reference")->set_arr("structure", ref_scene.first.s);
    settings_blk.get_block("input")->get_block("synthetic_reference")->set_arr("params", ref_scene.second.p);
    if (fixed_structure)
      settings_blk.get_block("optimization")->get_block("start")->set_arr("structure", ref_scene.first.s);

    return reconstruct_sdf(settings_blk)[0];
  }

  SceneDesc scene_1_sphere()
  {
    SceneDesc desc;
    desc.first.s = {2,1};
    desc.second.p = {0,0,0,1};
    return desc;
  }

  SceneDesc scene_8_spheres()
  {
    SceneDesc desc;
    desc.first.s = {3,3,3,2,1,2,1,3,2,1,2,1,3,3,2,1,2,1,3,2,1,2,1};
    desc.second.p = {0.6,0.6,0.6,0.5, -0.6,0.6,0.6,0.5, 0.6,-0.6,0.6,0.5, -0.6,-0.6,0.6,0.5, 
                     0.6,0  ,0  ,0.5, -0.6,0  ,0  ,0.5, 0  ,-0.6,0  ,0.5, 0   ,0.6 ,0  ,0.5};
    return desc;
  }

  SceneDesc scene_1_box()
  {
    SceneDesc desc;
    desc.first.s = {2,4};
    desc.second.p = {0.2,-0.1,0,0.5,0.5,0.5};
    return desc;
  }

  SceneDesc scene_bubbles(int cnt_x, int cnt_z)
  {
    glm::vec3 p0(-1,1.2, 0.5);
    glm::vec3 p1(1, 1, -0.5);
    
    float base_r = std::min(abs(p1.x-p0.x)/(cnt_x+1), abs(p1.z-p0.z)/(cnt_z+1));

    std::vector<uint16_t> structure_inv;
    std::vector<float> params_inv;
    uint32_t num = 1;
    for (int i=0;i<cnt_z;i++)
    {
      for (int j=0;j<cnt_x;j++)
      {
        glm::vec3 p = glm::vec3(p0.x + (j+1+urand(-0.5,0.5))*(p1.x-p0.x)/(cnt_x+1), 
                                p1.y,
                                p0.z + (i+1+urand(-0.5,0.5))*(p1.z-p0.z)/(cnt_z+1));
        float rnd = urand(0.5,1);
        float r = rnd*base_r;
        structure_inv.push_back(1);
        structure_inv.push_back(2);
        params_inv.push_back(r);
        params_inv.push_back(p.x);
        params_inv.push_back(p.y);
        params_inv.push_back(p.z);

        uint32_t S = 1;
        while (num>= S && ((num & S) == 0))
        {
          structure_inv.push_back(3);
          S = S << 1;
        }
        num++;
      }
    }

    structure_inv.push_back(4);
    structure_inv.push_back(2);
    structure_inv.push_back(3);
    params_inv.push_back(0.5*abs(p0.x-p1.x));
    params_inv.push_back(0.5*abs(p0.y-p1.y));
    params_inv.push_back(0.5*abs(p0.z-p1.z));
    params_inv.push_back(0.5*(p0.x+p1.x));
    params_inv.push_back(0.5*(p0.y+p1.y));
    params_inv.push_back(0.5*(p0.z+p1.z));

    std::vector<uint16_t> structure = structure_inv;
    std::vector<float> params = params_inv;
    for (int i=0;i<structure.size();i++)
      structure[i] = structure_inv[structure.size()-i-1];
    for (int i=0;i<params.size();i++)
      params[i] = params_inv[params.size()-i-1];
    
    SceneDesc desc;
    desc.first.s = structure;
    desc.second.p = params;
    return desc;
  }

  std::vector<UPGReconstructionResult> benchmark_for_optimizer(std::string optimizer_name, bool fixed_structure)
  {
    std::map<std::string, SceneDesc> scenes;
    scenes["1 Sphere"] = scene_1_sphere();
    scenes["8 Spheres"] = scene_8_spheres();
    scenes["1 Box"] = scene_1_box();
    scenes["8 Bubbles"] = scene_bubbles(4, 2);
    scenes["32 Bubbles"] = scene_bubbles(8, 4);

    int max_iters = 200'000;
    std::vector<UPGReconstructionResult> results;
    std::vector<int> timings;
    std::chrono::steady_clock::time_point t1, t2;
    
    for (auto &scene : scenes)
    {
      t1 = std::chrono::steady_clock::now();
      results.push_back(benchmark_single(optimizer_name, max_iters, scene.second, fixed_structure));
      t2 = std::chrono::steady_clock::now();
      timings.push_back(std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count());
    }


    debug("Performed benchmark for optimizer %s\n",optimizer_name.c_str());
    int j = 0;
    for (auto &scene : scenes)
    {
      debug("Scene %d: %s\n", j, scene.first.c_str());
      debug("Value: %.1f*1e-6\n", 1e6*results[j].loss_optimizer);
      debug("Time %d m %d s\n", timings[j]/60, timings[j] % 60);
      j++;
    }

    return results;
  }

  void benchmark_sdf_complex_optimization()
  {
    benchmark_for_optimizer("particle_swarm", true);
  }

  void perform_benchmarks(const Block &blk)
  {
    benchmark_sdf_complex_optimization();
  }
}