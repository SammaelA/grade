#include "upg.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/image_metrics.h"
#include "graphics_utils/modeling.h"
#include "sdf_rendering.h"
#include "generation.h"
#include "sdf_node.h"
#include <unistd.h>
#include <functional>
#include <chrono>
#include "neural_SDF/neural_sdf.h"
#include "neuralCore/neural_network.h"
#include "neuralCore/dataset.h"

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
            points_count:i = 50000
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
            constructive_reconstruction:b = true
            population_size:i = 1000
            tournament_size:i = 50
            local_opt_count:i = 10
            local_opt_iters:i = 500
            local_learning_rate:r = 0.005
            good_soulution_thr:r = 1.0
            part_based_memetic:b = false
            finish_threshold:r = -1.0
        }
        //step_1 {
        //    learning_rate:r = 0.003
        //    iterations:i = 500
        //    verbose:b = false
        //}
    }
    results {
        check_image_quality:b = false
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);

    settings_blk.get_block("optimization")->get_block("step_0")->set_string("optimizer_name", optimizer_name);
    settings_blk.get_block("optimization")->get_block("step_0")->set_int("iterations", total_iterations);
    settings_blk.get_block("input")->get_block("synthetic_reference")->set_arr("structure", ref_scene.first.s);
    settings_blk.get_block("input")->get_block("synthetic_reference")->set_arr("params", ref_scene.second.p);
    if (fixed_structure)
    {
      settings_blk.get_block("optimization")->get_block("start")->set_arr("structure", ref_scene.first.s);
      settings_blk.get_block("optimization")->get_block("step_0")->set_arr("structure", ref_scene.first.s);
      settings_blk.get_block("optimization")->get_block("step_0")->set_arr("params", ref_scene.second.p);
    }

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

  SceneDesc scene_8_rboxes()
  {
    SceneDesc desc;
    desc.first.s = {3,3,3,2,6,2,6,3,2,6,2,6,3,3,2,6,2,6,3,2,6,2,6};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3,0.1, -0.6,0.6,0.6,0.3,0.3,0.3,0.1, 0.6,-0.6,0.6,0.3,0.3,0.3,0.1, -0.6,-0.6,0.6,0.3,0.3,0.3,0.1,
                     0.6,0  ,0  ,0.3,0.3,0.3,0.1, -0.6,0  ,0  ,0.3,0.3,0.3,0.1, 0  ,-0.6,0  ,0.3,0.3,0.3,0.1, 0   ,0.6 ,0  ,0.3,0.3,0.3,0.1};
    return desc;
  }

  SceneDesc scene_1_box()
  {
    SceneDesc desc;
    desc.first.s = {2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3};
    return desc;
  }

  SceneDesc scene_2_boxes()
  {
    SceneDesc desc;
    desc.first.s = {3,2,4,2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3, -0.6,0.6,0.6,0.3,0.3,0.3};
    return desc;
  }

  SceneDesc scene_4_boxes()
  {
    SceneDesc desc;
    desc.first.s = {3,3,2,4,2,4,3,2,4,2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3, -0.6,0.6,0.6,0.3,0.3,0.3, 0.6,-0.6,0.6,0.3,0.3,0.3, -0.6,-0.6,0.6,0.3,0.3,0.3};
    return desc;
  }

  SceneDesc scene_8_boxes()
  {
    SceneDesc desc;
    desc.first.s = {3,3,3,2,4,2,4,3,2,4,2,4,3,3,2,4,2,4,3,2,4,2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3, -0.6,0.6,0.6,0.3,0.3,0.3, 0.6,-0.6,0.6,0.3,0.3,0.3, -0.6,-0.6,0.6,0.3,0.3,0.3,
                     0.6,0  ,0  ,0.3,0.3,0.3, -0.6,0  ,0  ,0.3,0.3,0.3, 0  ,-0.6,0  ,0.3,0.3,0.3, 0   ,0.6 ,0  ,0.3,0.3,0.3};
    return desc;
  }

  SceneDesc scene_bubbles(int cnt_x, int cnt_z)
  {
    glm::vec3 p0(-2,1.2, 0.6);
    glm::vec3 p1(1, 1, -0.6);
    
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

  SceneDesc scene_1_cone()
  {
    SceneDesc desc;
    desc.first.s = {2,8};
    desc.second.p = {0.2,-0.1,0,0.2,sqrtf(1-0.2*0.2),1};
    return desc;
  }

  SceneDesc scene_chair()
  {
    SceneDesc desc;
    desc.first.s = {3, 3,2,4,2,4, 3,3,2,5,2,5,3,2,5,2,5};
    desc.second.p = {0,0.0,0, 0.5,0.05,0.5, 0,0.5,-0.45, 0.5,0.5,0.05,
                     0.4,-0.3,0.4, 0.3,0.05,  -0.4,-0.3,0.4, 0.3,0.05,  0.4,-0.3,-0.4, 0.3,0.05,  -0.4,-0.3,-0.4, 0.3,0.05};
    return desc;    
  }

  std::vector<UPGReconstructionResult> benchmark_for_optimizer(std::string optimizer_name, bool fixed_structure)
  {
    srand(time(NULL));
    std::map<std::string, SceneDesc> scenes;
    //scenes["1 Sphere"] = scene_1_sphere();
    //scenes["8 Spheres"] = scene_8_spheres();
    //scenes["1 Box"] = scene_1_box();
    //scenes["8 Bubbles"] = scene_bubbles(4, 2);
    //scenes["32 Bubbles"] = scene_bubbles(8, 4);
    //scenes["8 Rounded Boxes"] = scene_8_rboxes();
    //scenes["1 Box"] = scene_1_box();
    //scenes["2 Boxes"] = scene_2_boxes();
    //scenes["4 Boxes"] = scene_4_boxes();
    //scenes["8 Boxes"] = scene_8_boxes();
    scenes["Chair"] = scene_chair();

    int max_iters = 200'000;
    int tries = 1;
    std::vector<std::vector<UPGReconstructionResult>> results;
    std::vector<int> timings;
    std::chrono::steady_clock::time_point t1, t2;
    
    for (auto &scene : scenes)
    {
      t1 = std::chrono::steady_clock::now();
      results.emplace_back();
      for (int i=0;i<tries;i++)
      results.back().push_back(benchmark_single(optimizer_name, max_iters, scene.second, fixed_structure));
      t2 = std::chrono::steady_clock::now();
      timings.push_back(std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count());
    }


    debug("Performed benchmark for optimizer %s\n",optimizer_name.c_str());
    int j = 0;
    for (auto &scene : scenes)
    {
      float average_loss = 0;
      float average_quality = 0;
      for (auto &r : results[j])
      {
        average_loss += r.loss_optimizer;
        average_quality += r.quality_ir;
      }
      average_loss /= tries;
      average_quality /= tries;
      timings[j] /= tries;
      debug("Scene %d: %s\n", j, scene.first.c_str());
      debug("Average loss: %.1f*1e-6\n", 1e6*average_loss);
      debug("Average quality: %.1f*1e-6\n", 1e6*average_quality);
      debug("Time %d m %d s\n", timings[j]/60, timings[j] % 60);
      j++;
    }

    return results[0];
  }

  void benchmark_sdf_complex_optimization()
  {
    //benchmark_for_optimizer("CC", true);
    //benchmark_for_optimizer("DE", true);
    benchmark_for_optimizer("memetic", true);
    //benchmark_for_optimizer("particle_swarm", true);
    //benchmark_for_optimizer("iterative_fitting", true);
  }

  void neural_sdf_test()
  {
    /*
    nn::TensorProcessor::init("GPU");
    unsigned count = 25000;

    SceneDesc scene = scene_4_boxes();
    SdfGenInstance gen(scene.first);
    ProceduralSdf sdf = gen.generate(scene.second.p);
    std::vector<glm::vec3> positions;
    std::vector<float> distances;
    sdf_to_point_cloud_with_dist(sdf, count, &positions, &distances);
    AABB bbox = sdf.root->get_bbox();
    std::vector<float> positions_f;
    positions_f.reserve(positions.size()*3);
    for (auto &p : positions)
    {
      positions_f.push_back(p.x);
      positions_f.push_back(p.y);
      positions_f.push_back(p.z);
    }

    nn::Siren network(nn::Siren::Type::SDF, 2, 32);
    network.train(positions_f, distances, 512, 5000);
  
    CameraSettings camera;
    camera.origin = glm::vec3(0,0,3);
    camera.target = glm::vec3(0,0,0);
    camera.up = glm::vec3(0,1,0);

    Texture t = render_neural_sdf(network, bbox, camera, 256, 256, 9, true);
    engine::textureManager->save_png(t, "NN demo");
    */
  }

  void voxNet_create_dataset(int size)
  {
    unsigned samples = 16;
    unsigned vox_size = 32;
    unsigned models_count = size;
    std::vector<float> voxels(vox_size*vox_size*vox_size*models_count, 0.0f);
    std::vector<float> labels(2*models_count, 0.0f);

    CameraSettings camera;
    camera.origin = glm::vec3(0,0,3);
    camera.target = glm::vec3(0,0,0);
    camera.up = glm::vec3(0,1,0);

    for (int mod=0;mod<models_count;mod++)
    {
      UPGStructure s;
      UPGParametersRaw p;
      if (urand() < 0.5)
      {
        s.s = {SdfNode::OR, SdfNode::MOVE, SdfNode::SPHERE, SdfNode::MOVE, SdfNode::BOX};
        p.p = {(float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),  (float)urand(0.3, 0.6),
               (float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),  (float)urand(0.3, 0.6),(float)urand(0.3, 0.6),(float)urand(0.3, 0.6)};
        labels[2*mod + 0] = 1;
        labels[2*mod + 1] = 0;
      }
      else
      {
        s.s = {SdfNode::MOVE, SdfNode::BOX, SdfNode::MOVE, SdfNode::BOX};
        p.p = {(float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),  (float)urand(0.3, 0.6),(float)urand(0.3, 0.6),(float)urand(0.3, 0.6),
               (float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),(float)urand(-0.4, 0.4),  (float)urand(0.3, 0.6),(float)urand(0.3, 0.6),(float)urand(0.3, 0.6)};
        labels[2*mod + 0] = 0;
        labels[2*mod + 1] = 1;
      }
      SdfGenInstance gen(s);
      ProceduralSdf sdf = gen.generate(p.p);

      //Texture t = render_sdf(sdf, camera, 256, 256, 4, true);
      //engine::textureManager->save_png(t, "dataset_image_"+std::to_string(mod));

      for (int i=0;i<vox_size;i++)
      {
        for (int j=0;j<vox_size;j++)
        {
          for (int k=0;k<vox_size;k++)
          {
            voxels[mod*vox_size*vox_size*vox_size + i*vox_size*vox_size + j*vox_size + k] = 0;
            for (int s=0;s<samples;s++)
            {
              glm::vec3 p = {(k+urand())/vox_size, (j+urand())/vox_size, (i+urand())/vox_size};
              p = 2.0f*p - 1.0f;
              voxels[mod*vox_size*vox_size*vox_size + i*vox_size*vox_size + j*vox_size + k] += sdf.get_distance(p) < 0;
            }
            voxels[mod*vox_size*vox_size*vox_size + i*vox_size*vox_size + j*vox_size + k] /= samples;
            //debug("%.2f ", voxels[mod*vox_size*vox_size*vox_size + i*vox_size*vox_size + j*vox_size + k]);
          }
          //debugnl();
        }
        //debugnl();
        //debugnl();
      }
      //debug("####\n");
      if (mod % 100 == 0)
        debug("generating %d/%d\n", mod, models_count);
    }

    nn::Dataset dataset;
    dataset.train_data = voxels;
    dataset.train_labels = labels;
    dataset.train_elements = models_count;
    dataset.label_size = 2;
    dataset.element_size = {vox_size, vox_size, vox_size};

    nn::save_dataset("saves/voxNet2D_test_dataset.bin", &dataset);
  }

  void voxNet2D_test()
  {
    
    nn::Dataset dataset;
    nn::load_dataset("saves/voxNet2D_test_dataset.bin", &dataset);

    nn::TensorProcessor::init("GPU");
    nn::NeuralNetwork net;
    net.add_layer(std::make_shared<nn::Conv2DLayer>(32,32,32, 64, 5), nn::Initializer::GlorotNormal);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    net.add_layer(std::make_shared<nn::MaxPoolingLayer>(28, 28, 64));
    net.add_layer(std::make_shared<nn::Conv2DLayer>(14,14,64, 64), nn::Initializer::GlorotNormal);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    net.add_layer(std::make_shared<nn::MaxPoolingLayer>(12, 12, 64));
    net.add_layer(std::make_shared<nn::FlattenLayer>(6,6,64));
    net.add_layer(std::make_shared<nn::DenseLayer>(6*6*64, 512), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    net.add_layer(std::make_shared<nn::DenseLayer>(512, 512), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    net.add_layer(std::make_shared<nn::DenseLayer>(512, 2), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::SoftMaxLayer>());
    net.train(dataset.train_data.data(), dataset.train_labels.data(), dataset.train_elements, 64, 500, true, nn::OptimizerAdam(0.001f), 
              nn::Loss::CrossEntropy, nn::Metric::Accuracy, true);
  }

  void benchmark_sdf_rendering(int image_size, int spp)
  {
    std::map<std::string, SceneDesc> scenes;
    scenes["1 Cone"] = scene_1_cone();
    scenes["1 Sphere"] = scene_1_sphere();
    scenes["8 Spheres"] = scene_8_spheres();
    scenes["1 Box"] = scene_1_box();
    scenes["8 Bubbles"] = scene_bubbles(4, 2);
    scenes["32 Bubbles"] = scene_bubbles(8, 4);
    scenes["8 Rounded Boxes"] = scene_8_rboxes();
    scenes["Chair"] = scene_chair();

    CameraSettings camera;
    camera.origin = glm::vec3(0,0,3);
    camera.target = glm::vec3(0,0,0);
    camera.up = glm::vec3(0,1,0);

    std::chrono::steady_clock::time_point t1, t2;

    for (auto &scene : scenes)
    {
      SdfGenInstance gen(scene.second.first);
      ProceduralSdf sdf = gen.generate(scene.second.second.p);
      t1 = std::chrono::steady_clock::now();
      Texture t = render_sdf(sdf, camera, image_size, image_size, spp, true);
      t2 = std::chrono::steady_clock::now();
      float time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
      debug("%s rendered in %.1f s. %d kRays/s\n", scene.first.c_str(), 1e-3*time_ms, (int)((image_size*image_size*spp)/time_ms));
      engine::textureManager->save_png(t, scene.first+" demo");
    }

  }

  void perform_benchmarks(const Block &blk)
  {
    std::string name = blk.get_string("name", "rendering");
    if (name == "rendering")
      benchmark_sdf_rendering(512, 1);
    else if (name == "nsdf")
      nsdf::neural_SDF_test();
    else if (name == "vox_net")
    {
      //voxNet_create_dataset(5000);
      voxNet2D_test();
      return;      
    }
    else
      benchmark_sdf_complex_optimization();
  }
}