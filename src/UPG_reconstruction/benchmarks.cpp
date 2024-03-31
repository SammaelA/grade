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
#include "sdf_grid.h"
#include "optimization.h"
#include "density_field_process.h"
#include "sdf_scene_convert.h"
#include "sdf_octree.h"
#include "graphics_utils/modeling.h"
#include "graphics_utils/render_wireframe.h"
#include "LiteRT/IRenderer.h"
#include "LiteMath/Image2d.h"

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
    desc.second.p = {0,0,0,0.8};
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
    float3 p0(-2,1.2, 0.6);
    float3 p1(1, 1, -0.6);
    
    float base_r = std::min(abs(p1.x-p0.x)/(cnt_x+1), abs(p1.z-p0.z)/(cnt_z+1));

    std::vector<uint16_t> structure_inv;
    std::vector<float> params_inv;
    uint32_t num = 1;
    for (int i=0;i<cnt_z;i++)
    {
      for (int j=0;j<cnt_x;j++)
      {
        float3 p = float3(p0.x + (j+1+urand(-0.5,0.5))*(p1.x-p0.x)/(cnt_x+1), 
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
    desc.second.p = {0,0.0,0, 0.5,0.07,0.5, 0,0.45,-0.45, 0.5,0.45,0.07,
                     0.4,-0.3,0.4, 0.3,0.07,  -0.4,-0.3,0.4, 0.3,0.07,  0.4,-0.3,-0.4, 0.3,0.07,  -0.4,-0.3,-0.4, 0.3,0.07};
    return desc;    
  }

  SceneDesc scene_complex_chair()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::CROTATE, SdfNodeType::CHAIR};
    desc.second.p = {1, 0, 1, 0.15, 0.15, 0.15, 0.05, 0.3, 0.5, 0.05, 0.5, 0.5};
    return desc;    
  }

  SceneDesc scene_extrusion()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::ROTATE, SdfNodeType::EXTRUSION, SdfNodeType::CIRCLE};
    desc.second.p = {0.0, 0.0, 1.0, 0.7, 0.5};
    return desc;    
  }

  SceneDesc scene_subtraction()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::SUBTRACT, SdfNodeType::MOVE, SdfNodeType::SPHERE, SdfNodeType::MOVE, SdfNodeType::SPHERE};
    desc.second.p = {0,0,0,0.6, 0.5,0.0,0.5,0.4};
    return desc;    
  }

  SceneDesc scene_CSG_1()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::SCALE, SdfNodeType::ROTATE,
                    SdfNodeType::SUBTRACT, 
                    SdfNodeType::AND, SdfNodeType::OR, SdfNodeType::MOVE, SdfNodeType::BOX, SdfNodeType::MOVE, SdfNodeType::BOX,
                                      SdfNodeType::MOVE, SdfNodeType::SPHERE,
                    SdfNodeType::AND, SdfNodeType::MOVE, SdfNodeType::CYLINDER,
                                      SdfNodeType::MOVE, SdfNodeType::BOX};
    desc.second.p = {0.5, PI/2,0,-PI/6,
                    -1,1.5  ,0, 1,1.5,0.5,
                     1.5,1  ,0, 1.5,1,0.5,
                     0  ,0  ,0, 3,
                     0,0.5,0.5, 2,0.5,
                    -1  ,1  ,0, 0.75,1,3};
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
    ProceduralSdf sdf(scene.first);
    sdf.set_parameters(scene.second.p);
    std::vector<float3> positions;
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
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    Texture t = render_neural_sdf(network, bbox, camera, 256, 256, 9, true);
    engine::textureManager->save_png(t, "NN demo");
    */
  }

  void add_circle(UPGStructure &s, UPGParametersRaw &p)
  {
    s.s.push_back(SdfNodeType::SPHERE);
    p.p.push_back((float)urand(0.3, 0.6));
  }

  void add_box(UPGStructure &s, UPGParametersRaw &p)
  {
    s.s.push_back(SdfNodeType::BOX);
    p.p.push_back((float)urand(0.3, 0.6));
    p.p.push_back((float)urand(0.3, 0.6));
    p.p.push_back((float)urand(0.3, 0.6));    
  }

  void add_cylinder(UPGStructure &s, UPGParametersRaw &p)
  {
    s.s.push_back(SdfNodeType::CYLINDER);
    p.p.push_back((float)urand(0.3, 0.6));
    p.p.push_back((float)urand(0.3, 0.6));
    p.p.push_back((float)urand(0.3, 0.6));    
  }

  void add_r_box(UPGStructure &s, UPGParametersRaw &p)
  {
    s.s.push_back(SdfNodeType::ROUNDED_BOX);
    p.p.push_back((float)urand(0.3, 0.6));
    p.p.push_back((float)urand(0.3, 0.6));
    p.p.push_back((float)urand(0.3, 0.6));  
    p.p.push_back((float)urand(0.01, 0.3));   
  }

  void voxNet_create_dataset(int size)
  {
    unsigned samples = 8;
    unsigned vox_size = 32;
    unsigned models_count = size;
    std::vector<float> voxels(vox_size*vox_size*vox_size*models_count, 0.0f);
    std::vector<float> labels(2*models_count, 0.0f);

    CameraSettings camera;
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    #pragma omp parallel for
    for (int mod=0;mod<models_count;mod++)
    {
      UPGStructure s;
      UPGParametersRaw p;

      bool have_sphere = false;
      int cnt = urand(1.1, 9.9);
      for (int i=0; i<cnt; i++)
      {
        if (i != cnt-1)
          s.s.push_back(SdfNodeType::OR);
        
        s.s.push_back(SdfNodeType::MOVE);
        p.p.push_back((float)urand(-0.5, 0.5));
        p.p.push_back((float)urand(-0.5, 0.5));
        p.p.push_back((float)urand(-0.5, 0.5));

        s.s.push_back(SdfNodeType::ROTATE);
        p.p.push_back((float)urand(-PI, PI));
        p.p.push_back((float)urand(-PI, PI));
        p.p.push_back((float)urand(-PI, PI));

        float r = urand();
        if (r < 0.15)
        {
          have_sphere = true;
          add_circle(s, p);
        }
        else if (r < 0.3)
          add_cylinder(s, p);
        else if (r < 0.45)
          add_r_box(s, p);
        else 
          add_box(s, p);
      }

      if (have_sphere)
      {
        labels[2*mod + 0] = 1;
        labels[2*mod + 1] = 0;
      }
      else
      {
        labels[2*mod + 0] = 0;
        labels[2*mod + 1] = 1;
      }

      ProceduralSdf sdf(s);
      sdf.set_parameters(p.p);

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
              float3 p = {(k+urand())/vox_size, (j+urand())/vox_size, (i+urand())/vox_size};
              p = 2.0f*p - 1.0f;
              float d = sdf.get_distance(p);
              if (d > 2.0/vox_size)
              {
                voxels[mod*vox_size*vox_size*vox_size + i*vox_size*vox_size + j*vox_size + k] = 0;
                break;
              }
              else if (d < -2.0/vox_size)
              {
                voxels[mod*vox_size*vox_size*vox_size + i*vox_size*vox_size + j*vox_size + k] = samples;
                break;
              } 
              else
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

    nn::save_dataset("saves/voxNet2D_diverse_dataset.bin", &dataset);
  }

  void voxNet2D_test()
  {
    
    nn::Dataset dataset;
    nn::load_dataset("saves/voxNet2D_diverse_dataset.bin", &dataset);
    for (int i=0;i<dataset.train_elements;i++)
      logerr("label %d %f %f", i, dataset.train_labels[2*i+0], dataset.train_labels[2*i+1]);

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

  void CNN_2_5D_create_dataset(int size)
  {
    unsigned samples = 8;
    unsigned image_size = 64;
    unsigned layers = 8;
    unsigned models_count = size;
    std::vector<float> images(image_size*image_size*layers*models_count, 0.0f);
    std::vector<float> labels(2*models_count, 0.0f);

    std::vector<CameraSettings> cameras(layers);
    for (int i=0;i<layers;i++)
    {
      float phi = urand(0, 2*PI);
      float psi = urand(-PI/2, PI/2);
      float d = 4;
      cameras[i].origin = float3(d*cos(phi)*cos(psi), d*sin(psi), d*sin(phi)*cos(psi));
      cameras[i].target = float3(0,0,0);
      cameras[i].up = float3(-d*cos(phi)*sin(psi), d*cos(psi), -d*sin(phi)*sin(psi));
      cameras[i].z_near = 1;
      cameras[i].z_far = 7;
      cameras[i].orthographic = true;
    }

    #pragma omp parallel for
    for (int mod=0;mod<models_count;mod++)
    {
      UPGStructure s;
      UPGParametersRaw p;

      bool have_sphere = false;
      int cnt = urand(1.1, 5.9);
      for (int i=0; i<cnt; i++)
      {
        if (i != cnt-1)
          s.s.push_back(SdfNodeType::OR);
        
        s.s.push_back(SdfNodeType::MOVE);
        p.p.push_back((float)urand(-0.5, 0.5));
        p.p.push_back((float)urand(-0.5, 0.5));
        p.p.push_back((float)urand(-0.5, 0.5));

        s.s.push_back(SdfNodeType::ROTATE);
        p.p.push_back((float)urand(-PI, PI));
        p.p.push_back((float)urand(-PI, PI));
        p.p.push_back((float)urand(-PI, PI));

        float r = urand();
        if (r < 0.2)
        {
          have_sphere = true;
          add_circle(s, p);
        }
        else if (r < 0.3)
          add_cylinder(s, p);
        else if (r < 0.45)
          add_r_box(s, p);
        else 
          add_box(s, p);
      }

      if (have_sphere)
      {
        labels[2*mod + 0] = 1;
        labels[2*mod + 1] = 0;
      }
      else
      {
        labels[2*mod + 0] = 0;
        labels[2*mod + 1] = 1;
      }

      ProceduralSdf sdf(s);
      sdf.set_parameters(p.p);

      for (int i=0;i<layers;i++)
      {
        std::span<float> im = std::span<float>(images.data()+image_size*image_size*layers*mod + image_size*image_size*i, image_size*image_size);
        render_sdf_to_array(im, sdf, cameras[i], image_size, image_size, 4, SDFRenderMode::INVERSE_LINEAR_DEPTH);
      
        if (mod < 0)
        {
          Texture t = render_sdf(sdf, cameras[i], image_size, image_size, 9, SDFRenderMode::INVERSE_LINEAR_DEPTH);
          engine::textureManager->save_png(t, "2.5D_dataset_image_"+std::to_string(mod)+"_view_"+std::to_string(i));
        }
      }
      if (mod % 100 == 0)
        debug("generating %d/%d\n", mod, models_count);
    }

    nn::Dataset dataset;
    dataset.train_data = images;
    dataset.train_labels = labels;
    dataset.train_elements = models_count;
    dataset.label_size = 2;
    dataset.element_size = {image_size, image_size, layers};

    nn::save_dataset("saves/CNN_2.5D_dataset_25000.bin", &dataset);
  }
  
  void CNN_2_5D_get_net(nn::NeuralNetwork &net)
  {
    net.add_layer(std::make_shared<nn::Conv2DLayer>(64,64,8,8), nn::Initializer::GlorotNormal);
    net.add_layer(std::make_shared<nn::MaxPoolingLayer>(62, 62, 8));
    net.add_layer(std::make_shared<nn::FlattenLayer>(31,31,8));

    net.add_layer(std::make_shared<nn::DenseLayer>(31*31*8, 512), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());

    net.add_layer(std::make_shared<nn::DenseLayer>(512, 256), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());

    net.add_layer(std::make_shared<nn::DenseLayer>(256, 128), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());

    net.add_layer(std::make_shared<nn::DenseLayer>(128, 2), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::SoftMaxLayer>());
  }

  void CNN_2_5D_train()
  {
    nn::Dataset dataset;
    nn::load_dataset("saves/CNN_2.5D_dataset_25000.bin", &dataset);
    nn::rebalance_classes(&dataset);

    nn::TensorProcessor::init("GPU");
    nn::NeuralNetwork net;
    CNN_2_5D_get_net(net);

    net.train(dataset.train_data.data(), dataset.train_labels.data(), dataset.train_elements, 100, 50, true,
              nn::OptimizerRMSProp(3e-4f, 0.999, 1e-6, true), 
              nn::Loss::CrossEntropy, nn::Metric::Accuracy, true);
    net.save_weights_to_file("./saves/CNN_2.5D_4p_weights.bin");
  }

  void CNN_2_5D_test()
  {
    nn::Dataset dataset;
    nn::load_dataset("saves/CNN_2.5D_dataset_10000.bin", &dataset);

    nn::TensorProcessor::init("GPU");
    nn::NeuralNetwork net;
    CNN_2_5D_get_net(net);
    net.initialize_from_file("./saves/CNN_2.5D_4p_weights.bin");

    std::vector<float> out_labels(dataset.train_labels.size(), 0.0f);

    net.evaluate(dataset.train_data.data(), out_labels.data(), dataset.train_elements);
    float acc = net.calculate_metric(out_labels.data(), dataset.train_labels.data(), dataset.train_elements, nn::Metric::Accuracy);
    logerr("CNN_2_5D_train Accuracy: %f (%u elements)\n", acc, dataset.train_elements);
  }

  void sdf_grid_test()
  {
    SceneDesc scene = scene_chair();
    ProceduralSdf sdf(scene.first);
    sdf.set_parameters(scene.second.p);

    unsigned vox_size = 32;
    unsigned samples = 4;
    std::vector<float> data(vox_size*vox_size*vox_size, 0.0f);
    AABB bbox = AABB({-1,-1,-1},{1,1,1});
    GridSdfNode grid(SdfNodeType::GRID_32, vox_size);
    grid.set_param_span(data, 0);

    bool direct = false;
    if (direct)
    {
      for (int i=0;i<vox_size;i++)
      {
        for (int j=0;j<vox_size;j++)
        {
          for (int k=0;k<vox_size;k++)
          {
            float d = 0.0;
            for (int s=0;s<samples;s++)
            {
              float3 p = {(k+0.5)/vox_size, (j+0.5)/vox_size, (i+0.5)/vox_size};
              p = bbox.size()*p + bbox.min_pos;
              d += sdf.get_distance(p);
            }
            grid.set_voxel(uint3(k,j,i), d/samples);
          }
        }
      }
    }
    else
    {
      unsigned params_cnt = vox_size*vox_size*vox_size;
      std::vector<float> ddist_dpos(3, 0.0f);
      std::vector<float> ddist_dp(params_cnt, 0.0f);
      std::vector<float> x_grad(params_cnt, 0.0f);
      std::vector<float> V(params_cnt, 0.0f);
      std::vector<float> S(params_cnt, 0.0f);
      std::vector<float> &X = data;
      float alpha = 0.01;
      float beta_1 = 0.9;
      float beta_2 = 0.999;
      float eps = 1e-8;
      bool verbose = true;

      unsigned samples = params_cnt;
      unsigned iterations = 300;
      
      for (int iter=0; iter< iterations; iter++)
      {
        std::fill_n(x_grad.begin(), x_grad.size(), 0.0f);
        double total_loss = 0.0;
        for (int i=0;i<samples;i++)
        {
          int x = i % vox_size;
          int y = i / vox_size % vox_size;
          int z = i / (vox_size*vox_size);
          float3 p = {(x+urand())/vox_size, (y+urand())/vox_size, (z+urand())/vox_size};
          p = bbox.size()*p + bbox.min_pos;

          std::fill_n(ddist_dp.begin(), ddist_dp.size(), 0.0);
          float d = grid.get_distance(p, ddist_dp, ddist_dpos) - sdf.get_distance(p);
          total_loss += d*d;
          for (int j=0;j<ddist_dp.size();j++)
            x_grad[j] += 2*d*ddist_dp[j] / samples;
        }
        float val = total_loss/samples;
        //debug("grad  [");
        for (int i=0;i<params_cnt;i++)
        {
          //debug("%f ",x_grad[i]);
          float g = x_grad[i];
          V[i] = beta_1 * V[i] + (1-beta_1)*g;
          float Vh = V[i] / (1 - pow(beta_1, iter+1)); 
          S[i] = beta_2 * S[i] + (1-beta_2)*g*g;
          float Sh = S[i] / (1 - pow(beta_2, iter+1)); 
          X[i] -= alpha*Vh/(sqrt(Sh) + eps);
        }
        //debug("]\n");
        if ((iter % 5 == 0) && verbose)
          debug("Adam iter %3d  val = %.8f\n", iter, val);
      }
    }

    int steps = 15;
    ProceduralSdf g_sdf({{SdfNodeType::GRID_32}});
    g_sdf.set_parameters(data);
    for (int i=0;i<steps;i++)
    {
      CameraSettings camera;
      camera.origin = float3(3*cos((2.0f*PI*i)/steps),0,3*sin((2.0f*PI*i)/steps));
      camera.target = float3(0,0,0);
      camera.up = float3(0,1,0);
      Texture t = render_sdf(g_sdf, camera, 512, 512, 16, SDFRenderMode::LAMBERT);
      engine::textureManager->save_png(t, "dataset_image_grid_"+std::to_string(i));
    }
  }

  void primitive_SDF_to_grid(ProceduralSdf &sdf, const AABB &bbox, float *grid, unsigned vox_size,
                             unsigned samples)
  {
    int offset = 0;
    for (int i = 0; i < vox_size; i++)
    {
      for (int j = 0; j < vox_size; j++)
      {
        for (int k = 0; k < vox_size; k++)
        {
          float3 p = {(k + 0.5) / vox_size, (j + 0.5) / vox_size, (i + 0.5) / vox_size};
          p = bbox.size() * p + bbox.min_pos;
          float d = sdf.get_distance(p);
          grid[offset] = d;
          offset++;
        }
      }
    }
  }

  void grid_demonstrate_different_sizes()
  {
    auto scene = scene_chair();
    ProceduralSdf reference_sdf(scene.first);
    reference_sdf.set_parameters(scene.second.p);
    AABB bbox = reference_sdf.get_bbox();
    float max_s = 0.5f*MAX(bbox.size().x, MAX(bbox.size().y, bbox.size().z));
    float3 bbox_center = 0.5f*(bbox.min_pos + bbox.max_pos);
    bbox = AABB(bbox_center - float3(max_s), bbox_center + float3(max_s)).expand(1.1f);
    bbox = AABB({-1,-1,-1},{1,1,1});

    std::vector<uint16_t> grid_types = {SdfNodeType::GRID_16, SdfNodeType::GRID_32, SdfNodeType::GRID_64, SdfNodeType::GRID_128};
    std::vector<uint16_t> grid_sizes = {16,32,64,128,256};
    std::vector<uint16_t> full_structure = {3, 2, 3, 3,2,4,2,4, 3,3,2,5,2,5,3,2,5,2,5};
    std::vector<float> full_params = {3,0,0,
                                      0,0.0,0, 0.5,0.07,0.5, 0,0.45,-0.45, 0.5,0.45,0.07,
                                      0.4,-0.3,0.4, 0.3,0.07,  -0.4,-0.3,0.4, 0.3,0.07,  
                                      0.4,-0.3,-0.4, 0.3,0.07,  -0.4,-0.3,-0.4, 0.3,0.07};
    std::vector<float3> moves = {float3(-3,0,0), float3(-1.5,0,0), float3(0,0,0), float3(1.5,0,0)};

    for (int i=0;i<grid_types.size();i++)
    {
      UPGStructure structure = {{SdfNodeType::MOVE, SdfNodeType::SCALE, grid_types[i]}};
      std::vector<float> parameters(4 + grid_sizes[i]*grid_sizes[i]*grid_sizes[i], 0.0f);

      primitive_SDF_to_grid(reference_sdf, bbox, parameters.data()+4, grid_sizes[i], 1);

      parameters[0] = moves[i].x;
      parameters[1] = moves[i].y;
      parameters[2] = moves[i].z;
      parameters[3] = 1;
      if (i!=grid_types.size()-1)
        full_structure.push_back(SdfNodeType::OR);
      full_structure.insert(full_structure.end(), structure.s.begin(), structure.s.end());
      full_params.insert(full_params.end(), parameters.begin(), parameters.end());
    }

    ProceduralSdf sdf({full_structure});
    sdf.set_parameters(full_params);
    full_params.clear();

    unsigned steps = 15;
    for (int i=0;i<steps;i++)
    {
      CameraSettings camera;
      camera.origin = float3(7*cos((2.0f*PI*i)/steps),2,7*sin((2.0f*PI*i)/steps));
      camera.target = float3(0,0,0);
      camera.up = float3(0,1,0);
      Texture t = render_sdf(sdf, camera, 2048, 2048, 9, SDFRenderMode::LAMBERT);
      engine::textureManager->save_png(t, "Grid SDFs demo "+std::to_string(i));
    }
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
    scenes["Complex Chair"] = scene_complex_chair();
    scenes["Subtraction"] = scene_subtraction();
    scenes["Extrusion"] = scene_extrusion();

    CameraSettings camera;
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    std::chrono::steady_clock::time_point t1, t2;

    for (auto &scene : scenes)
    {
      ProceduralSdf sdf(scene.second.first);
      sdf.set_parameters(scene.second.second.p);
      t1 = std::chrono::steady_clock::now();
      Texture t = render_sdf(sdf, camera, image_size, image_size, spp, SDFRenderMode::LAMBERT);
      t2 = std::chrono::steady_clock::now();
      float time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
      debug("%s rendered in %.1f s. %d kRays/s\n", scene.first.c_str(), 1e-3*time_ms, (int)((image_size*image_size*spp)/time_ms));
      engine::textureManager->save_png(t, scene.first+" demo");
    }

  }

  /*
    the argument 'cuts_quant' in the 'get_bbox_list' (see below) function shows
    how many time we want to cut each bbox in half along some axis
  */
  std::vector<AABB> get_bbox_list(std::function<float(const float3 &)> sdf,
                                  const AABB &sdf_bbox, int bbox_count, int cuts_quant){
    std::vector<AABB> bbox_list;
    
    int bbox_segment_quant = pow(2, cuts_quant);
    float3 size_bbox = sdf_bbox.size();
    float3 segment_bbox = size_bbox / bbox_segment_quant;

    // Thanks to the point_xyz variable, we move between the "bboxes"
    float3 point_xyz = segment_bbox / 2;
    
    int semidiag_len = sqrt(pow(point_xyz.x, 2) + pow(point_xyz.y, 2) +
                            pow(point_xyz.z, 2));

    for(unsigned int axis_x = 0; axis_x < bbox_segment_quant; axis_x += 1){
      for(unsigned int axis_y = 0; axis_y < bbox_segment_quant; axis_y += 1){
        for(unsigned int axis_z = 0; axis_z < bbox_segment_quant; axis_z += 1){
          if (sdf(point_xyz) >= semidiag_len){
            bbox_list.push_back(AABB(point_xyz - segment_bbox / 2, point_xyz +
                                     segment_bbox / 2));
          }
          point_xyz.z += segment_bbox.z;
        }
        point_xyz.z = 0;
        point_xyz.y += segment_bbox.y;
      }
      point_xyz.y = 0;
      point_xyz.x += segment_bbox.x;
    }

    return bbox_list;
  }
  
  void sdf_scene_test()
  {
    SceneDesc s = scene_CSG_1();
    {
      SdfScene scene = create_sdf_scene(s.first, s.second);
      save_sdf_scene(scene, "saves/test_scene.bin");
      save_sdf_scene_hydra(scene, "../LiteRT/scenes/02_sdf_scenes", "csg_new");
    }
    SdfScene scene;
    load_sdf_scene(scene, "saves/test_scene.bin");

    CameraSettings camera;
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    std::chrono::steady_clock::time_point t1, t2;

    t1 = std::chrono::steady_clock::now();
    Texture t = render_sdf(scene, camera, 512, 512, 16, SDFRenderMode::LAMBERT);
    t2 = std::chrono::steady_clock::now();
    float time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    debug("%s rendered in %.1f s. %d kRays/s\n", "SDF Scene", 1e-3 * time_ms, (int)((512 * 512 * 16) / time_ms));
    engine::textureManager->save_png(t, "SDF Scene demo");
  }

  static std::vector<float4> octree_visualize_palette = {
    float4(0.5,0.5,0.5,1), //0-gray
    float4(1,0,0,1),       //1-red
    float4(0,1,0,1),       //2-green
    float4(0,0,1,1),       //3-blue
    float4(1,1,0,1),       //4-yellow
    float4(0,1,1,1),       //5-cyan
    float4(1,0,1,1)        //6-magenta
  };

  void sdf_octree_to_model_rec(const SparseOctree &octree, Mesh *m, unsigned level_from, unsigned level_to, bool border_only,
                               unsigned idx, unsigned level, float3 p)
  {
    if (level == level_to)
      return;
    if (octree.get_node(idx).offset > 0)
    {
      for (int i=0;i<8;i++)
      {
        float3 dp = float3((i & 4) >> 2, (i & 2) >> 1, i & 1);
        sdf_octree_to_model_rec(octree, m, level_from, level_to, border_only, octree.get_node(idx).offset + i, level+1, 2.0f*p + dp);
      }
    }

    bool last_level = level >= level_from && (octree.get_node(idx).offset == 0 || level == level_to-1);
    
    if (last_level &&
       (octree.get_node(idx).value < 1e-6 || !border_only))
    {
      //printf("add level %d\n", level);
      float d = 1.0f/pow(2.0f, level);
      float3 pos = 2.0f*(p*d) - 1.0f;
      Box b = Box(pos, float3(2*d, 0, 0), float3(0, 2*d, 0), float3(0, 0, 2*d));
      visualizer::body_to_model(&b, m, true, octree_visualize_palette[std::max(level,(unsigned)octree_visualize_palette.size()-1)]);
    }
  }

  Texture visualize_sdf_octree(const SparseOctree &octree, CameraSettings cam, unsigned level_from, unsigned level_to)
  {
    //Box b = Box(shift + float3(scale.x * x, scale.y * y, scale.z * z), float3(dot_size, 0, 0), float3(0, dot_size, 0), float3(0, 0, dot_size));
    Model *m = new Model();
    sdf_octree_to_model_rec(octree, m, level_from, level_to, true,  0, 0,float3(0,0,0));
    m->update();
    WireframeRenderer wr;
    Texture t = wr.render(*m, cam.get_viewProj(false), 4096, 4096);
    delete m;
    return t;
  }

  void sdf_octree_test()
  {
    SceneDesc s = scene_complex_chair();
    ProceduralSdf sdf(s.first);
    sdf.set_parameters(s.second.p);
    std::vector<float> data = {0};

    ProceduralSdf g_sdf({{SdfNodeType::OCTREE}});
    g_sdf.set_parameters(data);
    dynamic_cast<OctreeSdfNode*>(g_sdf.root)->construct( [&sdf](const float3 &p) {
      /*logerr("sample %f %f %f",p.x,p.y,p.z);*/ return sdf.get_distance(p); });
    for (int i=0;i<0;i++)
    {
      float3 p = float3(urand(-1,1), 1, urand(-1,1));
      float d = dynamic_cast<OctreeSdfNode*>(g_sdf.root)->octree.sample(p);
      float d2 = sdf.get_distance(p);
      float d3 = dynamic_cast<OctreeSdfNode*>(g_sdf.root)->octree.sample_mip_skip_3x3(p);
      //if (d*d2<0)
        printf("%f %f %f -- d= %f %f %f\n",p.x,p.y,p.z, d,d2,d3);
    }
    //return;

    SparseOctree &octree = dynamic_cast<OctreeSdfNode*>(g_sdf.root)->octree;
    printf("Octree size %d\n", (int)(octree.get_nodes().size()));
    for (int i=3;i<9;i++)
    {
      CameraSettings camera;
      camera.origin = float3(3,0,3);
      camera.target = float3(0,0,0);
      camera.up = float3(0,1,0);
      Texture t = visualize_sdf_octree(octree, camera, 3, i+1);
      engine::textureManager->save_png(t, "octree_level_"+std::to_string(i));
    }

    //return;
    int steps = 1;
    for (int i=0;i<steps;i++)
    {
      CameraSettings camera;
      //camera.origin = float3(3*cos((2.0f*PI*i)/steps),0,3*sin((2.0f*PI*i)/steps));
      camera.origin = float3(3,0,3);
      camera.target = float3(0,0,0);
      camera.up = float3(0,1,0);
      Texture t = render_sdf(g_sdf, camera, 1024, 1024, 1, SDFRenderMode::LAMBERT);
      engine::textureManager->save_png(t, "octree_test_"+std::to_string(i));
    }
  }

  void liteRT_render_test()
  {
    SceneDesc s = scene_CSG_1();
    SdfScene scene = create_sdf_scene(s.first, s.second);

    CameraSettings camera;
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    unsigned W = 1024, H = 1024;
    LiteImage::Image2D<uint32_t> image(W, H);

    auto pRender = MakeEyeRayShooterRenderer("GPU");
    pRender->SetAccelStruct(CreateSceneRT("BVH2Common", "cbvh_embree2", "SuperTreeletAlignedMerged4"));
    pRender->SetViewport(0,0, W, H);
    pRender->UpdateCamera(camera.get_view(), camera.get_proj(false));

    pRender->GetAccelStruct()->ClearGeom();
    pRender->GetAccelStruct()->AddGeom_SdfScene(scene);
    pRender->GetAccelStruct()->ClearScene();
    pRender->GetAccelStruct()->AddInstance(0, LiteMath::float4x4());
    pRender->GetAccelStruct()->AddInstance(0, LiteMath::translate4x4(float3(-100,-100,-100)));
    pRender->GetAccelStruct()->CommitScene();
    pRender->CommitDeviceData();

    auto t1 = std::chrono::steady_clock::now();
    float timings[4] = {0,0,0,0};

    pRender->Clear(W, H, "color");
    pRender->Render(image.data(), W, H, "color"); 
    pRender->GetExecutionTime("CastRaySingleBlock", timings);

    auto t2 = std::chrono::steady_clock::now();
    float time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    debug("%s rendered in %.1f ms. %d kRays/s\n", "SDF Scene", time_ms, (int)((W * H) / time_ms));
    debug("CastRaySingleBlock took %.1f ms\n", timings[0]);

    LiteImage::SaveImage<uint32_t>("saves/liteRT_render_test.bmp", image);
  }

  void liteRT_grid_test()
  {
    auto scene = scene_complex_chair();
    ProceduralSdf reference_sdf(scene.first);
    reference_sdf.set_parameters(scene.second.p);
    AABB bbox = AABB({-1,-1,-1},{1,1,1});

    unsigned grid_size = 128;
    std::vector<float> grid_data(grid_size*grid_size*grid_size, 0);
    primitive_SDF_to_grid(reference_sdf, bbox, grid_data.data(), grid_size, 1);

    CameraSettings camera;
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    unsigned W = 1024, H = 1024;
    LiteImage::Image2D<uint32_t> image(W, H);

    auto pRender = MakeEyeRayShooterRenderer("CPU");
    pRender->SetAccelStruct(CreateSceneRT("BVH2Common", "cbvh_embree2", "SuperTreeletAlignedMerged4"));
    pRender->SetViewport(0,0, W, H);
    pRender->UpdateCamera(camera.get_view(), camera.get_proj(false));

    pRender->GetAccelStruct()->ClearGeom();
    pRender->GetAccelStruct()->AddGeom_SdfGrid({uint3(grid_size, grid_size, grid_size), grid_data.data()});
    pRender->GetAccelStruct()->ClearScene();
    pRender->GetAccelStruct()->AddInstance(0, LiteMath::float4x4());
    pRender->GetAccelStruct()->AddInstance(0, LiteMath::translate4x4(float3(-100,-100,-100)));
    pRender->GetAccelStruct()->CommitScene();
    pRender->CommitDeviceData();

    auto t1 = std::chrono::steady_clock::now();
    float timings[4] = {0,0,0,0};

    pRender->Clear(W, H, "color");
    pRender->Render(image.data(), W, H, "color"); 
    pRender->GetExecutionTime("CastRaySingleBlock", timings);

    auto t2 = std::chrono::steady_clock::now();
    float time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    debug("%s rendered in %.1f ms. %d kRays/s\n", "SDF 64x64x64 Grid", time_ms, (int)((W * H) / time_ms));
    debug("CastRaySingleBlock took %.1f ms\n", timings[0]);

    LiteImage::SaveImage<uint32_t>("saves/liteRT_grid_test.bmp", image);
  }

  void liteRT_octree_test()
  {
    auto scene = scene_complex_chair();
    ProceduralSdf reference_sdf(scene.first);
    reference_sdf.set_parameters(scene.second.p);
    AABB bbox = AABB({-1,-1,-1},{1,1,1});

    std::vector<float> data = {0};
    ProceduralSdf g_sdf({{SdfNodeType::OCTREE}});
    g_sdf.set_parameters(data);
    dynamic_cast<OctreeSdfNode*>(g_sdf.root)->construct( [&reference_sdf](const float3 &p) {
      /*logerr("sample %f %f %f",p.x,p.y,p.z);*/ return reference_sdf.get_distance(p); });

    SparseOctree &octree = dynamic_cast<OctreeSdfNode*>(g_sdf.root)->octree;

    CameraSettings camera;
    camera.origin = float3(0,0,3);
    camera.target = float3(0,0,0);
    camera.up = float3(0,1,0);

    unsigned W = 1024, H = 1024;
    LiteImage::Image2D<uint32_t> image(W, H);

    auto pRender = MakeEyeRayShooterRenderer("GPU");
    pRender->SetAccelStruct(CreateSceneRT("BVH2Common", "cbvh_embree2", "SuperTreeletAlignedMerged4"));
    pRender->SetViewport(0,0, W, H);
    pRender->UpdateCamera(camera.get_view(), camera.get_proj(false));

    pRender->GetAccelStruct()->ClearGeom();
    pRender->GetAccelStruct()->AddGeom_SdfOctree({(unsigned)(octree.get_nodes().size()), 
                                                  (const SdfOctreeNode *)octree.get_nodes().data()});
    pRender->GetAccelStruct()->ClearScene();
    pRender->GetAccelStruct()->AddInstance(0, LiteMath::float4x4());
    pRender->GetAccelStruct()->AddInstance(0, LiteMath::translate4x4(float3(-100,-100,-100)));
    pRender->GetAccelStruct()->CommitScene();
    pRender->CommitDeviceData();

    auto t1 = std::chrono::steady_clock::now();
    float timings[4] = {0,0,0,0};

    pRender->Clear(W, H, "color");
    pRender->Render(image.data(), W, H, "color"); 
    pRender->GetExecutionTime("CastRaySingleBlock", timings);

    auto t2 = std::chrono::steady_clock::now();
    float time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    debug("%s rendered in %.1f ms. %d kRays/s\n", "SDF 64x64x64 Grid", time_ms, (int)((W * H) / time_ms));
    debug("CastRaySingleBlock took %.1f ms\n", timings[0]);

    LiteImage::SaveImage<uint32_t>("saves/liteRT_octree_test.bmp", image);
  }

  void perform_benchmarks(const Block &blk)
  {
    //liteRT_grid_test();
    //liteRT_octree_test();
    //return;
    //sdf_octree_test();
    //return;
    std::string name = blk.get_string("name", "rendering");
    if (name == "rendering")
      benchmark_sdf_rendering(512, 16);
    else if (name == "nsdf")
      nsdf::neural_SDF_test();
    else if (name == "vox_net")
    {
      //voxNet_create_dataset(2500);
      //voxNet2D_test(); 
      //CNN_2_5D_create_dataset(25000);
      //CNN_2_5D_train();
      CNN_2_5D_test();
    }
    else if (name == "grid")
    {
      sdf_grid_test();
    }
    else if (name == "density_field")
    {
      // std::string obj_path = "./resources/mitsuba_data/meshes/sphere.obj";
      // auto model = dgen::load_obj(obj_path);
      
      // std::vector<float> sdf_model = df::pipeline(model);

      ////  Save sdf for next use
      // df::save_sdf(sdf_model, "sphere.sdf");

      auto sdf_model = df::readFile("sphere.sdf");

      int steps = 15;
      ProceduralSdf g_sdf({{SdfNodeType::GRID_32}});
      g_sdf.set_parameters(sdf_model);
      
      for (int i=0;i<steps;i++)
      {
        CameraSettings camera;
        camera.origin = float3(3*cos((2.0f*PI*i)/steps),0,3*sin((2.0f*PI*i)/steps));
        camera.target = float3(0,0,0);
        camera.up = float3(0,1,0);
        Texture t = render_sdf(g_sdf, camera, 512, 512, 4, SDFRenderMode::LAMBERT);
        engine::textureManager->save_png(t, "reconstructed_image_grid_"+std::to_string(i));
      }
    }
    else if (name == "QR")
    {
      std::vector<float> M(64 * 64, 0);

      for (int i = 0; i < 64 * 64; i++)
      {
        M[i] = (int)(1 + 10.f * rand() / (float)RAND_MAX);
      }

      std::vector<float> Q(64 * 64, 0), R(64 * 64, 0);

      interpolation::QR(M, 64, Q, R);

      std::vector<float> G = interpolation::mul_qr(Q, R, 64);

      std::cout << std::endl << "Deviation of the resulting decomposition from the original matrix: " << interpolation::matrix_norm(M, G) << std::endl << std::endl;
    }
    else
      benchmark_sdf_complex_optimization();
  }
}