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
    desc.second.p = {0,0.0,0, 0.5,0.05,0.5, 0,0.45,-0.45, 0.5,0.45,0.05,
                     0.4,-0.3,0.4, 0.3,0.05,  -0.4,-0.3,0.4, 0.3,0.05,  0.4,-0.3,-0.4, 0.3,0.05,  -0.4,-0.3,-0.4, 0.3,0.05};
    return desc;    
  }

  SceneDesc scene_subtraction()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::SUBTRACT, SdfNodeType::MOVE, SdfNodeType::SPHERE, SdfNodeType::MOVE, SdfNodeType::SPHERE};
    desc.second.p = {0,0,0,0.6, 0.5,0.0,0.5,0.4};
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
    camera.origin = glm::vec3(0,0,3);
    camera.target = glm::vec3(0,0,0);
    camera.up = glm::vec3(0,1,0);

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
              glm::vec3 p = {(k+urand())/vox_size, (j+urand())/vox_size, (i+urand())/vox_size};
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
      cameras[i].origin = glm::vec3(d*cos(phi)*cos(psi), d*sin(psi), d*sin(phi)*cos(psi));
      cameras[i].target = glm::vec3(0,0,0);
      cameras[i].up = glm::vec3(-d*cos(phi)*sin(psi), d*cos(psi), -d*sin(phi)*sin(psi));
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

    nn::save_dataset("saves/CNN_2.5D_dataset.bin", &dataset);
  }
  
  void CNN_2_5D_test()
  {
    
    nn::Dataset dataset;
    nn::load_dataset("saves/CNN_2.5D_dataset.bin", &dataset);

    /*
    for (int i=0;i<100*64*64*8;i+=64*64*8)
    {
      for (int j=0;j<64;j+=2)
      {
        debug("[");
        for (int k=0;k<64;k+=2)
        {
          debug("%.2f ", dataset.train_data[i + 64*j + k]);
        }   
        debug("]\n");     
      }
      debug("\n\n");
    }
    return;
    */
    //for (int i=0;i<dataset.train_elements;i++)
    //  logerr("label %d %f %f", i, dataset.train_labels[2*i+0], dataset.train_labels[2*i+1]);

    nn::TensorProcessor::init("GPU");
    nn::NeuralNetwork net;
/*
    net.add_layer(std::make_shared<nn::Conv2DLayer>(64,64,8, 8, 5), nn::Initializer::GlorotNormal);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    //net.add_layer(std::make_shared<nn::BatchNormLayer>(), nn::Initializer::BatchNorm);
    net.add_layer(std::make_shared<nn::MaxPoolingLayer>(60, 60, 8));

    net.add_layer(std::make_shared<nn::Conv2DLayer>(30,30,8, 16), nn::Initializer::GlorotNormal);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    //net.add_layer(std::make_shared<nn::BatchNormLayer>(), nn::Initializer::BatchNorm);
    net.add_layer(std::make_shared<nn::MaxPoolingLayer>(28, 28, 16));

    net.add_layer(std::make_shared<nn::Conv2DLayer>(14,14,16, 16), nn::Initializer::GlorotNormal);
    net.add_layer(std::make_shared<nn::ReLULayer>());
    //net.add_layer(std::make_shared<nn::BatchNormLayer>(), nn::Initializer::BatchNorm);
    net.add_layer(std::make_shared<nn::MaxPoolingLayer>(12, 12, 16));

    net.add_layer(std::make_shared<nn::FlattenLayer>(6,6,16));

    net.add_layer(std::make_shared<nn::DenseLayer>(6*6*16, 256), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());

    net.add_layer(std::make_shared<nn::DenseLayer>(256, 256), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());

    net.add_layer(std::make_shared<nn::DenseLayer>(256, 128), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::ReLULayer>());

    net.add_layer(std::make_shared<nn::DenseLayer>(128, 2), nn::Initializer::He);
    net.add_layer(std::make_shared<nn::SoftMaxLayer>());
*/
    net.add_layer(std::make_shared<nn::Conv2DLayer>(64,64,8,8), nn::Initializer::GlorotNormal);
    //net.add_layer(std::make_shared<nn::ReLULayer>());
    //net.add_layer(std::make_shared<nn::BatchNormLayer>(), nn::Initializer::BatchNorm);
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
    net.train(dataset.train_data.data(), dataset.train_labels.data(), dataset.train_elements, 100, 250, true,
              nn::OptimizerRMSProp(1e-3f, 0.999, 1e-6, true), 
              nn::Loss::CrossEntropy, nn::Metric::Recall, true);
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
    GridSdfNode grid(vox_size, bbox);
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
              glm::vec3 p = {(k+0.5)/vox_size, (j+0.5)/vox_size, (i+0.5)/vox_size};
              p = bbox.size()*p + bbox.min_pos;
              d += sdf.get_distance(p);
            }
            grid.set_voxel(glm::uvec3(k,j,i), d/samples);
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
          glm::vec3 p = {(x+urand())/vox_size, (y+urand())/vox_size, (z+urand())/vox_size};
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
    ProceduralSdf g_sdf({{SdfNodeType::GRID}});
    g_sdf.set_parameters(data);
    for (int i=0;i<steps;i++)
    {
      CameraSettings camera;
      camera.origin = glm::vec3(3*cos((2.0f*PI*i)/steps),0,3*sin((2.0f*PI*i)/steps));
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(g_sdf, camera, 512, 512, 16, SDFRenderMode::LAMBERT);
      engine::textureManager->save_png(t, "dataset_image_grid_"+std::to_string(i));
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
    scenes["Subtraction"] = scene_subtraction();

    CameraSettings camera;
    camera.origin = glm::vec3(0,0,3);
    camera.target = glm::vec3(0,0,0);
    camera.up = glm::vec3(0,1,0);

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

  void perform_benchmarks(const Block &blk)
  {
    std::string name = blk.get_string("name", "rendering");
    if (name == "rendering")
      benchmark_sdf_rendering(512, 16);
    else if (name == "nsdf")
      nsdf::neural_SDF_test();
    else if (name == "vox_net")
    {
      //voxNet_create_dataset(2500);
      //voxNet2D_test(); 
      //CNN_2_5D_create_dataset(2000);
      CNN_2_5D_test();
    }
    else if (name == "grid")
    {
      sdf_grid_test();
    }
    else
      benchmark_sdf_complex_optimization();
  }
}