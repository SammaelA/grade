#include "upg.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/image_metrics.h"
#include "graphics_utils/modeling.h"
#include "generation.h"
#include "sdf_node.h"
#include <unistd.h>
#include <functional>

namespace upg
{
  //TEST 1 SDF NODES
  //tests the most basic SDF nodes functions
  //such as distance, merging and moving nodes
  //an putting derivative is right (root,left,right) order
  void sdf_test_1()
  {
    ProceduralSdf one_circle({std::vector<uint16_t>{1}});
    ProceduralSdf moved_circle({std::vector<uint16_t>{2,1}});
    ProceduralSdf two_circles({std::vector<uint16_t>{3,2,1,2,1}});

    debug("TEST 1. SDF NODES\n");
    {
    int pcnt_1 = one_circle.desc.get_total_params_count();
    int pcnt_2 = moved_circle.desc.get_total_params_count();
    int pcnt_3 = two_circles.desc.get_total_params_count();
    debug("  1.1. %-64s", "SDF instances are created with expected number of parameters ");
    if (pcnt_1 == 1 && pcnt_2 == 4 && pcnt_3 == 8)
      debug("passed\n");
    else
      debug("FAILED %d %d %d\n", pcnt_1, pcnt_2, pcnt_3);
    }
    {
      std::vector<float> params = {0.5};
      one_circle.set_parameters(params);
      std::vector<float> ddist,dpos = {0,0,0};
      std::vector<float> ddist_ref = {-1}, dpos_ref = {0,1,0};
      float d1 = one_circle.get_distance({0,1,0},&ddist,&dpos);
      float d2 = one_circle.get_distance({1,0,0});

      debug("  1.2. %-64s", "Distance to circle correct");
      if (abs(d1 - (0.5)) < 1e-6 && abs(d2 - (0.5)) < 1e-6)
        debug("passed\n");
      else
        debug("FAILED %f %f\n", d1, d2);
      
      float dist_1=0,dist_2=0;
      for (int i=0;i<std::min(ddist.size(), ddist_ref.size());i++)
        dist_1 += std::abs(ddist[i]-ddist_ref[i]);
      for (int i=0;i<std::min(dpos.size(), dpos_ref.size());i++)
        dist_2 += std::abs(dpos[i]-dpos_ref[i]);
      
      debug("  1.3. %-64s", "Derivatives correct");
      if (ddist.size() == 1 && dist_1 < 1e-6 && dpos.size() == 3 && dist_2 < 1e-6)
        debug("passed\n");
      else
        debug("FAILED %d %d %d %d\n", ddist.size() == 1, dist_1 < 1e-6, dpos.size() == 3, dist_2 < 1e-6);
    }
  }

    //TEST 2 SPHERE SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_2()
  {
    srand(0);
    debug("TEST 2. SPHERE SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0,0,0,1}
            structure:arr = {2,1}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,0.2,-0.1,0.7}    
            structure:arr = {2,1} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug("  2.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  2.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 1)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  2.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 4)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  2.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  2.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 3 SPHERES SDF RECONSTRUCTION MEMETIC
  void sdf_test_3()
  {
    srand(0);
    debug("TEST 3. SPHERES SDF RECONSTRUCTION MEMETIC\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.6,0.6,0.6,0.5, -0.6,0.6,0.6,0.5, 0.6,-0.6,0.6,0.5, -0.6,-0.6,0.6,0.5, 
                          0.6,0,0,0.5, -0.6,0,0,0.5, 0,-0.6,0,0.5, 0,0.6,0,0.5}
            structure:arr = {3,3,3,2,1,2,1,3,2,1,2,1,3,3,2,1,2,1,3,2,1,2,1}
        } 
    }
    generator {

    }
    optimization {
        start {
            //params:arr = {0.5,0.1,-0.1,0.5, -0.5,0.06,-0.09,0.54, 0.1,-0.66,-0.09,0.45, 0.1,0.57,0,0.51}
            params:arr = {0.6,0.6,0.6,0.5, -0.6,0.6,0.6,0.5, 0.6,-0.6,0.6,0.5, -0.6,-0.6,0.6,0.5, 
                          0.6,0,0,0.5, -0.6,0,0,0.5, 0,-0.6,0,0.5, 0,0.6,0,0.5}
            structure:arr = {3,3,3,2,1,2,1,3,2,1,2,1,3,3,2,1,2,1,3,2,1,2,1}
        }
        step_0 {
            optimizer_name:s = "memetic"
            verbose:b = false
        }
        step_1 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug("  3.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    bool str_eq = true;
    std::vector<uint16_t> ref_struct = {3,3,3,2,1,2,1,3,2,1,2,1,3,3,2,1,2,1,3,2,1,2,1};
    for (int i=0;i<std::min(res[0].structure.s.size(), ref_struct.size());i++)
       str_eq = str_eq && (res[0].structure.s[i] == ref_struct[i]);
    debug("  3.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == ref_struct.size() && str_eq)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  3.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 4*8)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  3.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  3.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 4 BOX SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_4()
  {
    srand(time(NULL));
    debug("TEST 4. BOX SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0,0.5,0.5,0.5}
            structure:arr = {2,4}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,-0.15,-0.1,0.41,0.43,0.61}    
            structure:arr = {2,4} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug("  4.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  4.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 4)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  4.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 6)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  4.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  4.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }


  //TEST 5 COMPLEX DETAIL RECONSTRUCTION MEMETIC
  void sdf_test_5()
  {
    srand(time(NULL));
    debug("TEST 5. COMPLEX DETAIL RECONSTRUCTION MEMETIC\n");
    debug("TEMPORARY DISABLED\n");
    return;
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
        }
        step_1 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);

    int cnt_x = 8;
    int cnt_z = 4;
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

    settings_blk.get_block("input")->get_block("synthetic_reference")->set_arr("structure", structure);
    settings_blk.get_block("input")->get_block("synthetic_reference")->set_arr("params", params);
    settings_blk.get_block("optimization")->get_block("start")->set_arr("structure", structure);

    //return;
    auto res = reconstruct_sdf(settings_blk);
    
    debug("  5.1. %-64s", "Low optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
  }

  //TEST 6 ROTATING BODY RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void sdf_test_6()
  {
    srand(0);
    debug("TEST 6. ROTATING BODY MULTI-VIEW RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0.2, 0.21, 0.23, 0.26, 0.3, 0.35, 0.41, 0.48}
            structure:arr = {7}
        } 
        view_0 {
            camera.origin:p3 = 2.000000, 0.500000, 2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_1 {
            camera.origin:p3 = 2.000000, -0.500000, -2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_2 {
            camera.origin:p3 = -2.000000, 0.500000, 2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_3 {
            camera.origin:p3 = -2.000000, -0.500000, -2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4}
            structure:arr = {7}
        }
        step_0 {
            render_w:i = 512
            render_h:i = 512
            iterations:i = 1000
            verbose:b = false
            save_intermediate_images:b = false
            learning_rate:r = 0.003
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct(settings_blk);

    debug("  6.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug("  6.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug("  6.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 7 ROUND BOX SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_7()
  {
    srand(time(NULL));
    debug("TEST 7. ROUND BOX SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0,0.5,0.5,0.5, 0.1}
            structure:arr = {2,6}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.17,-0.15,-0.03,0.46,0.49,0.51, 0.1}    
            structure:arr = {2,6} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug("  7.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  7.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 6)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  7.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 7)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  7.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 2e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 2e-5);
    
    debug("  7.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 33)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 33.0f);
  }

  void sdf_test_8()
  {
    srand(time(NULL));
    debug("TEST 8. PRISM SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0,0.5,0.5}
            structure:arr = {2,7}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,-0.15,-0.1,0.41,0.43}    
            structure:arr = {2,7} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug("  8.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  8.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 7)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  8.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 5)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  8.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  8.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void sdf_test_9()
  {
    srand(time(NULL));
    debug("TEST 9. CYLINDER SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0,0.5,0.5}
            structure:arr = {2,5}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,-0.15,-0.1,0.41,0.43}    
            structure:arr = {2,5} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug("  9.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  9.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 5)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  9.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 5)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug("  9.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  9.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void sdf_test_10()
  {
    srand(time(NULL));
    debug("TEST 10. CONE SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0,0.5,0.5,1}
            structure:arr = {2,8}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,-0.15,-0.1,0.51,0.43,1.07}    
            structure:arr = {2,8} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug(" 10.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug(" 10.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 8)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug(" 10.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 6)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug(" 10.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 2e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 2e-5);
    
    debug(" 10.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 30)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 30.0f);
  }


  //TEST 11 INTERSECTION SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_11()
  {
    srand(time(NULL));
    debug("TEST 11. INTERSECTION SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0,0,0,0.6,0.3,0.6,  0,0,0,0.5}
            structure:arr = {9, 2,4, 2,1}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {-0.08, 0.1, 0.1, 0.53, 0.27, 0.67, 0.11, 0.08, -0.09, 0.51}    
            structure:arr = {9, 2,4, 2,1}
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug(" 11.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug(" 11.2. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 12 SUBTRACT SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_12()
  {
    srand(time(NULL));
    debug("TEST 12. SUBTRACT SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0,1,0,0.6,0.3,0.6,  0,1,0,0.5}
            structure:arr = {10, 2,4, 2,1}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {-0.08, 0.91, 0.1, 0.53, 0.27, 0.67, 0.11, 0.98, -0.09, 0.51}    
            structure:arr = {10, 2,4, 2,1}
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);

    debug(" 12.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug(" 12.2. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void sdf_test_13()
  {
    srand(0);
    debug("TEST 13. COMPLEX MULTI-ROTATION\n");
    std::string settings = R""""(
    {
    params:arr = {1, 0, 0, 0, 1, 1, 0,0,0, 1,0,0, 0,1,-1}    
    structure:arr = {9, 1} 
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    ComplexModel m;
    bool res;
    debug(" 13.1. %-64s", "Creating model from params and structure ");
    if ((res = create_model_from_block(settings_blk, m)))
    {
      debug("passed\n");
    }
    else
    {
      debug("FAILED %d\n", res);
      return;
    }
    debug(" 13.2. %-64s", "Compare results and expectations ");
    Model *model = m.models[0];
    std::vector<float> p = model->positions;
    res = (p.size() == 9 * MESH_REPEATS);
    if (res)
    {
      debug("passed\n");
    }
    else
    {
      debug("FAILED %d != %d\n", p.size(), 9 * MESH_REPEATS);
    }
  }

  void sdf_test_14()
  {
    srand(0);
    debug("TEST 14. PRISM MULTI-VIEW RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0.23, 0.26, 0.3, 0.35, 0.41, 0.2}
            structure:arr = {13}
        } 
        view_0 {
            camera.origin:p3 = 2.000000, 0.500000, 2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_1 {
            camera.origin:p3 = 2.000000, -0.500000, -2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_2 {
            camera.origin:p3 = -2.000000, 0.500000, 2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_3 {
            camera.origin:p3 = -2.000000, -0.500000, -2.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1, 0.1, 0.3, 0.3, 0.25, 0.15}
            structure:arr = {13}
        }
        step_0 {
            render_w:i = 512
            render_h:i = 512
            iterations:i = 1000
            verbose:b = false
            save_intermediate_images:b = false
            learning_rate:r = 0.003
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct(settings_blk);

    debug(" 14.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug(" 14.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug(" 14.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void sdf_test_15()
  {
    ProceduralSdf rot_box({std::vector<uint16_t>{11, 4}});

    debug("TEST 15. SDF ROTATION NODE\n");
    {
    int pcnt = rot_box.desc.get_total_params_count();
    debug(" 15.1. %-64s", "SDF instances are created with expected number of parameters ");
    if (pcnt == 6)
      debug("passed\n");
    else
      debug("FAILED %d\n", pcnt);
    }
    {
      std::vector<float> params = {0, 0, PI/4, 1, 1, 1};
      rot_box.set_parameters(params);
      std::vector<float> ddist,dpos = {0,0,0};
      std::vector<float> ddist_ref = {-1}, dpos_ref = {0,1,0};
      float d1 = rot_box.get_distance({0,sqrt(2),0},&ddist,&dpos);
      float d2 = rot_box.get_distance({0,sqrt(0.5),sqrt(0.5)});
      float d3 = rot_box.get_distance({1,0.5,0.5});

      debug(" 15.2. %-64s", "Distance to box correct");
      if (abs(d1) < 1e-6 && abs(d2) < 1e-6 && abs(d3) < 1e-6)
        debug("passed\n");
      else
        debug("FAILED %f %f %f\n", d1, d2, d3);
    }
  }

  //TEST 16 ROTATED BOX SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_16()
  {
    srand(time(NULL));
    debug("TEST 16. ROTATED BOX SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0, 0,0,0.7, 0.3,0.8,0.3}
            structure:arr = {2,11,4}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,-0.15,-0.1, 0.05,0.025,0.65, 0.2,0.9,0.27}    
            structure:arr = {2,11,4} 
        }
        step_0 {
            learning_rate:r = 0.003
            iterations:i = 1000
            verbose:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct_sdf(settings_blk);
    
    debug(" 16.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug(" 16.2. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void perform_tests_sdf_reconstruction(const std::vector<int> &test_ids)
  {
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      sdf_test_1,  sdf_test_2,  sdf_test_3,  sdf_test_4,  sdf_test_5,
      sdf_test_6,  sdf_test_7,  sdf_test_8,  sdf_test_9,  sdf_test_10,
      sdf_test_11, sdf_test_12, sdf_test_13, sdf_test_14, sdf_test_15,
      sdf_test_16
    };

    if (tests.empty())
    {
      tests.resize(test_functions.size());
      for (int i=0;i<test_functions.size();i++)
        tests[i] = i+1;
    }

    for (int i=0;i<80;i++)
      debug("#");
    debug("\nSDF RECONSTRUCTION TESTS\n");
    for (int i=0;i<80;i++)
      debug("#");
    debug("\n");
    
    for (int i : tests)
    {
      assert(i > 0 && i <= test_functions.size());
      test_functions[i-1]();
    }
  }
};