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


  //TEST 5 SCALED BOX SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_5()
  {
    srand(time(NULL));
    debug("TEST 5. SCALED BOX SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 50000
            params:arr = {0.2,-0.1,0, 0.5, 1,1,1}
            structure:arr = {2,12,4}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,-0.15,-0.1,0.41, 1,1,1}    
            structure:arr = {2,12,4} 
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
    
    debug("  5.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  5.2. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 40)
      debug("passed\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 6 ROTATED BOX SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_6()
  {
    srand(time(NULL));
    debug("TEST 6. ROTATED BOX SDF RECONSTRUCTION\n");
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
    
    debug("  6.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug("  6.2. %-64s", "Perfect multi-view PSNR ");
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
    ProceduralSdf rot_box({std::vector<uint16_t>{11, 4}});

    debug("TEST 13. SDF ROTATION NODE\n");
    {
    int pcnt = rot_box.desc.get_total_params_count();
    debug(" 13.1. %-64s", "SDF instances are created with expected number of parameters ");
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

      debug(" 13.2. %-64s", "Distance to box correct");
      if (abs(d1) < 1e-6 && abs(d2) < 1e-6 && abs(d3) < 1e-6)
        debug("passed\n");
      else
        debug("FAILED %f %f %f\n", d1, d2, d3);
    }
  }

  void test_14()
  {
    debug("TEST 14. SDF CHAIR NODE\n");
    ProceduralSdf chair({std::vector<uint16_t>{SdfNodeType::CHAIR}});
    {
      int pcnt = chair.desc.get_total_params_count();
      debug(" 14.1. %-64s", "SDF instances are created with expected number of parameters ");
      if (pcnt == 6)
        debug("PASSED\n");
      else
        debug("FAILED %d\n", pcnt);
    }
    {
      std::vector<float> params = {0.2, 1, 0.8, 0.1, 0.8, 1.5};
      chair.set_parameters(params);
      std::vector<float> ddist,dpos = {0,0,0};
      debug(" 14.2. %-64s", "distances getting ");
      float d1 = chair.get_distance({0,sqrt(2),0},&ddist,&dpos);
      float d2 = chair.get_distance({0,sqrt(0.5),sqrt(0.5)});
      float d3 = chair.get_distance({1,0.5,0.5});
      debug("PASSED\n");

      /*debug(" 34.2. %-64s", "Distance to box correct");
      if (abs(d1) < 1e-6 && abs(d2) < 1e-6 && abs(d3) < 1e-6)
        debug("PASSED\n");
      else
        debug("FAILED %f %f %f\n", d1, d2, d3);*/
      
      /*float dist_1=0,dist_2=0;
      for (int i=0;i<std::min(ddist.size(), ddist_ref.size());i++)
        dist_1 += std::abs(ddist[i]-ddist_ref[i]);
      for (int i=0;i<std::min(dpos.size(), dpos_ref.size());i++)
        dist_2 += std::abs(dpos[i]-dpos_ref[i]);
      
      debug(" 20.3. %-64s", "Derivatives correct");
      if (ddist.size() == 1 && dist_1 < 1e-6 && dpos.size() == 3 && dist_2 < 1e-6)
        debug("PASSED\n");
      else
        debug("FAILED %d %d %d %d\n", ddist.size() == 1, dist_1 < 1e-6, dpos.size() == 3, dist_2 < 1e-6);*/
    }
  }

  //TEST 15 FIELD-BASED SPHERE SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void sdf_test_15()
  {
    srand(0);
    debug("TEST 15. FIELD-BASED SPHERE SDF RECONSTRUCTION\n");
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
            field:b = true
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

    debug(" 15.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("passed\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug(" 15.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 2 && res[0].structure.s[0] == 2 && res[0].structure.s[1] == 1)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug(" 15.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 4)
      debug("passed\n");
    else
      debug("FAILED\n");
    
    debug(" 15.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("passed\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug(" 15.5. %-64s", "Perfect multi-view PSNR ");
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
      sdf_test_11, sdf_test_12, sdf_test_13, sdf_test_15, sdf_test_15,
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