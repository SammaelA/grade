#include "upg.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/image_metrics.h"
#include "graphics_utils/modeling.h"
#include <unistd.h>

namespace upg
{

  //TEST 1 ONE TRIANGLE RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 80+ PSNR)
  void test_1()
  {
    srand(0);
    debug("TEST 1. ONE TRIANGLE SINGLE-VIEW RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 256
            reference_image_h:i = 256
            params:arr = {0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {1}
        } 
        view_0 {
            camera.origin:p3 = 0.000000, 0.000000, 3.000000
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
            params:arr = {0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81}    
            structure:arr = {1} 
        }
        step_0 {
            render_w:i = 256
            render_h:i = 256
            iterations:i = 100
            verbose:b = false
            save_intermediate_images:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
        //save_folder:s = "upg_triangle"
        //save_turntable:b = true
        save_turntable_hydra_settings {
            save_filename:s = "upg_triangle/result"
            image_count:i = 16
            rays_per_pixel:i = 512
            image_size:i2 = 1024, 1024
            distance:r = 2
            height:r = 0.5
            render_terrain:b = false
        }
        //save_reference_turntable:b = true
        save_reference_turntable_hydra_settings {
            save_filename:s = "upg_triangle/reference"
            image_count:i = 16
            rays_per_pixel:i = 512
            image_size:i2 = 1024, 1024
            distance:r = 2
            height:r = 0.5
            render_terrain:b = false
        }
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct(settings_blk);

    debug("  1.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("PASSED\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  1.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 1 && res[0].structure.s[0] == 1)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug("  1.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 9)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug("  1.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug("  1.5. %-64s", "Perfect one-view PSNR ");
    if (res[0].quality_ir > 80)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 80);

    debug("  1.6. %-64s", "Adequate turntable PSNR ");
    if (res[0].quality_synt > 10)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 10);

    ComplexModel mod;
    bool cm_res = create_model(res[0].structure, res[0].parameters, mod);
    debug("  1.7. %-64s", "Reconstructed model can be created ");
    if (cm_res)
      debug("PASSED\n");
    else
      debug("FAILED\n");

    debug("  1.8. %-64s", "Reconstructed model is valid ");
    if (mod.materials.size() == 1 && mod.models.size() == 1 && mod.models[0] && mod.models[0]->positions.size() == 3*3 && 
        mod.models[0]->normals.size() == 3*3 && mod.models[0]->colors.size() == 4*3)
      debug("PASSED\n");
    else
      debug("FAILED %d %d %d %d %d %d\n", mod.materials.size() == 1, mod.models.size() == 1, mod.models[0], 
            mod.models[0] ? (mod.models[0]->positions.size() == 3*3) : 0, 
            mod.models[0] ? (mod.models[0]->normals.size() == 3*3) : 0,
            mod.models[0] ? (mod.models[0]->colors.size() == 4*3) : 0);    
  }

  //TEST 2 ONE TRIANGLE RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_2()
  {
    srand(0);
    debug("TEST 2. ONE TRIANGLE MULTI-VIEW RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {1}
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
            params:arr = {0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81}    
            structure:arr = {1} 
        }
        step_0 {
            render_w:i = 512
            render_h:i = 512
            iterations:i = 500
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

    debug("  2.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug("  2.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug("  2.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 3 ONE TRIANGLE RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_3()
  {
    srand(0);
    debug("TEST 3. ONE TRIANGLE RECONSTRUCTION FROM MASK\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 256
            reference_image_h:i = 256
            params:arr = {0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {1}
        } 
        view_0 {
            camera.origin:p3 = 0.000000, 0.000000, 3.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
            mask:s = "saves/tests/test_3_input.png"
        }
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81}    
            structure:arr = {1} 
        }
        step_0 {
            render_w:i = 256
            render_h:i = 256
            iterations:i = 100
            verbose:b = false
            save_intermediate_images:b = false
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

    debug("  3.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("PASSED\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug("  3.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 1 && res[0].structure.s[0] == 1)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug("  3.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 9)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug("  3.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug("  3.5. %-64s", "Perfect one-view PSNR ");
    if (res[0].quality_ir > 80)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 80);
  }

  //image metrics
  void test_4()
  {
    srand(0);
    debug("TEST 4. IMAGE METRICS\n");
    Texture t = engine::textureManager->load_unnamed_tex("saves/tests/test_3_input.png", 1);
    
    float mae = ImageMetric::get(t, t, ImageMetric::MAE);
    debug("  4.1. %-64s", "MAE ");
    if (mae < 1e-6)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n",mae, 1e-6);
    
    float mse = ImageMetric::get(t, t, ImageMetric::MSE);
    debug("  4.2. %-64s", "MSE ");
    if (mse < 1e-6)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n",mse, 1e-6);

    float psnr = ImageMetric::get(t, t, ImageMetric::PSNR);
    debug("  4.3. %-64s", "PSNR ");
    if (psnr > 90 - 1e-6)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n",psnr, 90 - 1e-6);
    
    float iou = ImageMetric::get(t, t, ImageMetric::IOU);
    debug("  4.4. %-64s", "IOU ");
    if (iou > 1 - 1e-6)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n",iou, 1 - 1e-6);
  }

  //hydra visualization of reconstructed model
  void test_5()
  {
    srand(0);
    debug("TEST 5. HYDRA VISUALIZATION OF RECONSTRUCTION RESULTS\n");
    
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 256
            reference_image_h:i = 256
            params:arr = {0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {1}
        } 
        view_0 {
            camera.origin:p3 = 0.000000, 0.000000, 3.000000
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
            params:arr = {0,0,0, -1,0,0, 0,1,-1}   
            structure:arr = {1} 
        }
        step_0 {
            render_w:i = 256
            render_h:i = 256
            iterations:i = 1
            verbose:b = false
            save_intermediate_images:b = false
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
        save_folder:s = "tests/test_5"
        save_turntable:b = true
        save_turntable_hydra_settings {
            save_filename:s = "tests/test_5/result"
            image_count:i = 16
            rays_per_pixel:i = 256
            image_size:i2 = 256, 256
            distance:r = 2
            height:r = 0.5
            render_terrain:b = false
        }
        save_reference_turntable:b = true
        save_reference_turntable_hydra_settings {
            save_filename:s = "tests/test_5/reference"
            image_count:i = 16
            rays_per_pixel:i = 256
            image_size:i2 = 256, 256
            distance:r = 2
            height:r = 0.5
            render_terrain:b = false
        }
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct(settings_blk);

    sleep(1);

    bool same_size = true;
    float min_psnr_ref = 90;
    float min_psnr_res = 90;
    for (int i=0;i<16;i++)
    {
      char path[1024];
      sprintf(path, "%s-%04d.png", "saves/tests/test_5_ref/reference", i);
      Texture t1 = engine::textureManager->load_unnamed_tex(std::string(path), 1);
      sprintf(path, "%s-%04d.png", "saves/tests/test_5/reference", i);
      Texture t2 = engine::textureManager->load_unnamed_tex(std::string(path), 1);
      sprintf(path, "%s-%04d.png", "saves/tests/test_5_ref/result", i);
      Texture t3 = engine::textureManager->load_unnamed_tex(std::string(path), 1);
      sprintf(path, "%s-%04d.png", "saves/tests/test_5/result", i);
      Texture t4 = engine::textureManager->load_unnamed_tex(std::string(path), 1);

      if (t2.get_W() == 256 && t2.get_H() == 256 && t4.get_W() == 256 && t4.get_H() == 256)
      {
        min_psnr_ref = MIN(min_psnr_ref, ImageMetric::get(t1, t2, ImageMetric::PSNR));
        min_psnr_res = MIN(min_psnr_res, ImageMetric::get(t3, t4, ImageMetric::PSNR));
      }
      else
        same_size = false;
    }

    debug("  5.1. %-64s", "Images rendered ");
    if (same_size)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug("  5.2. %-64s", "Reference images match ");
    if (min_psnr_ref > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n",min_psnr_ref, 40);
    
    debug("  5.3. %-64s", "Result images match ");
    if (min_psnr_res > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n",min_psnr_res, 40);
  }

  //TEST 6 PARAMETERS PRESERVATION + OBJ SAVING
  //Create model from synthetic reference. Optimize with 0 learning rate
  //and starting parameters identical to reference. Save result to obj.
  //Assure that it saved exactly the same model as input one and got
  //perfect metrics for reconstruction
  void test_6()
  {
    srand(0);
    debug("TEST 6. PARAMETERS PRESERVATION + OBJ SAVING\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 256
            reference_image_h:i = 256
            params:arr = {0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {1}
        } 
        view_0 {
            camera.origin:p3 = 0.000000, 0.000000, 3.000000
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
            params:arr = {0,0,0, -1,0,0, 0,1,-1}    
            structure:arr = {1} 
        }
        step_0 {
            render_w:i = 256
            render_h:i = 256
            iterations:i = 100
            verbose:b = false
            save_intermediate_images:b = false
            learning_rate:r = 0.0
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
        save_model:b = true
        save_folder:s = "tests/test_6"
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct(settings_blk);
    
    debug("  6.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug("  6.2. %-64s", "Perfect one-view PSNR ");
    if (res[0].quality_ir > 80.0)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 80.0);

    debug("  6.3. %-64s", "Perfect turntable PSNR ");
    if (res[0].quality_synt > 80.0)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 80.0);

    sleep(1);
    Model *m;
    std::vector<float> ref_positions = {0,0,0, -1,0,0, 0,1,-1};
    std::vector<float> ref_normals = {-0.000000, -0.707107, -0.707107, -0.000000, -0.707107, -0.707107, -0.000000, -0.707107, -0.707107};
    std::vector<float> ref_tc = {0,0,0,1, 1,0,0,1, 0,1,0,1};
    auto match = [](std::vector<float> &v1, std::vector<float> &v2) -> bool
    {
      for (int i=0;i<v1.size();i++)
      {
        if (abs(v2[i] - v1[i]) > 1e-6)
          return false;
      }
      return true;
    };
    {
    m = model_loader::load_model_from_obj_directly("saves/tests/test_6/reconstructed_model.obj");
    debug("  6.4. %-64s", "Reconstructed model saved to obj file ");
    if (m)
      debug("PASSED\n");
    else 
      debug("FAILED\n");
    bool pos_match = m && (m->positions.size() == ref_positions.size()) && match(m->positions, ref_positions);
    bool norm_match = m && (m->normals.size() == ref_normals.size()) && match(m->normals, ref_normals);
    bool tc_match = m && (m->colors.size() == ref_tc.size()) && match(m->colors, ref_tc);
    debug("  6.5. %-64s", "Reconstructed model preserved ");
    if (pos_match && norm_match && tc_match)
      debug("PASSED\n");
    else 
      debug("FAILED %d(sz %d, %d) %d %d\n", pos_match, m->positions.size(), ref_positions.size(), norm_match, tc_match);  
    delete m;
    }
    {
    m = model_loader::load_model_from_obj_directly("saves/tests/test_6/reference_model.obj");
    debug("  6.6. %-64s", "Reference model saved to obj file ");
    if (m)
      debug("PASSED\n");
    else 
      debug("FAILED\n");
    bool pos_match = m && (m->positions.size() == ref_positions.size()) && match(m->positions, ref_positions);
    bool norm_match = m && (m->normals.size() == ref_normals.size()) && match(m->normals, ref_normals);
    bool tc_match = m && (m->colors.size() == ref_tc.size()) && match(m->colors, ref_tc);
    debug("  6.7. %-64s", "Reference model preserved ");
    if (pos_match && norm_match && tc_match)
      debug("PASSED\n");
    else 
      debug("FAILED %d %d %d\n", pos_match, norm_match, tc_match);  
    delete m;
    }
  }

  //TEST 7 USING OBJ REFERENCE MODEL
  void test_7()
  {
    srand(0);
    debug("TEST 7. USING OBJ REFERENCE MODEL\n");
    std::string settings = R""""(
    {
    input {
        model_reference {
            reference_image_w:i = 256
            reference_image_h:i = 256
            obj_filename:s = "saves/tests/test_7_ref/reference_model.obj"
        } 
        view_0 {
            camera.origin:p3 = 0.000000, 0.000000, 3.000000
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
            params:arr = {0,0,0, -1,0,0, 0,1,-1}    
            structure:arr = {1} 
        }
        step_0 {
            render_w:i = 256
            render_h:i = 256
            iterations:i = 100
            verbose:b = false
            save_intermediate_images:b = false
            learning_rate:r = 0.0
        }
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
        save_model:b = true
        save_folder:s = "tests/test_7"
    }
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    auto res = reconstruct(settings_blk);
    
    debug("  7.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug("  7.2. %-64s", "Perfect one-view PSNR ");
    if (res[0].quality_ir > 80.0)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 80.0);

    debug("  7.3. %-64s", "Perfect turntable PSNR ");
    if (res[0].quality_synt > 80.0)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 80.0);

    sleep(1);
    Model *m, *m_ref;
    auto match = [](std::vector<float> &v1, std::vector<float> &v2) -> bool
    {
      for (int i=0;i<v1.size();i++)
      {
        if (abs(v2[i] - v1[i]) > 1e-6)
          return false;
      }
      return true;
    };
    {
    m = model_loader::load_model_from_obj_directly("saves/tests/test_7/reconstructed_model.obj");
    m_ref = model_loader::load_model_from_obj_directly("saves/tests/test_7_ref/reconstructed_model.obj");
    debug("  7.4. %-64s", "Reconstructed model saved to obj file ");
    if (m)
      debug("PASSED\n");
    else 
      debug("FAILED\n");
    bool pos_match = m && (m->positions.size() == m_ref->positions.size()) && match(m->positions, m_ref->positions);
    bool norm_match = m && (m->normals.size() == m_ref->normals.size()) && match(m->normals, m_ref->normals);
    bool tc_match = m && (m->colors.size() == m_ref->colors.size()) && match(m->colors, m_ref->colors);
    debug("  7.5. %-64s", "Reconstructed model preserved ");
    if (pos_match && norm_match && tc_match)
      debug("PASSED\n");
    else 
      debug("FAILED %d(sz %d, %d) %d %d\n", pos_match, m->positions.size(), m_ref->positions.size(), norm_match, tc_match);  
    delete m;
    delete m_ref;
    }
    {
    m = model_loader::load_model_from_obj_directly("saves/tests/test_7/reference_model.obj");
    m_ref = model_loader::load_model_from_obj_directly("saves/tests/test_7_ref/reference_model.obj");
    debug("  7.6. %-64s", "Reference model saved to obj file ");
    if (m)
      debug("PASSED\n");
    else 
      debug("FAILED\n");
    bool pos_match = m && (m->positions.size() == m_ref->positions.size()) && match(m->positions, m_ref->positions);
    bool norm_match = m && (m->normals.size() == m_ref->normals.size()) && match(m->normals, m_ref->normals);
    bool tc_match = m && (m->colors.size() == m_ref->colors.size()) && match(m->colors, m_ref->colors);
    debug("  7.7. %-64s", "Reference model preserved ");
    if (pos_match && norm_match && tc_match)
      debug("PASSED\n");
    else 
      debug("FAILED %d %d %d\n", pos_match, norm_match, tc_match);  
    delete m;
    delete m_ref;
    }
  }

  void perform_tests()
  {
    test_1();
    test_2();
    test_3();
    test_4();
    test_5();
    test_6();
    test_7();
  }
};