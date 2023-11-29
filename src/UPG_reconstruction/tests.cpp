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

  void test_8()
  {
    srand(0);
    debug("TEST 8. SIMPLE TRIANGLE SHIFTING GEN TREE\n");
    std::string settings = R""""(
    {
    params:arr = {-1,1,3, 0,0,0, -1,0,0, 0,1,-1}    
    structure:arr = {3, 1} 
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    std::vector<float> tr = {0, 0, 0, -1, 0, 0, 0, 1, -1};
    std::vector<float> shift = {-1, 1, 3};
    ComplexModel m;
    bool res;
    debug("  8.1. %-64s", "Creating model from params and structure ");
    if ((res = create_model_from_block(settings_blk, m)))
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d\n", res);
      return;
    }
    debug("  8.2. %-64s", "Compare results and expectations ");
    Model *model = m.models[0];
    std::vector<float> p = model->positions;
    res = (p.size() == tr.size());
    if (!res)
    {
      debug("FAILED %d != %d\n", p.size(), tr.size());
      return;
    }
    for (int i = 0; res && i < 9; ++i)
    {
      res = (p[i] == tr[i] + shift[i % 3]);
    }
    if (res)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %f %f %f, %f %f %f, %f %f %f\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
    }
  }

  void test_9()
  {
    srand(0);
    debug("TEST 9. SIMPLE TRIANGLE SCALING GEN TREE\n");
    std::string settings = R""""(
    {
    params:arr = {2,0.5,1, 0,0,0, -1,0,0, 0,1,-1}    
    structure:arr = {2, 1} 
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    std::vector<float> tr = {0, 0, 0, -1, 0, 0, 0, 1, -1};
    std::vector<float> sc = {2, 0.5, 1};
    ComplexModel m;
    bool res;
    debug("  9.1. %-64s", "Creating model from params and structure ");
    if ((res = create_model_from_block(settings_blk, m)))
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d\n", res);
      return;
    }
    debug("  9.2. %-64s", "Compare results and expectations ");
    Model *model = m.models[0];
    std::vector<float> p = model->positions;
    res = (p.size() == tr.size());
    if (!res)
    {
      debug("FAILED %d != %d\n", p.size(), tr.size());
      return;
    }
    for (int i = 0; res && i < 9; ++i)
    {
      res = (p[i] == tr[i] * sc[i % 3]);
    }
    if (res)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %f %f %f, %f %f %f, %f %f %f\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
    }
  }

  void test_10()
  {
    srand(0);
    debug("TEST 10. SIMPLE TRIANGLE ROTATING GEN TREE\n");//180 graduses
    std::string settings = R""""(
    {
    params:arr = {1,0,0,3.14159265, 0,0,0, -1,0,0, 0,1,-1}    
    structure:arr = {4, 1} 
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    std::vector<float> tr = {0, 0, 0, -1, 0, 0, 0, -1, 1};
    ComplexModel m;
    bool res;
    debug(" 10.1. %-64s", "Creating model from params and structure ");
    if ((res = create_model_from_block(settings_blk, m)))
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d\n", res);
      return;
    }
    debug(" 10.2. %-64s", "Compare results and expectations ");
    Model *model = m.models[0];
    std::vector<float> p = model->positions;
    res = (p.size() == tr.size());
    if (!res)
    {
      debug("FAILED %d != %d\n", p.size(), tr.size());
      return;
    }
    for (int i = 0; res && i < 9; ++i)
    {
      res = ((p[i] - tr[i] < 1e-5 && p[i] >= tr[i]) || (tr[i] - p[i] < 1e-5 && p[i] <= tr[i]));
    }
    if (res)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %f %f %f, %f %f %f, %f %f %f\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
    }
  }

  void test_11()
  {
    srand(0);
    debug("TEST 11. COMPLEX TRIANGLE ROTATING GEN TREE\n");
    std::string settings = R""""(
    {
    params:arr = {0.1,0.2,0.3,4, 0.1,0.2,0.3,-4, 0,0,0, -1,0,0, 0,1,-1}    
    structure:arr = {4, 4, 1} 
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    std::vector<float> tr = {0, 0, 0, -1, 0, 0, 0, 1, -1};
    ComplexModel m;
    bool res;
    debug(" 11.1. %-64s", "Creating model from params and structure ");
    if ((res = create_model_from_block(settings_blk, m)))
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d\n", res);
      return;
    }
    debug(" 11.2. %-64s", "Compare results and expectations ");
    Model *model = m.models[0];
    std::vector<float> p = model->positions;
    res = (p.size() == tr.size());
    if (!res)
    {
      debug("FAILED %d != %d\n", p.size(), tr.size());
    }
    for (int i = 0; res && i < 9; ++i)
    {
      res = ((p[i] - tr[i] < 1e-5 && p[i] >= tr[i]) || (tr[i] - p[i] < 1e-5 && p[i] <= tr[i]));
    }
    if (res)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %f %f %f, %f %f %f, %f %f %f\n", p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
    }
  }

  void test_12()
  {
    srand(0);
    debug("TEST 12. COMPLEX GEN TREE\n");
    std::string settings = R""""(
    {
    params:arr = {0.5,2,1.5, 0,0,0, -1,0,0, 0,1,-1, 0.5,0,0.2,3, -0.5,1.5,0, 0,0,0, -1,0,0, 0,1,-1}    
    structure:arr = {6, 2, 1, 4, 3, 1} 
    }
      )"""";
    Block settings_blk;
    load_block_from_string(settings, settings_blk);
    ComplexModel m;
    bool res;
    debug(" 12.1. %-64s", "Creating model from params and structure ");
    if ((res = create_model_from_block(settings_blk, m)))
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d\n", res);
      return;
    }
    debug(" 12.2. %-64s", "Compare results and expectations ");
    Model *model = m.models[0];
    std::vector<float> p = model->positions;
    res = (p.size() == 18);
    if (res)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d != %d\n", p.size(), 18);
    }
  }

  void test_13()
  {
    debug("TEST 13. TREE JACOBIAN CALCULATION\n");

    upg::UPGStructure structure;
    structure.s = {3, 2, 1};
    upg::UPGParametersRaw params;
    params.p = {10000,10000,10000,  10,10,10,  1,2,3,4,5,6,7,8,9};
    upg::MeshGenInstance gen(structure);
    upg::UniversalGenJacobian jac;
    auto mesh = gen.generate(params.p, &jac);

    bool res;
    debug(" 13.1. %-64s", "Model and jacobian created, have right size ");
    if (mesh.pos.size() == 9 && jac.get_xn() == 9 && jac.get_yn() == 15)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d %d %d\n", mesh.pos.size(), jac.get_xn(), jac.get_yn());
    }

    std::vector<float> reference_jac = {
      1.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 
      0.0000, 1.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 1.0000, 0.0000, 
      0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 1.0000, 0.0000, 0.0000, 1.0000,
      1.0000, 0.0000, 0.0000, 4.0000, 0.0000, 0.0000, 7.0000, 0.0000, 0.0000, 
      0.0000, 2.0000, 0.0000, 0.0000, 5.0000, 0.0000, 0.0000, 8.0000, 0.0000, 
      0.0000, 0.0000, 3.0000, 0.0000, 0.0000, 6.0000, 0.0000, 0.0000, 9.0000,
      10.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
      0.0000, 10.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
      0.0000, 0.0000, 10.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
      0.0000, 0.0000, 0.0000, 10.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
      0.0000, 0.0000, 0.0000, 0.0000, 10.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 10.0000, 0.0000, 0.0000, 0.0000, 
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 10.0000, 0.0000, 0.0000, 
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 10.0000, 0.0000, 
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 10.0000};

    float diff = 0;
    for (int i=0;i<jac.get_yn();i++)
    {
      for (int j=0;j<jac.get_xn();j++)
        diff += abs(reference_jac[i*jac.get_xn() + j] - jac.at(i,j));
    }
    debug(" 13.2. %-64s", "Jacobian is correct ");
    if (diff < 1e-7)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %f > %f\n", diff, 1e-7);
    }
  }

  //TEST 14 MOVE TRIANGLE RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_14()
  {
    srand(0);
    debug("TEST 14. MOVE TRIANGLE RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0.2,0.2,0.2, 0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {3, 1}
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
            params:arr = {0.15,0.23,0.25, 0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81}    
            structure:arr = {3, 1} 
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

    debug(" 14.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug(" 14.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug(" 14.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 15 SCALE TRIANGLE RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_15()
  {
    srand(0);
    debug("TEST 15. SCALE TRIANGLE RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0.8,1.0,1.2, 0,0,0, -1,0,0, 0,1,-1}
            structure:arr = {2, 1}
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
            params:arr = {1.2,1,0.8, 0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81}    
            structure:arr = {2, 1} 
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

    debug(" 15.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug(" 15.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug(" 15.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void test_16()
  {
    debug("TEST 16. STACKED CUBES MODEL AND JACOBIAN\n");

    upg::UPGStructure structure;
    structure.s = {6,  6, 3,2,5, 3,2,5,  3,2,5};
    upg::UPGParametersRaw params;
    params.p = {0,0,0, 1,1,1,   0,1,0, 0.5,0.5,0.5,  0,1.5,0, 0.25,0.25,0.25};
    upg::MeshGenInstance gen(structure);
    upg::UniversalGenJacobian jac;
    auto mesh = gen.generate(params.p, &jac);

    bool res;
    debug(" 16.1. %-64s", "Model and jacobian created, have right size ");
    if (mesh.pos.size() == 3*6*2*3*3 && jac.get_xn() == mesh.pos.size() && jac.get_yn() == params.p.size())
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d %d %d\n", mesh.pos.size(), jac.get_xn(), jac.get_yn());
    }

    std::vector<float> reference_model = {
      0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 1.000, 0.000, 0.000, 
      1.000, 1.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000, 0.000, 
      0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 1.000, 0.000, 
      0.000, 1.000, 1.000, 0.000, 1.000, 0.000, 0.000, 0.000, 1.000, 
      0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, 1.000, 
      1.000, 0.000, 1.000, 0.000, 0.000, 1.000, 1.000, 0.000, 0.000, 
      1.000, 1.000, 1.000, 0.000, 1.000, 1.000, 1.000, 0.000, 1.000, 
      0.000, 0.000, 1.000, 1.000, 0.000, 1.000, 0.000, 1.000, 1.000, 
      1.000, 1.000, 1.000, 1.000, 0.000, 1.000, 1.000, 1.000, 0.000, 
      1.000, 0.000, 0.000, 1.000, 1.000, 0.000, 1.000, 0.000, 1.000, 
      1.000, 1.000, 1.000, 1.000, 1.000, 0.000, 0.000, 1.000, 1.000,
      0.000, 1.000, 0.000, 0.000, 1.000, 1.000, 1.000, 1.000, 0.000,
      
      0.000, 1.000, 0.000, 0.000, 1.500, 0.000, 0.500, 1.000, 0.000,
      0.500, 1.500, 0.000, 0.500, 1.000, 0.000, 0.000, 1.500, 0.000, 
      0.000, 1.000, 0.000, 0.000, 1.000, 0.500, 0.000, 1.500, 0.000, 
      0.000, 1.500, 0.500, 0.000, 1.500, 0.000, 0.000, 1.000, 0.500, 
      0.000, 1.000, 0.000, 0.500, 1.000, 0.000, 0.000, 1.000, 0.500, 
      0.500, 1.000, 0.500, 0.000, 1.000, 0.500, 0.500, 1.000, 0.000, 
      0.500, 1.500, 0.500, 0.000, 1.500, 0.500, 0.500, 1.000, 0.500, 
      0.000, 1.000, 0.500, 0.500, 1.000, 0.500, 0.000, 1.500, 0.500, 
      0.500, 1.500, 0.500, 0.500, 1.000, 0.500, 0.500, 1.500, 0.000, 
      0.500, 1.000, 0.000, 0.500, 1.500, 0.000, 0.500, 1.000, 0.500, 
      0.500, 1.500, 0.500, 0.500, 1.500, 0.000, 0.000, 1.500, 0.500, 
      0.000, 1.500, 0.000, 0.000, 1.500, 0.500, 0.500, 1.500, 0.000, 
      
      0.000, 1.500, 0.000, 0.000, 1.750, 0.000, 0.250, 1.500, 0.000, 
      0.250, 1.750, 0.000, 0.250, 1.500, 0.000, 0.000, 1.750, 0.000, 
      0.000, 1.500, 0.000, 0.000, 1.500, 0.250, 0.000, 1.750, 0.000, 
      0.000, 1.750, 0.250, 0.000, 1.750, 0.000, 0.000, 1.500, 0.250, 
      0.000, 1.500, 0.000, 0.250, 1.500, 0.000, 0.000, 1.500, 0.250, 
      0.250, 1.500, 0.250, 0.000, 1.500, 0.250, 0.250, 1.500, 0.000, 
      0.250, 1.750, 0.250, 0.000, 1.750, 0.250, 0.250, 1.500, 0.250, 
      0.000, 1.500, 0.250, 0.250, 1.500, 0.250, 0.000, 1.750, 0.250, 
      0.250, 1.750, 0.250, 0.250, 1.500, 0.250, 0.250, 1.750, 0.000, 
      0.250, 1.500, 0.000, 0.250, 1.750, 0.000, 0.250, 1.500, 0.250, 
      0.250, 1.750, 0.250, 0.250, 1.750, 0.000, 0.000, 1.750, 0.250, 
      0.000, 1.750, 0.000, 0.000, 1.750, 0.250, 0.250, 1.750, 0.000};
    float diff = 0;
    for (int i=0;i<mesh.pos.size();i++)
      diff += abs(mesh.pos[i] - reference_model[i]);
    debug(" 16.2. %-64s", "Model is correct ");
    if (diff < 1e-6)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", diff, 1e-6);

    debug(" 16.3. %-64s", "Jacobian has correct structure ");
    //the jacobian is too large to actually make sure that it is correct
    //so we chech that it has the right structure:
    //[A 0 0 
    // 0 A 0
    // 0 0 A]
    int bx = mesh.pos.size()/3;
    int by = params.p.size()/3;
    for (int bi=0;bi<3;bi++)
    {
      for(int bj=0;bj<3;bj++)
      {
        if (bi != bj)
        {
          //it should be zero
          for (int i=0;i<by;i++)
            for (int j=0;j<bx;j++)
              if (abs(jac.at(bi*by + i, bj*bx + j)) > 1e-6)
              {
                logerr("block %d %d is not zero", bi,bj);
                goto fail;
              }
        }
        else
        {
          //it should be equal to zero block
          for (int i=0;i<by;i++)
            for (int j=0;j<bx;j++)
              if (abs(jac.at(bi*by + i, bj*bx + j) - jac.at(i, j)) > 1e-6)
              {
                logerr("blocks 0 0 and %d %d are not equal!", bi,bj);
                goto fail;
              }
        }
      }
    }
    if (true)
      debug("PASSED\n");
    else
fail: debug("FAILED\n");
  }

  //TEST 17 STACKED CUBES RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_17()
  {
    srand(0);
    debug("TEST 17. STACKED CUBES RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0,0,0, 1,1,1,   0,1,0, 0.5,0.5,0.5,  0,1.5,0, 0.25,0.25,0.25}
            structure:arr = {6,  6, 3,2,5, 3,2,5,  3,2,5}
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
            params:arr = {0.05,0.09,-0.07, 1.1,0.94,0.98,   0.04,0.96,0.08, 0.45,0.51,0.53,  -0.05,1.6,0.1, 0.23,0.26,0.28}
            structure:arr = {6,  6, 3,2,5, 3,2,5,  3,2,5}
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

    debug(" 17.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug(" 17.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug(" 17.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  void test_18()
  {
    debug("TEST 18. ROTATE JACOBIAN CALCULATION\n");

    upg::UPGStructure structure;
    structure.s = {4, 1};
    upg::UPGParametersRaw params;
    params.p = {1, 0, 0, 0,  0,0,0,0,0,1,0,1,0};
    upg::MeshGenInstance gen(structure);
    upg::UniversalGenJacobian jac;
    auto mesh = gen.generate(params.p, &jac);

    bool res;
    debug(" 18.1. %-64s", "Model and jacobian created, have right size ");
    if (mesh.pos.size() == 9 && jac.get_xn() == 9 && jac.get_yn() == 13)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %d %d %d\n", mesh.pos.size(), jac.get_xn(), jac.get_yn());
    }

    std::vector<float> reference_jac = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, -1, 0, 0, 0, 1,
                                        1, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 1, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 1, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 1, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 1, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 1, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 1, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 1, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 1};

    float diff = 0;
    for (int i=0;i<jac.get_yn();i++)
    {
      for (int j=0;j<jac.get_xn();j++)
      {
        diff += abs(reference_jac[i*jac.get_xn() + j] - jac.at(i,j));
      }
    }
    debug(" 18.2. %-64s", "Jacobian is correct ");
    if (diff < 1e-7)
    {
      debug("PASSED\n");
    }
    else
    {
      debug("FAILED %f > %f\n", diff, 1e-7);
    }
  }

  //TEST 19 STACKED CUBES RECONSTRUCTION
  //It uses Memetic+adam optimizer with no initial parameters set
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_19()
  {
    srand(0);
    debug("TEST 19. STACKED CUBES RECONSTRUCTION WITH MEMETIC\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            reference_image_w:i = 512
            reference_image_h:i = 512
            params:arr = {0,0,0, 0.5,0.5,0.5,   0,0.5,0, 0.4,0.4,0.4,  0,-0.5,0, 0.4,0.4,0.4}
            structure:arr = {6,  6, 3,2,5, 3,2,5,  3,2,5}
        } 
        view_0 {
            camera.origin:p3 = 3.000000, 0.500000, 3.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_1 {
            camera.origin:p3 = 3.000000, -0.500000, -3.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_2 {
            camera.origin:p3 = -3.000000, 0.500000, 3.000000
            camera.target:p3 = 0.000000, 0.000000, 0.000000
            camera.up:p3 = 0.000000, 1.000000, 0.000000
            camera.z_near:r = 0.100000
            camera.z_far:r = 100.000000
            camera.fov_rad:r = 1.00000
            camera.fixed:b = true
        }
        view_3 {
            camera.origin:p3 = -3.000000, -0.500000, -3.000000
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
            params:arr = {0.05,0.09,-0.07, 1.1,0.94,0.98,   0.04,0.96,0.08, 0.45,0.51,0.53,  -0.05,1.6,0.1, 0.23,0.26,0.28}
            structure:arr = {6,  6, 3,2,5, 3,2,5,  3,2,5}
        }
        step_0 {
            optimizer_name:s = "memetic"
            render_w:i = 256
            render_h:i = 256
            verbose:b = false
            save_intermediate_images:b = false
        }
        step_1 {
            optimizer_name:s = "adam"
            render_w:i = 512
            render_h:i = 512
            iterations:i = 500
            verbose:b = true
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

    debug(" 19.1. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);

    debug(" 19.2. %-64s", "Extremely high PSNR on given views ");
    if (res[0].quality_ir > 50)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_ir, 50);

    debug(" 19.3. %-64s", "Extremely high turntable PSNR ");
    if (res[0].quality_synt > 40)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 40);
  }

  //TEST 20 SDF NODES
  //tests the most basic SDF nodes functions
  //such as distance, merging and moving nodes
  //an putting derivative is right (root,left,right) order
  void test_20()
  {
    SdfGenInstance one_circle({std::vector<uint16_t>{1}});
    SdfGenInstance moved_circle({std::vector<uint16_t>{2,1}});
    SdfGenInstance two_circles({std::vector<uint16_t>{3,2,1,2,1}});

    debug("TEST 20. SDF NODES\n");
    {
    int pcnt_1 = one_circle.desc.get_total_params_count();
    int pcnt_2 = moved_circle.desc.get_total_params_count();
    int pcnt_3 = two_circles.desc.get_total_params_count();
    debug(" 20.1. %-64s", "SDF instances are created with expected number of parameters ");
    if (pcnt_1 == 4 && pcnt_2 == 7 && pcnt_3 == 14)
      debug("PASSED\n");
    else
      debug("FAILED %d %d %d\n", pcnt_1, pcnt_2, pcnt_3);
    }
    {
      std::vector<float> params = {0,1,0, 2};
      ProceduralSdf sdf = one_circle.generate(params);
      std::vector<float> ddist,dpos = {0,0,0};
      std::vector<float> ddist_ref = {0,1,0,-1}, dpos_ref = {0,-1,0};
      float d1 = sdf.get_distance({0,0,0},&ddist,&dpos);
      float d2 = sdf.get_distance({1,0,0});

      debug(" 20.2. %-64s", "Distance to circle correct");
      if (abs(d1 - (-1)) < 1e-6 && abs(d2 - (sqrtf(2)-2)) < 1e-6)
        debug("PASSED\n");
      else
        debug("FAILED %f %f\n", d1, d2);
      
      float dist_1=0,dist_2=0;
      for (int i=0;i<std::min(ddist.size(), ddist_ref.size());i++)
        dist_1 += std::abs(ddist[i]-ddist_ref[i]);
      for (int i=0;i<std::min(dpos.size(), dpos_ref.size());i++)
        dist_2 += std::abs(dpos[i]-dpos_ref[i]);
      
      debug(" 20.3. %-64s", "Derivatives correct");
      if (ddist.size() == 4 && dist_1 < 1e-6 && dpos.size() == 3 && dist_2 < 1e-6)
        debug("PASSED\n");
      else
        debug("FAILED %d %d %d %d\n", ddist.size() == 4, dist_1 < 1e-6, dpos.size() == 3, dist_2 < 1e-6);
    }
    {
      std::vector<float> params = {1,0,0 ,0,0,0,1.9,  -1,0,0 ,0,0,0,1};
      ProceduralSdf sdf = two_circles.generate(params);
      std::vector<float> ddist,dpos = {0,0,0};
      float d1 = sdf.get_distance({-1,1.1,0},&ddist,&dpos);

      debug(" 20.4. %-64s", "Distance to two triangles correct");
      if (abs(d1 - 0.1) < 1e-6)
        debug("PASSED\n");
      else
        debug("FAILED %f\n", d1);
      
      debug(" 20.5. %-64s", "Sdf merge derivatives correct");
      if (abs(ddist[0])+abs(ddist[1])+abs(ddist[2])+abs(ddist[3])+abs(ddist[4])+abs(ddist[5])+abs(ddist[6]) < 1e-6 &&
          abs(ddist[7])+abs(ddist[8])+abs(ddist[9])+abs(ddist[10])+abs(ddist[11])+abs(ddist[12])+abs(ddist[13]) > 1e-6)
        debug("PASSED\n");
      else
        debug("FAILED\n");
    }
  }

    //TEST 21 SPHERE SDF RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly
  void test_21()
  {
    srand(0);
    debug("TEST 21. SPHERE SDF RECONSTRUCTION\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 1000
            params:arr = {0,0,0,1}
            structure:arr = {1}
        } 
    }
    generator {

    }
    optimization {
        start {
            params:arr = {0.1,0.2,-0.1,0.7}    
            structure:arr = {1} 
        }
        step_0 {
            iterations:i = 100
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

    debug(" 21.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("PASSED\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    debug(" 21.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == 1 && res[0].structure.s[0] == 1)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug(" 21.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 4)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug(" 21.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug(" 21.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 80)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 80);
  }

  //TEST 22 SPHERES SDF RECONSTRUCTION MEMETIC
  void test_22()
  {
    srand(0);
    debug("TEST 22. SPHERES SDF RECONSTRUCTION MEMETIC\n");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            points_count:i = 1000
            params:arr = {0.6,0,0,0.5, -0.6,0,0,0.5, 0,-0.6,0,0.5, 0,0.6,0,0.5}
            structure:arr = {3,3,1,1,3,1,1}
        } 
    }
    generator {

    }
    optimization {
        start {
            structure:arr = {3,3,1,1,3,1,1}
        }
        step_0 {
            optimizer_name:s = "memetic"
            verbose:b = false
        }
        step_1 {
            iterations:i = 100
            verbose:b = true
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

    debug(" 22.1. %-64s", "ReconstructionResult size ");
    if (res.size() == 1)
      debug("PASSED\n");
    else
      debug("FAILED %d != %d\n", res.size(), 1);
    
    bool str_eq = true;
    std::vector<uint16_t> ref_struct = {3,3,1,1,3,1,1};
    for (int i=0;i<std::min(res[0].structure.s.size(), ref_struct.size());i++)
       str_eq = str_eq && (res[0].structure.s[i] == ref_struct[i]);
    debug(" 22.2. %-64s", "Preserved structure ");
    if (res[0].structure.s.size() == ref_struct.size() && str_eq)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug(" 22.3. %-64s", "Preserved parameters count ");
    if (res[0].parameters.p.size() == 4*4)
      debug("PASSED\n");
    else
      debug("FAILED\n");
    
    debug(" 22.4. %-64s", "Perfect optimization loss ");
    if (res[0].loss_optimizer < 1e-5)
      debug("PASSED\n");
    else
      debug("FAILED %f > %f\n", res[0].loss_optimizer, 1e-5);
    
    debug(" 22.5. %-64s", "Perfect multi-view PSNR ");
    if (res[0].quality_synt > 80)
      debug("PASSED\n");
    else
      debug("FAILED %f < %f\n", res[0].quality_synt, 80);
  }

  void perform_tests(const Block &blk)
  {

    Block *tests_blk = blk.get_block("tests");
    if (!tests_blk)
    {
      logerr("UPG Tests: tests block should exist in configuration");
      return;
    }

    std::vector<int> tests;
    tests_blk->get_arr("tests_num", tests);
    std::vector<std::function<void(void)>> test_functions = {
      test_1,  test_2,  test_3,  test_4,  test_5,
      test_6,  test_7,  test_8,  test_9,  test_10,
      test_11, test_12, test_13, test_14, test_15,
      test_16, test_17, test_18, test_19, test_20,
      test_21, test_22
    };

    for (int i : tests)
    {
      assert(i > 0 && i <= test_functions.size());
      test_functions[i-1]();
    }
  }
};