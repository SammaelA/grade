#include "upg.h"

namespace upg
{
  //TEST 1 ONE TRIANGLE RECONSTRUCTION
  //It uses Adam optimizer with initial state close to target one
  //Reconstruction should perform perfectly (like 90 PSNR)
  void test_1()
  {
    debug("TEST 1. ONE TRIANGLE RECONSTRUCTION: ");
    std::string settings = R""""(
    {
    input {
        synthetic_reference {
            tex_w:i = 256
            tex_h:i = 256
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
    if (res.empty())
      debug("FAILED. No reconstruction result returned");
    else if (res[0].quality_ir < 80)
      debug("FAILED. Low PSNR (%.1f < 80.0(target))", res[0].quality_ir);
    else
      debug("PASSED.");
    debugnl();
  }


  void perform_tests()
  {
    test_1();
  }
};