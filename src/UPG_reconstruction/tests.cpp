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
        start_parameters:arr = {0.1,0.1,0.1, -0.9,-0.1,-0.05, 0.07,0.85,-0.81}
        start_structure:arr = {1}
        render_w:i = 256
        render_h:i = 256
        iterations:i = 100
        verbose:b = false
        save_intermediate_images:b = false
    }
    results {
        check_image_quality:b = true
        check_model_quality:b = true
        save_folder:b = upg_triangle_reconstruction
        save_turntable:b = true
        save_reference_turntable:b = true
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