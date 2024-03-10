#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "myomp.h"

#include "common_utils/LiteMath_ext.h"

#include <cassert>
#include <iomanip>

#include "dmesh.h"
#include "functions.h"
#include "raytrace.h"

#include "optimizer.h"
#include "scenes.h"

#include "qmc.h"
#include "drender.h"

#include "utils.h"
#include "fin_diff.h"
#include "Image2d.h"
#include "tests.h"
//#include "enzyme.h"
#include <chrono>

constexpr static int  SAM_PER_PIXEL = 16;
double __enzyme_autodiff(void*, ...);

double square(double x) {
    return x * x;
}
double dsquare(double x) {
    // This returns the derivative of square or 2 * x
    return __enzyme_autodiff((void*) square, x);
}

using namespace diff_render;

#ifdef USE_CUSTOM_DIFF_RENDER
int custom_diff_render_main(int argc, char *argv[])
#else
int main(int argc, char *argv[]) //
#endif
{
  prepare_and_clear_directory("output");
  prepare_and_clear_directory("output/mitsuba_images");
  prepare_and_clear_directory("output/rendered");
  prepare_and_clear_directory("output/rendered_opt0");
  prepare_and_clear_directory("output/rendered_opt1");
  prepare_and_clear_directory("output/rendered_opt2");
  prepare_and_clear_directory("output/fin_diff");

  if (argc > 1 && ::std::string(argv[1]) == "-tests")
  {
    std::cout << "begin testing ... " << std::endl;
    Tester t;
    t.test_base_derivatives();
    
    #ifdef USE_MITSUBA
    std::cout << "mitsuba tests ... " << std::endl;
    t.test_3_1_mitsuba_triangle();
    t.test_3_2_mitsuba_sphere();
    t.test_3_3_mitsuba_teapot();
    t.test_3_4_mitsuba_cube();
    #endif

    std::cout << "2.1--2.4 tests ... " << std::endl;
    t.test_2_1_triangle();
    t.test_2_2_pyramid();
    t.test_2_3_sphere();
    t.test_2_4_pyramid_vcol();

    std::cout << "2.5--2.11 tests ... " << std::endl;
    t.test_2_5_teapot_diffuse(); // Embree crashes
    t.test_2_7_mesh_on_static_scene();
    t.test_2_8_instancing();
    t.test_2_9_transform();
    t.test_2_10_multiple_meshes();
    t.test_2_11_restricted_transforms();
    
    std::cout << "end testing ... " << std::endl;
    return 0;
  }

  Img img(256, 256);
  
  constexpr int camsNum = 3;
  CamInfo cameras[camsNum] = {};
  for(int i=0;i<camsNum;i++) {
    cameras[i].width  = float(img.width());
    cameras[i].height = float(img.height());
    cameras[i].mWorldView.identity();
    cameras[i].mProj.identity();
  }

  float4x4 mProj = LiteMath::perspectiveMatrix(45.0f, cameras[0].width / cameras[0].height, 0.1f, 100.0f);

  cameras[0].mProj      = mProj;
  cameras[0].mWorldView = LiteMath::translate4x4(float3(0,0,-3));

  cameras[1].mProj      = mProj;
  cameras[1].mWorldView = LiteMath::translate4x4(float3(0,0,-3))*LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*120.0f)*LiteMath::rotate4x4X(LiteMath::DEG_TO_RAD*45.0f);

  cameras[2].mProj      = mProj;
  cameras[2].mWorldView = LiteMath::translate4x4(float3(0,0,-3))*LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*(-120.0f))*LiteMath::rotate4x4X(LiteMath::DEG_TO_RAD*(-45.0f));

  for(int i=0;i<camsNum;i++)
    cameras[i].commit();

  auto g_uniforms = cameras[0];

  Scene initialScene, targetScene;
  SHADING_MODEL mode = SHADING_MODEL::LAMBERT;
  {
    TriangleMesh initialMesh, targetMesh;
    TriangleMesh initialMesh2, targetMesh2;
    //scn01_TwoTrisFlat(initialMesh, targetMesh);
    //scn02_TwoTrisSmooth(initialMesh, targetMesh);
    //scn03_Triangle3D_White(initialMesh, targetMesh);
    //scn04_Triangle3D_Colored(initialMesh, targetMesh); // bad
    //scn05_Pyramid3D(initialMesh, targetMesh);
    //scn06_Cube3D_VColor(initialMesh, targetMesh);      // bad     
    //scn08_Cube3D_Textured(initialMesh, targetMesh);
    scn09_Sphere3D_Textured(initialMesh, targetMesh);
    initialScene.add_mesh(initialMesh, {LiteMath::translate4x4(float3(0,-0.1,0))});
    initialScene.add_mesh(initialMesh, {LiteMath::translate4x4(float3(0.5,0.5,0)), LiteMath::translate4x4(float3(-0.5,0.5,0))});
    
    targetScene.add_mesh(targetMesh, {LiteMath::float4x4()});
    targetScene.add_mesh(initialMesh, {LiteMath::translate4x4(float3(0.5,0.5,0)), LiteMath::translate4x4(float3(-0.5,0.5,0))});
  }

  auto pDRender = MakeDifferentialRenderer({mode, SAM_PER_PIXEL});
  
  Img targets[camsNum];
  for(int i=0;i<camsNum;i++) {
    targets[i].resize(img.width(), img.height());
    targets[i].clear(float3{0,0,0});
  }

  pDRender->commit(targetScene);
  pDRender->render(targetScene, cameras, targets, camsNum);

  for(int i=0;i<camsNum;i++) {
    ::std::stringstream strOut;
    strOut  << "output/rendered_opt" << i << "/z_target.bmp";
    auto temp = strOut.str();
    LiteImage::SaveImage(temp.c_str(), targets[i]);
  }

  IOptimizer* pOpt = CreateSimpleOptimizer();

  OptimizerParameters op = OptimizerParameters(OptimizerParameters::GD_Adam);
  op.position_lr = 0.0;
  op.textures_lr = 0.2;
  op.transforms_lr = 0.1;
  pOpt->Init(initialScene, pDRender, cameras, targets, 3, op);

  float error = 0;
  int iters = 300;
  ::std::chrono::steady_clock::time_point t1 = ::std::chrono::steady_clock::now();
  Scene res_scene = pOpt->Run(iters, error);
  ::std::chrono::steady_clock::time_point t2 = ::std::chrono::steady_clock::now();
  float r = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t2 - t1).count();
  logerr("optimization took %.1f ms per iteration\n", r/iters);
  
  //img.clear(float3{0,0,0});
  //pDRender->commit(mesh3);
  //pDRender->render(mesh3, cameras, &img, 1);
  //LiteImage::SaveImage("output/rendered_opt/z_target2.bmp", img);
  
  delete pOpt; pOpt = nullptr;
  return 0;
}
