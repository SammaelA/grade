#include "tests.h"
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "myomp.h"
#include <functional>
#include <chrono>

#include "common_utils/LiteMath_ext.h"
using namespace LiteMath;

#include <cassert>
#include <iomanip>

#include "dmesh.h"
#include "drender.h"
#include "scenes.h"
#include "optimizer.h"
#include "drender_mitsuba.h"
namespace diff_render
{
#define TESTER_PSNR_PASS 35
#define TESTER_FIN_DIFF_PASS 0.5
#define TESTER_FIN_DIFF_TEXTURE_PASS 0.01
#define TESTER_MITSUBA_PASS 0.1

void Tester::test_base_derivatives()
{

  constexpr int IMAGE_W = 256;
  constexpr int IMAGE_H = 256;
  constexpr int SILHOUETTE_SPP = 16;
  constexpr int BASE_SPP = 16;
  
  CamInfo camera;
  camera.width  = float(IMAGE_W);
  camera.height = float(IMAGE_H);
  camera.mWorldView = LiteMath::translate4x4(float3(0,0,-3));
  camera.mProj = LiteMath::perspectiveMatrix(45.0f, camera.width / camera.height, 0.1f, 100.0f);
  camera.commit();

  {
    Scene initialScene, targetScene;
    TriangleMesh initialMesh, targetMesh;
    scn03_Triangle3D_White(initialMesh, targetMesh);
    initialScene.add_mesh(initialMesh);
    targetScene.add_mesh(targetMesh);
    auto res = test_derivatives(initialScene, targetScene, camera, {SHADING_MODEL::SILHOUETTE, SILHOUETTE_SPP}, 100, 0);

    bool pass = res.pos_error < TESTER_FIN_DIFF_PASS;
    printf("%s TEST 1.1: EDGE SAMPLING TRIANGLE with error %.3f\n", pass ? "    PASSED:" : "FAILED:    ", res.pos_error);
  }

  {
    Scene initialScene, targetScene;
    TriangleMesh initialMesh, targetMesh;
    scn05_Pyramid3D(initialMesh, targetMesh);
    initialScene.add_mesh(initialMesh);
    targetScene.add_mesh(targetMesh);
    auto res = test_derivatives(initialScene, targetScene, camera, {SHADING_MODEL::SILHOUETTE, SILHOUETTE_SPP}, 100, 0);

    bool pass = res.pos_error < TESTER_FIN_DIFF_PASS;
    printf("%s TEST 1.2: EDGE SAMPLING PYRAMID with error %.3f\n", pass ? "    PASSED:" : "FAILED:    ", res.pos_error);
  }

  {
    Scene initialScene, targetScene;
    TriangleMesh initialMesh, targetMesh;
    scn09_Sphere3D_Textured(initialMesh, targetMesh);
    LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.2f,0.0f));
    transform(initialMesh, mTranslate);
    initialScene.add_mesh(initialMesh); initialScene.restore_meshes(false, false, true);
    targetScene.add_mesh(targetMesh); targetScene.restore_meshes(false, false, true);
    auto res = test_derivatives(initialScene, targetScene, camera, {SHADING_MODEL::SILHOUETTE, SILHOUETTE_SPP}, 100, 0);

    bool pass = res.pos_error < TESTER_FIN_DIFF_PASS;
    printf("%s TEST 1.3: EDGE SAMPLING SPHERE with error %.3f\n", pass ? "    PASSED:" : "FAILED:    ", res.pos_error);
  }

  {
    Scene initialScene, targetScene;
    TriangleMesh initialMesh, targetMesh;
    scn05_Pyramid3D(initialMesh, targetMesh);
    initialScene.add_mesh(initialMesh);
    targetScene.add_mesh(targetMesh);
    auto res = test_derivatives(initialScene, targetScene, camera, {SHADING_MODEL::VERTEX_COLOR, SILHOUETTE_SPP}, 100, 0);

    bool pass = res.color_error < TESTER_FIN_DIFF_PASS;
    printf("%s TEST 1.4: VCOLOR DERIVATIVES with error %.3f\n", pass ? "    PASSED:" : "FAILED:    ", res.color_error);
  }

  {
    Scene initialScene, targetScene;
    TriangleMesh initialMesh, targetMesh;
    scn09_Sphere3D_Textured(initialMesh, targetMesh);
    initialScene.add_mesh(initialMesh);
    targetScene.add_mesh(targetMesh);
    auto res = test_derivatives(initialScene, targetScene, camera, {SHADING_MODEL::TEXTURE_COLOR, SILHOUETTE_SPP}, 0, 250);

    bool pass = res.texture_error < TESTER_FIN_DIFF_TEXTURE_PASS;
    printf("%s TEST 1.5: TEXTURE DERIVATIVES with error %.3f\n", pass ? "    PASSED:" : "FAILED:    ", res.texture_error);
  }
}

::std::vector<CamInfo> create_cameras_around(int cam_num, int sensor_w, int sensor_h)
{
  ::std::vector<CamInfo> cameras;
  float4x4 mProj = LiteMath::perspectiveMatrix(45.0f, sensor_w/sensor_h, 0.1f, 100.0f);
  float rot_y = (2*M_PI)/cam_num;
  float rot_x = (M_PI_2)/cam_num;
  int st = cam_num/2;
  for(int i=0;i<cam_num;i++)
  {
    auto rot = LiteMath::rotate4x4Y(rot_y*(i-st))*LiteMath::rotate4x4Y(rot_x*(i-st));
    cameras.push_back(CamInfo(float3(-5*sin(rot_y*(i-st)),0,-5*cos(rot_y*(i-st))), float3(0,0,0), float3(0,1,0), sensor_w, sensor_h));
  }

  return cameras;
}

void optimization_test(const ::std::string &test_name,
                       ::std::function<void(Scene&, Scene&)> create_scene,
                       const DiffRenderSettings &diff_render_settings,
                       const OptimizerParameters &opt_parameters,
                       int opt_steps = 300,
                       bool test_by_steps = false,
                       int cameras_count = 3,
                       int image_w = 256,
                       int image_h = 256,
                       ::std::function<::std::vector<CamInfo>(int, int, int)> create_cameras = create_cameras_around)
{
  ::std::chrono::steady_clock::time_point t1 = ::std::chrono::steady_clock::now();

  auto cameras = create_cameras(cameras_count, image_w, image_h);
  Scene initialScene, targetScene;
  create_scene(initialScene, targetScene);

  ::std::chrono::steady_clock::time_point t2 = ::std::chrono::steady_clock::now();

  auto pDRender = MakeDifferentialRenderer(diff_render_settings);

  ::std::chrono::steady_clock::time_point t3 = ::std::chrono::steady_clock::now();

  Img img(image_w, image_h);
  Img targets[cameras_count];
  for (int i = 0; i < cameras_count; i++)
  {
    targets[i].resize(img.width(), img.height());
    targets[i].clear(float3{0, 0, 0});
  }

  pDRender->commit(targetScene);
  pDRender->render(targetScene, cameras.data(), targets, cameras_count);

  ::std::chrono::steady_clock::time_point t4 = ::std::chrono::steady_clock::now();

  for(int i=0;i<cameras_count;i++) 
  {
    ::std::stringstream strOut;
    strOut  << "output/rendered_opt" << i << "/z_target.bmp";
    auto temp = strOut.str();
    LiteImage::SaveImage(temp.c_str(), targets[i]);
  }
  ::std::chrono::steady_clock::time_point t5 = ::std::chrono::steady_clock::now();

  ::std::chrono::steady_clock::time_point t6 = t5, t7 = t5;
  IOptimizer *pOpt = CreateSimpleOptimizer();
  float error = 1e9;
  if (!test_by_steps)
  {
    pOpt->Init(initialScene, pDRender, cameras.data(), targets, cameras_count, opt_parameters);
    t6 = ::std::chrono::steady_clock::now();
    Scene res_scene = pOpt->Run(opt_steps, error);
    t7 = ::std::chrono::steady_clock::now();
  }
  else
  {
    ::std::vector<Scene> iter_scenes;
    pOpt->Init(initialScene, pDRender, cameras.data(), targets, cameras_count, opt_parameters);
    Scene res_scene = pOpt->Run(opt_steps, error, &iter_scenes);
    for (auto &s : iter_scenes)
    {
      auto res = Tester::test_derivatives(s, targetScene, cameras[2], diff_render_settings, 100, 0);
      printf("derivatives error %f %f\n",res.pos_error, res.texture_error);      
    }
  }

  float psnr = -10*log10(max(1e-9f,error));
  bool pass = psnr > TESTER_PSNR_PASS;
  printf("%s %s with PSNR %.3f ", pass ? "    PASSED:" : "FAILED:    ", test_name.c_str(), psnr);

  float scene_ms = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t2 - t1).count();
  float dr_ms = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t3 - t2).count();
  float targ1_ms = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t4 - t3).count();
  float targ2_ms = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t5 - t4).count();
  float opt_ms = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t6 - t5).count();
  float main_ms = 1e-3 * ::std::chrono::duration_cast<::std::chrono::microseconds>(t7 - t6).count();
  printf(" %.1f + %.1f ms/iter [%.1f %.1f %.1f %.1f %.1f %.1f]\n", main_ms/opt_steps, (scene_ms+dr_ms+opt_ms)/opt_steps,
         scene_ms, dr_ms, targ1_ms, targ2_ms, opt_ms, main_ms);
}

void optimization_test(const ::std::string &test_name,
                       ::std::function<void(TriangleMesh&, TriangleMesh&)> create_mesh,
                       const DiffRenderSettings &diff_render_settings,
                       const OptimizerParameters &opt_parameters,
                       int opt_steps = 300,
                       bool test_by_steps = false,
                       int cameras_count = 3,
                       int image_w = 256,
                       int image_h = 256,
                       ::std::function<::std::vector<CamInfo>(int, int, int)> create_cameras = create_cameras_around)
{
  optimization_test(test_name, 
                    [&](Scene &initialScene, Scene &targetScene){
                      TriangleMesh initialMesh, targetMesh;
                      create_mesh(initialMesh, targetMesh);
                      initialScene.add_mesh(initialMesh);
                      targetScene.add_mesh(targetMesh);
                    },
                    diff_render_settings, opt_parameters, opt_steps, test_by_steps, cameras_count, image_w, image_h, create_cameras);

}

void Tester::test_2_1_triangle()
{
  optimization_test("TEST 2.1: ONE TRIANGLE SHAPE OPTIMIZATION",
                    scn03_Triangle3D_White,
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.05, 0.1},
                    300);
}

void Tester::test_2_2_pyramid()
{
  optimization_test("TEST 2.2: PYRAMID SHAPE OPTIMIZATION",
                    scn05_Pyramid3D,
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.05, 0.1},
                    300);
}

void Tester::test_2_3_sphere()
{
  optimization_test("TEST 2.3: SPHERE SHAPE OPTIMIZATION",
                    scn09_Sphere3D_Textured,
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.05, 0.1},
                    300);
}

void Tester::test_2_4_pyramid_vcol()
{
  optimization_test("TEST 2.4: PYRAMID VCOL+POS OPTIMIZATION",
                    scn05_Pyramid3D,
                    {SHADING_MODEL::VERTEX_COLOR, 16},
                    {OptimizerParameters::GD_Adam, 0.05, 0.0},
                    300);
}

void Tester::test_2_5_teapot_diffuse()
{
  optimization_test("TEST 2.5: TEAPOT TEXTURE_COLOR TEXTURE OPTIMIZATION",
                    scn10_Teapot3D_Textured,
                    {SHADING_MODEL::TEXTURE_COLOR, 16},
                    {OptimizerParameters::GD_Adam, 0.0, 0.05},
                    300);
}

void Tester::test_2_7_mesh_on_static_scene()
{
  optimization_test("TEST 2.7: SHAPE OPTIMIZATION ON STATIC SCENE",
                    [&](Scene &initialScene, Scene &targetScene){
                      TriangleMesh initialMesh, targetMesh;
                      scn05_Pyramid3D(initialMesh, targetMesh);
                      TriangleMesh _mm1, _mm2;
                      scn09_Sphere3D_Textured(_mm1, _mm2);
                      initialScene.add_mesh(initialMesh);
                      initialScene.add_mesh(_mm1, {LiteMath::translate4x4(float3(0.5,0.5,0)), LiteMath::translate4x4(float3(-0.5,0.5,0))});

                      targetScene.add_mesh(targetMesh);
                      targetScene.add_mesh(_mm1, {LiteMath::translate4x4(float3(0.5,0.5,0)), LiteMath::translate4x4(float3(-0.5,0.5,0))});
                    },
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.00, 0.0, 0.1},
                    300);
}

void Tester::test_2_8_instancing()
{
  optimization_test("TEST 2.8: SHAPE OPTIMIZATION WITH INSTANCING",
                    [&](Scene &initialScene, Scene &targetScene){
                      TriangleMesh initialMesh, targetMesh;
                      scn05_Pyramid3D(initialMesh, targetMesh);
                      ::std::vector<float4x4> initial_transforms;
                      ::std::vector<float4x4>  target_transforms;
                      int tr = 2;
                      float rot = (2*M_PI)/tr;
                      for (int i=0;i<tr;i++)
                      {
                        initial_transforms.push_back(LiteMath::scale4x4({0.5,0.5,0.5})*LiteMath::translate4x4({2.5f*(1-i),0.5f,0.0f}));
                        target_transforms.push_back(LiteMath::scale4x4({0.5,0.5,0.5})*LiteMath::translate4x4( {2.5f*(1-i),0.0f,0.0f}));
                        //initial_transforms.push_back(LiteMath::scale4x4({0.5,0.5,0.5})*LiteMath::rotate4x4Z(i*rot)*LiteMath::translate4x4({0,1,0}));
                        //target_transforms.push_back(LiteMath::scale4x4({0.5,0.5,0.5})*LiteMath::rotate4x4Z((i+0.5)*rot)*LiteMath::translate4x4({0,1,0}));
                      }
                      initialScene.add_mesh(targetMesh, initial_transforms);
                      targetScene.add_mesh(targetMesh, target_transforms);
                    },
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.00, 0.0, 0.1},
                    300);
}

void Tester::test_2_9_transform()
{
  optimization_test("TEST 2.9: TEAPOT TRANSFORM MATRIX OPTIMIZATION",
                    [&](Scene &initialScene, Scene &targetScene){
                      TriangleMesh initialMesh, targetMesh;
                      scn05_Pyramid3D(initialMesh, targetMesh);
                      initialScene.add_mesh(initialMesh, {LiteMath::rotate4x4Z(0.2)*LiteMath::rotate4x4Y(0.2)});
                      targetScene.add_mesh(initialMesh);
                    },
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.0, 0.0, 0.1},
                    300);
}

void Tester::test_2_10_multiple_meshes()
{
  optimization_test("TEST 2.10: SHAPE OPTIMIZATION MULTIPLE MESHES",
                    [&](Scene &initialScene, Scene &targetScene){
                      TriangleMesh initialMesh, targetMesh;
                      scn05_Pyramid3D(initialMesh, targetMesh);
                      TriangleMesh _mm1, _mm2;
                      scn03_Triangle3D_White(_mm1, _mm2);
                      initialScene.add_mesh(initialMesh, {LiteMath::translate4x4(float3(0.0,1,0))});
                      initialScene.add_mesh(_mm1, {LiteMath::translate4x4(float3(0.0,-1,0))});

                      targetScene.add_mesh(targetMesh, {LiteMath::translate4x4(float3(0.0,1,0))});
                      targetScene.add_mesh(_mm2, {LiteMath::translate4x4(float3(0.0,-1,0))});
                    },
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.0, 0.0, 0.05, {0,1}},
                    300);
}

void Tester::test_2_11_restricted_transforms()
{
  optimization_test("TEST 2.11: RESTRICTED TRANSFORMS OPTIMIZATION",
                    [&](Scene &initialScene, Scene &targetScene){
                      TriangleMesh initialMesh, targetMesh;
                      scn05_Pyramid3D(initialMesh, targetMesh);
                      initialScene.add_mesh(initialMesh, {TransformR({0,0,0}, {0.0,0.2,0.2}, 1)});
                      targetScene.add_mesh(initialMesh, {TransformR()});
                    },
                    {SHADING_MODEL::SILHOUETTE, 16},
                    {OptimizerParameters::GD_Adam, 0.0, 0.0, 0.05},
                    300);
}

void mitsuba_compare_test(const ::std::string &test_name,
                          ::std::function<void(TriangleMesh&, TriangleMesh&)> create_scene,
                          const DiffRenderSettings &diff_render_settings,
                          int cameras_count = 3,
                          int image_w = 256,
                          int image_h = 256,
                          ::std::function<::std::vector<CamInfo>(int, int, int)> create_cameras = create_cameras_around)
{
  auto cameras = create_cameras(cameras_count, image_w, image_h);

  Scene initialScene, targetScene;
  TriangleMesh initialMesh, targetMesh;
  create_scene(initialMesh, targetMesh);
  initialScene.add_mesh(initialMesh);
  targetScene.add_mesh(targetMesh);
  initialScene.restore_meshes(true, true, true);
  targetScene.restore_meshes(true, true, true);

  auto pDRender = new DiffRenderMitsuba();
  pDRender->init(diff_render_settings);
  auto pDRenderOurs = MakeDifferentialRenderer(diff_render_settings);

  Img img(image_w, image_h);
  Img targets[cameras_count];
  Img targetsOurs[cameras_count];
  for (int i = 0; i < cameras_count; i++)
  {
    targets[i].resize(img.width(), img.height());
    targets[i].clear(float3{0, 0, 0});
    targetsOurs[i].resize(img.width(), img.height());
    targetsOurs[i].clear(float3{0, 0, 0});
  }

  DScene gradMeshMitsuba;
  {
    pDRender->commit(targetScene);
    pDRender->render(targetScene, cameras.data(), targets, cameras_count);
    gradMeshMitsuba.reset(initialScene, pDRender->mode, {0});
    pDRender->d_render_and_compare(initialScene, cameras.data(), targets, cameras_count, img.width()*img.height(), gradMeshMitsuba);
  }

  DScene gradMeshOurs;
  {
    pDRenderOurs->commit(targetScene);
    pDRenderOurs->render(targetScene, cameras.data(), targetsOurs, cameras_count);
    gradMeshOurs.reset(initialScene, pDRenderOurs->mode, {0});
    pDRenderOurs->d_render_and_compare(initialScene, cameras.data(), targetsOurs, cameras_count, img.width()*img.height(), gradMeshOurs);
  }

  assert(gradMeshMitsuba.get_dmeshes().size() == gradMeshOurs.get_dmeshes().size());
  float diff = 0;
  int full_sz = 0;
  for (int m=0; m<gradMeshMitsuba.get_dmeshes().size(); m++)
  {
    int sz = gradMeshOurs.get_dmeshes()[m].full_size();
    full_sz += sz;
    double acc1 = 1e-12;
    double acc2 = 1e-12;
    for (int i=0;i<sz;i++)
      acc1 += gradMeshMitsuba.get_dmeshes()[m].full_data()[i];
    for (int i=0;i<sz;i++)
      acc2 += gradMeshOurs.get_dmeshes()[m].full_data()[i];
    for (int i=0;i<sz;i++)
    {
      //logerr("[%d]%f %f",i, sz*gradMeshMitsuba[i]/acc1, sz*gradMeshOurs[i]/acc2);
      diff += abs(gradMeshMitsuba.get_dmeshes()[m].full_data()[i]/acc1 - gradMeshOurs.get_dmeshes()[m].full_data()[i]/acc2);
    }
  }

  double image_diff = 0.0;
  for (int i=0;i<cameras_count;i++)
  {
    ::std::string path1 = "output/t"+::std::to_string(i)+"_1.png";
    ::std::string path2 = "output/t"+::std::to_string(i)+"_2.png";
    LiteImage::SaveImage(path1.c_str(), targets[i]);
    LiteImage::SaveImage(path2.c_str(), targetsOurs[i]);
    image_diff += LossAndDiffLoss(targets[i], targetsOurs[i], img);
  }
  image_diff /= cameras_count*img.width()*img.height();
  float psnr = -10*log10(max(1e-9,image_diff));
  diff /= full_sz;
  delete pDRender;
    
  printf("%s %s with image PSNR %.5f\n", (psnr > TESTER_PSNR_PASS) ? "    PASSED:" : "FAILED:    ", test_name.c_str(), psnr);
  printf("%s %s with derivatives difference %.5f\n", (diff < TESTER_MITSUBA_PASS) ? "    PASSED:" : "FAILED:    ", test_name.c_str(), diff);
}

void Tester::test_3_1_mitsuba_triangle()
{
  mitsuba_compare_test("TEST 3.1: TRIANGLE MITSUBA COMPARE",
                       scn03_Triangle3D_White,
                       {SHADING_MODEL::SILHOUETTE, 16},
                       1,
                       512,
                       512);
}

void Tester::test_3_2_mitsuba_sphere()
{
  mitsuba_compare_test("TEST 3.2: SPHERE MITSUBA COMPARE",
                       scn09_Sphere3D_Textured,
                       {SHADING_MODEL::SILHOUETTE, 16},
                       1,
                       512,
                       512);
}

void Tester::test_3_3_mitsuba_teapot()
{
  mitsuba_compare_test("TEST 3.3: TEAPOT MITSUBA COMPARE",
                       scn11_Teapot3D_Textured,
                       {SHADING_MODEL::SILHOUETTE, 16},
                       1,
                       512,
                       512);
}

void Tester::test_3_4_mitsuba_cube()
{
  mitsuba_compare_test("TEST 3.4: CUBE MITSUBA COMPARE",
                       scn06_Cube3D_VColor,
                       {SHADING_MODEL::SILHOUETTE, 16},
                       3,
                       512,
                       512);
}

void finDiff_param(float *param, GradReal *param_diff, Img &out_diffImage, float delta,
                   Img &img, const Img &target, const Img &MSEOrigin, const CamInfo& a_camData, const Scene &copy_scene,
                   ::std::shared_ptr<IDiffRender> a_pDRImpl, bool geom_changed)
{
  *param += delta;
  copy_scene.invalidate_prepared_scene();

  if (geom_changed)
    a_pDRImpl->commit(copy_scene);
  a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
  out_diffImage = (LiteImage::MSEImage(img,target) - MSEOrigin)/delta;   
  float3 summColor = SummOfPixels(out_diffImage); 
  *param_diff += GradReal(summColor.x + summColor.y + summColor.z);

  *param -= delta;
}

//cnt unique integers in range [from, to)
::std::vector<int> random_unique_indices(int from, int to, int cnt)
{
  if (cnt <= 0)
    return {};
  ::std::vector<int> res;
  if (to - from <= cnt)
  {
    res.resize(to-from, 0);
    for (int i=from;i<to;i++)
      res[i-from] = i;
  }
  else
  {
    int sh_cnt = min(4*cnt, to - from);
    int step = (to - from)/sh_cnt;
    ::std::vector<int> sh;
    sh.reserve((to - from)/step);
    if (step == 1)
    {
      for (int i=from;i<to;i+=step)
        sh.push_back(i);
    }
    else
    {
      for (int i=from;i<to;i+=step)
        sh.push_back(i + rand()%step);
    }

    ::std::shuffle(sh.begin(), sh.end(), std::default_random_engine{});
    res = ::std::vector<int>(sh.begin(), sh.begin() + cnt);
  }
  return res;
}

void Tester::test_fin_diff(const Scene &scene, const char* outFolder, const Img& origin, const Img& target, ::std::shared_ptr<IDiffRender> a_pDRImpl,
                           const CamInfo& a_camData, DScene &d_scene, int debug_mesh_id, int max_test_vertices, int max_test_texels,
                           ::std::vector<bool> &tested_mask)
{
  float dPos = 2.0f/float(origin.width());
  float dCol = 0.1f;
  bool save_images = true;

  Scene copy_scene = scene;
  a_pDRImpl->commit(copy_scene);
  TriangleMesh &mesh = copy_scene.get_mesh_modify(debug_mesh_id);
  d_scene.reset(copy_scene, a_pDRImpl->mode, {debug_mesh_id});
  DMesh &d_mesh = *d_scene.get_dmesh(debug_mesh_id);
  Img MSEOrigin = LiteImage::MSEImage(origin, target);
  Img img(origin.width(), origin.height());

  #define pFinDiff(param, dmesh_ptr, out_image, delta, geom_changed) \
          tested_mask[(int)((dmesh_ptr) - d_mesh.full_data())] = true; finDiff_param(param, dmesh_ptr, out_image, delta, img, target, MSEOrigin, a_camData, copy_scene, a_pDRImpl, geom_changed);
  
  ::std::vector<int> debug_vertex_ids = random_unique_indices(0, mesh.vertex_count(), max_test_vertices);
  Img pos_x, pos_y, pos_z;
  for(auto &i : debug_vertex_ids)
  {
    pFinDiff(mesh.vertices[i].M + 0, d_mesh.pos(i) + 0, pos_x, dPos, true);
    pFinDiff(mesh.vertices[i].M + 1, d_mesh.pos(i) + 1, pos_y, dPos, true);
    pFinDiff(mesh.vertices[i].M + 2, d_mesh.pos(i) + 2, pos_z, dPos, true);

    if (save_images)
    {
      Img diffImage(pos_x.width(), pos_x.height()); 
      for(int y=0;y<pos_x.height();y++)
        for(int x=0;x<pos_x.width();x++)
          diffImage[int2(x,y)] = float3(pos_x[int2(x,y)].x, pos_y[int2(x,y)].x, pos_z[int2(x,y)].x);

      if(outFolder != nullptr)
      {
        ::std::stringstream strOut;
        strOut << outFolder << "/" << "pos_xyz_" << i << ".bmp";
        auto path = strOut.str();
        LiteImage::SaveImage(path.c_str(), diffImage);
      }
      if(outFolder != nullptr)
      {
        ::std::stringstream strOut;
        strOut << outFolder << "/" << "orig_xyz_" << i << ".bmp";
        auto path = strOut.str();
        LiteImage::SaveImage(path.c_str(), img);
      }
    }

    if (mesh.colors.size() == mesh.vertices.size() && a_pDRImpl->mode == SHADING_MODEL::VERTEX_COLOR)
    {
      pFinDiff(mesh.colors[i].M + 0, d_mesh.color(i) + 0, pos_x, dCol, true);
      pFinDiff(mesh.colors[i].M + 1, d_mesh.color(i) + 1, pos_y, dCol, true);
      pFinDiff(mesh.colors[i].M + 2, d_mesh.color(i) + 2, pos_z, dCol, true);

      if (save_images)
      {
        Img diffImage(pos_x.width(), pos_x.height()); 
        for(int y=0;y<pos_x.height();y++)
          for(int x=0;x<pos_x.width();x++)
            diffImage[int2(x,y)] = float3(pos_x[int2(x,y)].x, pos_y[int2(x,y)].x, pos_z[int2(x,y)].x);

        if(outFolder != nullptr)
        {
          ::std::stringstream strOut;
          strOut << outFolder << "/" << "color_xyz_" << i << ".bmp";
          auto path = strOut.str();
          LiteImage::SaveImage(path.c_str(), diffImage);
        }
      }
    }
  }

  for (int tex_n = 0; tex_n < mesh.textures.size(); tex_n++)
  {
    CPUTexture &tex = mesh.textures[tex_n];
    ::std::vector<int2> debug_texel_ids;
    ::std::vector<int> dtids = random_unique_indices(0, tex.w*tex.h, max_test_texels);
    for(auto ind : dtids)
    {
      debug_texel_ids.push_back(int2(ind % tex.w, ind / tex.w));
    }
    
    for(auto &i : debug_texel_ids)
    {
      for (int ch =0; ch < tex.channels; ch++)
      {
        int off = tex.pixel_to_offset(i.x, i.y) + ch;
        pFinDiff(tex.data.data() + tex.pixel_to_offset(i.x, i.y) + ch, d_mesh.tex(tex_n, off), pos_x, 0.1, false);
      }
    }
  }
}

Tester::DerivativesTestResults Tester::test_derivatives(const Scene &initial_scene, const Scene &target_scene, const CamInfo& a_camData, 
                                                        const DiffRenderSettings &settings, int max_test_vertices, int max_test_texels, bool print)
{
  Img original = Img(a_camData.width, a_camData.height);
  Img target = Img(a_camData.width, a_camData.height);
  Img tmp = Img(a_camData.width, a_camData.height);

  auto Render = MakeDifferentialRenderer(settings);

  Render->commit(target_scene);
  Render->render(target_scene, &a_camData, &target, 1);

  Render->commit(initial_scene);
  Render->render(initial_scene, &a_camData, &original, 1);

  LossAndDiffLoss(original, target, tmp); 
  DerivativesTestResults r;
  int meshes_count = initial_scene.get_meshes().size();
  for (int i=0;i<meshes_count;i++)
  {
    DScene dMesh1 = DScene(initial_scene, settings.mode, {i});
    DScene dMesh2 = DScene(initial_scene, settings.mode, {i});
    ::std::vector<bool> mask(dMesh1.get_dmeshes()[0].full_size(), false);
    
    Render->commit(initial_scene);
    Render->d_render(initial_scene, &a_camData, &tmp, 1, target.width()*target.height(), dMesh1);
    test_fin_diff(initial_scene, "output", original, target, Render, a_camData, dMesh2, i, max_test_vertices, max_test_texels, mask);

    auto rm = PrintAndCompareGradients(dMesh1, dMesh2, mask, print);

    r.pos_error += rm.pos_error/meshes_count;
    r.color_error += rm.color_error/meshes_count;
    r.texture_error += rm.texture_error/meshes_count;
    r.average_error += rm.average_error/meshes_count;
  }

  return r;
}

Tester::DerivativesTestResults Tester::PrintAndCompareGradients(DScene& grad1_scene, DScene& grad2_scene, ::std::vector<bool> &tested_mask, bool print)
{
  double posError = 0.0;
  double colError = 0.0;
  double texError = 0.0;
  double posLengthL1 = 1e-12;
  double colLengthL1 = 1e-12;
  double texLengthL1 = 1e-12;

  assert(grad1_scene.get_dmeshes().size() == 1 && grad2_scene.get_dmeshes().size() == 1);
  DMesh grad1 = grad1_scene.get_dmeshes()[0];
  DMesh grad2 = grad2_scene.get_dmeshes()[0];

  for(size_t i=0; i<3*grad1.vertex_count(); i++) 
  {
    if (!tested_mask[i])
      continue;
    double diff = ::std::abs(double(grad1.pos(i/3)[i%3] - grad2.pos(i/3)[i%3]));
    posError    += diff;
    posLengthL1 += ::std::abs(grad1.pos(i/3)[i%3]) + ::std::abs(grad2.pos(i/3)[i%3]);

    if (print)
    {
      ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << grad1.pos(i/3)[i%3] << "\t";  
      ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << grad2.pos(i/3)[i%3] << ::std::endl;
    }
  }

  if (print)
    ::std::cout << "--------------------------" << ::std::endl;
  if (grad1.color(0))
  {
    for(size_t i=0; i<3*grad1.vertex_count(); i++)
    {
      if (!tested_mask[3*grad1.vertex_count()+i])
        continue;
      double diff = ::std::abs(double(grad1.color(i/3)[i%3] - grad2.color(i/3)[i%3]));
      colError   += diff;
      colLengthL1 += ::std::abs(grad1.color(i/3)[i%3]) + ::std::abs(grad2.color(i/3)[i%3]);
      if (print)
      {
        ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << grad1.color(i/3)[i%3] << "\t";  
        ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << grad2.color(i/3)[i%3] << ::std::endl;
      }
    }
  }
  
  assert(grad1.full_size() == grad2.full_size());

  if (grad1.tex_count() > 0)
  {
    for (int t_n=0;t_n<grad1.tex_count();t_n++)
    {
      int3 p = grad1.get_tex_info(t_n);
      int off = grad1.tex(t_n, 0) - grad1.full_data();
      for (int i=0;i<p.x*p.y*p.z;i++)
      {
        if (!tested_mask[off+i])
          continue;
        double diff = ::std::abs(double(*(grad1.tex(t_n, i)) - *(grad2.tex(t_n, i))));
        texError += diff;
        texLengthL1 += ::std::abs(*(grad1.tex(t_n, i))) + ::std::abs(*(grad2.tex(t_n, i)));
      }
    }
  }

  double totalError = posError + colError + texError;
  double totalLengthL1 = texLengthL1 + posLengthL1 + colLengthL1;

  if (print)
  {
  ::std::cout << "==========================" << ::std::endl;
  ::std::cout << "GradErr[L1](vpos   ) = " << ::std::setw(10) << ::std::setprecision(4) << posError/double(grad1.vertex_count()*3)    << "\t which is \t" << 100.0*(posError/posLengthL1) << "%" << ::std::endl;
  ::std::cout << "GradErr[L1](color  ) = " << ::std::setw(10) << ::std::setprecision(4) << colError/double(grad1.vertex_count()*3)    << "\t which is \t" << 100.0*(colError/colLengthL1) << "%" << ::std::endl;
  ::std::cout << "GradErr[L1](texture) = " << ::std::setw(10) << ::std::setprecision(4) << texError                               << "\t which is \t" << 100.0*(texError/texLengthL1) << "%" << ::std::endl;
  ::std::cout << "GradErr[L1](average) = " << ::std::setw(10) << ::std::setprecision(4) << totalError/double(grad1.full_size())        << "\t which is \t" << 100.0*(totalError/totalLengthL1) << "%" << ::std::endl;
  }

  return DerivativesTestResults{posError/posLengthL1, colError/colLengthL1, texError/texLengthL1, totalError/totalLengthL1};
}
}