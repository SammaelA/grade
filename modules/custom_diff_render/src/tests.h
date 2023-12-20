#pragma once
#include "utils.h"

namespace diff_render
{
struct DScene;
struct TriangleMesh;
struct Scene;
struct IDiffRender;
struct CamInfo;
struct DiffRenderSettings;

class Tester
{
public:
  struct DerivativesTestResults
  {
    double pos_error = 0;
    double color_error = 0;
    double texture_error = 0;
    double average_error = 0;
  };
  static void test_base_derivatives();
  static void test_2_1_triangle();
  static void test_2_2_pyramid();
  static void test_2_3_sphere();
  static void test_2_4_pyramid_vcol();
  static void test_2_5_teapot_diffuse();
  static void test_2_6_path_tracing();
  static void test_2_7_mesh_on_static_scene();
  static void test_2_8_instancing();
  static void test_2_9_transform();
  static void test_2_10_multiple_meshes();
  static void test_2_11_restricted_transforms();

  static void test_3_1_mitsuba_triangle();
  static void test_3_2_mitsuba_sphere();
  static void test_3_3_mitsuba_teapot();
  static void test_3_4_mitsuba_cube();

  static DerivativesTestResults test_derivatives(const Scene &initial_scene, const Scene &target_scene, const CamInfo& a_camData, const DiffRenderSettings &settings, 
                                                 int max_test_vertices = 100, int max_test_texels = 100, bool print = false);
  static void test_fin_diff(const Scene &mesh, const char* outFolder, const Img& origin, const Img& target, ::std::shared_ptr<IDiffRender> a_pDRImpl, const CamInfo& a_camData,
                            DScene &d_mesh, int debug_mesh_id, int max_test_vertices, int max_test_texels,
                            ::std::vector<bool> &tested_mask);
  static DerivativesTestResults PrintAndCompareGradients(DScene& grad1_scene, DScene& grad2_scene, ::std::vector<bool> &tested_mask, bool print);
};
}