#include <boost/algorithm/string.hpp>
#include "tinyEngine/engine.h"
#include "tinyEngine/image.h"
#include "diff_geometry_generation.h"
#include "diff_optimization.h"
#include "mitsuba_python_interaction.h"
#include "graphics_utils/model_texture_creator.h"
#include "common_utils/optimization/optimization_benchmark.h"
#include "depth_extract_compare.h"
#include "graphics_utils/bilateral_filter.h"
#include "graphics_utils/resize_image.h"
#include "graphics_utils/unsharp_masking.h"
#include "graphics_utils/silhouette.h"
#include "graphics_utils/voxelization/voxelization.h"
#include "compare_utils.h"
#include <cppad/cppad.hpp>
#include <thread>
#include <chrono>
#include <time.h>
#include <csignal>
#include "obj_utils.h"
#include "graphics_utils/simple_model_utils.h"
#include "iou3d.h"
#include "common_utils/blk.h"
#include "common_utils/distribution.h"


void render_normalized(MitsubaInterface &mi, const dgen::DFModel &model, const std::string &texture_name, const std::string &folder_name,
                       int image_size, int spp, int rotations, float camera_dist, float camera_y, CameraSettings camera,
                       Block *camera_presets = nullptr,
                       MitsubaInterface::ModelInfo *m_info = nullptr,
                       MitsubaInterface::RenderStyle render_style = MitsubaInterface::TEXTURED_DEMO,
                       float y_rot = 0, bool render_mask = true)
{
  MitsubaInterface::ModelInfo model_info = MitsubaInterface::ModelInfo::simple_mesh(texture_name, mi.get_default_material());
  if (m_info)
    model_info = *m_info;
  
  if (render_style == MitsubaInterface::TEXTURED_DEMO && render_mask) 
  {
  mi.init_scene_and_settings(MitsubaInterface::RenderSettings(image_size, image_size, 32, MitsubaInterface::CUDA, MitsubaInterface::SILHOUETTE),
                             model_info);
  for (int i = 0; i < rotations; i++)
  {
    float phi = (2 * PI * i) / rotations;
    float psi;
    if (i <= rotations/2)
      psi = ((PI * i) / rotations)*y_rot;
    else
      psi = (PI - (PI * i) / rotations)*y_rot;
    camera.target = glm::vec3(0, 0, 0);
    camera.up = glm::vec3(0, 1, 0);
    camera.origin = glm::vec3(camera_dist * sin(phi) * cos(psi), camera_y + camera_dist * sin(psi), camera_dist * cos(phi) * cos(psi));

    char path[1024];
    sprintf(path, "%s/mask-%04d.png", folder_name.c_str(), i);
    //logerr("ss %s", path);
    std::vector<float> no_transform_scene_params = {0, 0, 0, 0, 0, 0, 100, 1000, 100, 100, 100, 0.01, camera.fov_rad};
    mi.render_model_to_file(model, path, camera, no_transform_scene_params);
    
    Block camera_blk;
    visualizer::save_camera_settings(camera, camera_blk);
    camera_blk.set_string("textured", std::string(path));
    if (camera_presets)
      camera_presets->add_block("camera", &camera_blk);
  }
  }
  
  mi.init_scene_and_settings(MitsubaInterface::RenderSettings(image_size, image_size, spp, MitsubaInterface::CUDA, render_style),
                             model_info);
  for (int i = 0; i < rotations; i++)
  {
    float phi = (2 * PI * i) / rotations;
    float psi;
    if (i < rotations/2)
      psi = ((PI * i) / rotations)*y_rot;
    else
      psi = (PI - (PI * i) / rotations)*y_rot;
    camera.target = glm::vec3(0, 0, 0);
    camera.up = glm::vec3(0, 1, 0);
    camera.origin = glm::vec3(camera_dist * sin(phi) * cos(psi), camera_y + camera_dist * sin(psi), camera_dist * cos(phi) * cos(psi));

    char path[1024];
    sprintf(path, "%s/frame-%04d.png", folder_name.c_str(), i);
    //logerr("ss %s", path);
    std::vector<float> no_transform_scene_params = {0, 0, 0, 0, 0, 0, 100, 1000, 100, 100, 100, 0.01, camera.fov_rad};
    mi.render_model_to_file(model, path, camera, no_transform_scene_params);
    
    Block camera_blk;
    visualizer::save_camera_settings(camera, camera_blk);
    camera_blk.set_string("textured", std::string(path));
    if (camera_presets)
      camera_presets->add_block("camera", &camera_blk);
  }
}

void render_mygen_cup(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                      std::string texture_path, std::string save_dir)
{
  dgen::DFModel res;
  dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
  visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0, 1, 0)));

  auto bbox = visualizer::get_bbox(res.first);
  visualizer::normalize_model(res.first);
  bbox = visualizer::get_bbox(res.first);

  bool ok = prepare_directory(save_dir+"mygen_turntable") && prepare_directory(save_dir+"mygen_turntable_9") && 
            prepare_directory(save_dir+"monochrome_turntable") && prepare_directory(save_dir+"monochrome_turntable_9");
  if (!ok)
  {
    logerr("unable to save turntables!");
    return;
  }

  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"mygen_turntable", 1024, 64, 64, 3, 0, camera, nullptr, nullptr, MitsubaInterface::TEXTURED_DEMO);
  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"mygen_turntable_9", 256, 64, 9, 3, 0, camera, nullptr, nullptr, MitsubaInterface::TEXTURED_DEMO);
  
  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"monochrome_turntable", 1024, 64, 64, 3, 0, camera, nullptr, nullptr, MitsubaInterface::MONOCHROME_DEMO);
  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"monochrome_turntable_9", 256, 64, 9, 3, 0, camera, nullptr, nullptr, MitsubaInterface::MONOCHROME_DEMO);
}

void render_mygen_cup_demo(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                           std::string texture_path, std::string save_dir)
{
  dgen::DFModel res;
  dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
  visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0, 1, 0)));

  auto bbox = visualizer::get_bbox(res.first);
  visualizer::normalize_model(res.first);
  bbox = visualizer::get_bbox(res.first);

  bool ok = prepare_directory(save_dir+"demo_textured") && prepare_directory(save_dir+"demo_monochrome");
  if (!ok)
  {
    logerr("unable to save turntables!");
    return;
  }
                      
  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"demo_textured", 1024, 1024, 99, 3, 0, camera, nullptr, nullptr, MitsubaInterface::TEXTURED_DEMO, 0.8, false);
  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"demo_monochrome", 1024, 1024, 99, 3, 0, camera, nullptr, nullptr, MitsubaInterface::MONOCHROME_DEMO, 0.8, false);
}

void render_mygen_cup_not_diff(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                           std::string texture_path, std::string save_dir)
{
  dgen::DFModel res;
  dgen::not_diff_gen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
  visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0, 1, 0)));

  auto bbox = visualizer::get_bbox(res.first);
  visualizer::normalize_model(res.first);
  bbox = visualizer::get_bbox(res.first);

  bool ok = prepare_directory(save_dir+"demo_textured") && prepare_directory(save_dir+"demo_monochrome");
  if (!ok)
  {
    logerr("unable to save turntables!");
    return;
  }
                      
  render_normalized(mi, res, "../../"+texture_path,
                    save_dir+"demo_textured", 512, 256, 8, 3, 0, camera, nullptr, nullptr, MitsubaInterface::TEXTURED_DEMO, 0.8, false);
  //render_normalized(mi, res, "../../"+texture_path,
  //                  save_dir+"demo_monochrome", 1024, 1024, 99, 3, 0, camera, nullptr, nullptr, MitsubaInterface::MONOCHROME_DEMO, 0.8, false);
}

void render_random_cameras(MitsubaInterface &mi, const dgen::DFModel &model, const std::string &texture_name, const std::string &folder_name,
                           int image_size, int spp, int images, float camera_dist, float camera_y, CameraSettings camera,
                           float phi_min, float phi_max, float psi_min, float psi_max, float jitter)
{
  MitsubaInterface::ModelInfo model_info = MitsubaInterface::ModelInfo::simple_mesh(texture_name, mi.get_default_material());
  mi.init_scene_and_settings(MitsubaInterface::RenderSettings(image_size, image_size, spp, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_DEMO),
                             model_info);
  for (int i = 0; i < images; i++)
  {
    float phi = urand(phi_min, phi_max);
    float psi = urand(psi_min, psi_max);
    camera.target = glm::vec3(urand(-jitter, jitter), urand(-jitter, jitter), urand(-jitter, jitter));
    camera.up = glm::vec3(0, 1, 0);
    camera.origin = glm::vec3(camera_dist * sin(phi) * cos(psi), camera_y + camera_dist * sin(psi), camera_dist * cos(phi) * cos(psi));

    char path[1024];
    sprintf(path, "%s/frame-%04d.png", folder_name.c_str(), i);
    logerr("ss %s", path);
    std::vector<float> no_transform_scene_params = {0, 0, 0, 0, 0, 0, 100, 1000, 100, 100, 100, 0.01, camera.fov_rad};
    mi.render_model_to_file(model, path, camera, no_transform_scene_params);
  }
}

void compare_and_print(std::string mygen_name, std::string diff_sdf_name)
{
  logerr("\n%s DiffProcGen\n", mygen_name.c_str());
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/reference_turntable",
                                "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/mygen_turntable",
                                64);
  logerr("\n%s Zero1to3\n", mygen_name.c_str());
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/reference_turntable_9",
                                "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/zero1to3_turntable_9",
                                9);
  logerr("\n%s Instant-NGP 4\n", mygen_name.c_str());
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP4_turntable",
                                "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref",
                                16);
  logerr("\n%s Instant-NGP 16\n", mygen_name.c_str());
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP16_turntable",
                                "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref",
                                16);
  logerr("\n%s Instant-NGP 64\n", mygen_name.c_str());
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP64_turntable",
                                "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref",
                                16);
  logerr("\n%s DiffSDF 1\n", mygen_name.c_str());
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-1/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-1/warp/reference/turntable",
                                64);
  logerr("\n%s DiffSDF 2\n", mygen_name.c_str());
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-2/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-2/warp/reference/turntable",
                                64);
  logerr("\n%s DiffSDF 6\n", mygen_name.c_str());
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-6/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-6/warp/reference/turntable",
                                64);
  logerr("\n%s DiffSDF 12\n", mygen_name.c_str());
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-12/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/" + diff_sdf_name + "/diffuse-12/warp/reference/turntable",
                                64);
}

void cup_1_render_reference_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
  auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/cup_model_1/cup_1.obj");
  auto bbox = visualizer::get_bbox(model);
  visualizer::normalize_model(model);
  bbox = visualizer::get_bbox(model);

  Block cameras_64, cameras_9;
  dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part", 0}}};
  //render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
  //                  "../grade_resources/prezentations/spring_23_medialab/cup_model_1/reference_turntable", 1024, 64, 64, 3, 0, camera, &cameras_64);
  //save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/cup_model_1/reference_turntable/cameras.blk", cameras_64);
  render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                    "../grade_resources/prezentations/spring_23_medialab/cup_model_1/reference_turntable_9", 1024, 256, 9, 3, 0, camera, &cameras_9);
  save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/cup_model_1/reference_turntable_9/cameras.blk", cameras_9);
}

void cup_1_render_img2mesh_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
  auto model = dgen::load_obj("saves/cup1_img2mesh/model.obj");
  dgen::set_face_normals(model);
  auto bbox = visualizer::get_bbox(model);logerr("model bbox 1 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
      
  visualizer::normalize_model(model);
  bbox = visualizer::get_bbox(model);logerr("model bbox 1 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
      

  Block cameras_64, cameras_9;
  dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part", 0}}};
  render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                    "saves/cup1_img2mesh/turntable", 1024, 64, 64, 3, 0, camera, &cameras_64);
  save_block_to_file("../saves/cup1_img2mesh/turntable/cameras.blk", cameras_64);
  render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                    "saves/cup1_img2mesh/turntable_9", 1024, 256, 9, 3, 0, camera, &cameras_9);
  save_block_to_file("../saves/cup1_img2mesh/turntable_9/cameras.blk", cameras_9);
}

void cup_1_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      std::vector<float> params = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.074, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
      dgen::DFModel res;
      dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
      visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
      
      auto bbox = visualizer::get_bbox(res.first);
      visualizer::normalize_model(res.first);
      bbox = visualizer::get_bbox(res.first);

      render_normalized(mi, res, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/reconstructed_tex_complemented.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_1/mygen_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, res, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/reconstructed_tex_complemented.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_1/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}

void cup_1_multicam_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      std::vector<float> params = {3.329, 3.515, 3.716, 3.740, 3.814, 3.858, 3.917, 3.961, 4.006, 0.654, 1.000, 0.059, 0.547, 0.258, 0.262, 0.302, 0.317, 0.318, 0.308, 0.294, 0.281, 0.268, 0.259, 0.252, 0.249, 0.249, 0.250, 0.253, 0.257, 0.264, 0.256, 0.264, 0.328, 1.543, 1.370, 1.123, 0.990, 0.943, 0.956, 0.947, 0.961, 0.952, 0.951, 0.936, 0.928, 0.939, 0.938, 0.971, 0.940, 0.998, 1.261, 1.720, 1.974, -0.138, 0.042, 0.020, -0.000, 3.119, 0.001, -9.294, 6.318, 662.677, 1.000, 79.108, 0.200, 0.250};
      dgen::DFModel res;
      dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
      visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
      
      auto bbox = visualizer::get_bbox(res.first);
      visualizer::normalize_model(res.first);
      bbox = visualizer::get_bbox(res.first);

      render_normalized(mi, res, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1_multicam/reconstructed_tex_complemented.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_1_multicam/mygen_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, res, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1_multicam/reconstructed_tex_complemented.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_1_multicam/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}

void cup_4_render_reference_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/cup_model_4/cup_4.obj");
      auto bbox = visualizer::get_bbox(model);
      visualizer::normalize_model(model);
      bbox = visualizer::get_bbox(model);

      Block cameras_64, cameras_9;
      dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part",0}}};
      render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_4/cup_4_tex_inv.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable", 1024, 64, 64, 3, 0, camera, &cameras_64);
      save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable/cameras.blk", cameras_64);
      render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_4/cup_4_tex_inv.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable_9", 256, 64, 9, 3, 0, camera, &cameras_9);
      save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable_9/cameras.blk", cameras_9);
}

void cup_4_render_reference_turntable_2(MitsubaInterface &mi, CameraSettings &camera)
{
      auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/cup_model_4/cup_4.obj");
      auto bbox = visualizer::get_bbox(model);
      visualizer::normalize_model(model);
      bbox = visualizer::get_bbox(model);

      Block cameras_64, cameras_9;
      dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part",0}}};
      render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_4/cup_4_tex_inv.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable_9_1024", 1024, 64, 9, 3, 0, camera, &cameras_9);
      save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable_9_1024/cameras.blk", cameras_9);
}

void cup_4_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      std::vector<float> params = {2.681, 3.268, 3.711, 3.896, 4.045, 4.068, 4.218, 4.540, 4.965, 0.906, 1.000, 0.043, 0.519, 0.101, 0.141, 0.183, 0.257, 0.283, 0.280, 0.274, 0.263, 0.252, 0.240, 0.223, 0.200, 0.171, 0.162, 0.155, 0.157, 0.166, 0.179, 0.207, 0.291, 1.091, 0.578, 1.950, 0.871, 0.805, 0.640, 0.674, 0.688, 0.720, 0.783, 0.810, 0.791, 0.582, 0.700, 0.620, 0.500, 0.552, 0.783, 1.271, 1.629, 0.082, 0.572, 0.402, 0.162, -0.090, 0.003, 0.590, 12.159, 666.696, 1.000, 79.108, 0.200, 0.250};
      dgen::DFModel res;
      dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
      visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
      
      auto bbox = visualizer::get_bbox(res.first);
      visualizer::normalize_model(res.first);
      bbox = visualizer::get_bbox(res.first);

      render_normalized(mi, res, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_4/reconstructed_tex_complemented.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_4/mygen_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, res, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_4/reconstructed_tex_complemented.png",
                        "../grade_resources/prezentations/spring_23_medialab/cup_model_4/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}

void building_2_render_reference_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/test_building_2/original/original.obj");
      auto bbox = visualizer::get_bbox(model);
      logerr("model bbox 1 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
      visualizer::normalize_model(model);
      bbox = visualizer::get_bbox(model);
      logerr("model bbox 2 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);

      dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part",0}}};
      Block cameras_64, cameras_9;
      render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/test_building_2/original/tex_inv.png",
                        "../grade_resources/prezentations/spring_23_medialab/test_building_2/reference_turntable", 1024, 64, 64, 3, 0, camera, &cameras_64);
      save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/test_building_2/reference_turntable/cameras.blk", cameras_64);
      render_normalized(mi, df_model, "../../../grade_resources/prezentations/spring_23_medialab/test_building_2/original/tex_inv.png",
                        "../grade_resources/prezentations/spring_23_medialab/test_building_2/reference_turntable_9", 256, 64, 9, 3, 0, camera, &cameras_9);
      save_block_to_file("../../grade_resources/prezentations/spring_23_medialab/test_building_2/reference_turntable_9/cameras.blk", cameras_9);
}
void building_2_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("buildings_2").generator_description_blk_path, gen_info);
    Block &gen_mesh_parts = *gen_info.get_block("mesh_parts");
    MitsubaInterface::ModelInfo model_info;
    model_info.layout = dgen::ModelLayout(0, 3, 6, 8, 8);//default layout with pos, normals and tc

    for (int i=0;i<gen_mesh_parts.size();i++)
    {
      if (gen_mesh_parts.get_type(i) == Block::ValueType::STRING)
        model_info.parts.push_back({gen_mesh_parts.get_string(i), 
                                    "white.png", 
                                    MitsubaInterface::get_default_material()});
    }

    model_info.get_part("main_part")->texture_name = "wall3.png";
    model_info.get_part("interior")->texture_name = "concrete.png";
    model_info.get_part("windows")->material_name = "glass";
    model_info.get_part("wooden_parts")->texture_name = "wood6.png";
    model_info.get_part("metal_parts")->texture_name = "rusty_metal.png";
    model_info.get_part("roof")->texture_name = "roof1.png";
    std::vector<float> params = {1.000, 1.853, 0.040, 0.322, 0.005, 2.000, 5.000, 1.000, 5.000, 650.000, 1.000, 5.000, 341.000, 1.000, 0.080, 0.080, 0.100, 0.008, 0.410, 0.000, 0.024, 0.315, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 2.174, 1.000, 3.000, 0.000, 0.600, 0.486, 0.600, 1.484, 0.400, 0.500, 0.015, 0.050, 1.000, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 1.721, 0.064, 0.478, 0.737, 0.150, 0.100, 0.200, 0.416, 1.000, -0.160, -0.049, 0.954, 0.072, 0.523, 0.008, 0.000, 0.500, 10.000, 1.000, 100.000, 0.100, 0.300};
    dgen::DFModel res;
    dgen::dgen_test("buildings_2", params, res, false, dgen::ModelQuality(false, 3));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    auto res_bbox = visualizer::get_bbox(res.first);
    visualizer::normalize_model(res.first);
    render_normalized(mi, res, "../mitsuba_data/meshes/building/tex2.png",
                      "../grade_resources/prezentations/spring_23_medialab/test_building_2/mygen_turntable", 256, 64, 64, 3, 0, camera, nullptr, &model_info);
    //render_normalized(mi, res, "../mitsuba_data/meshes/building/tex2.png",
    //                  "../grade_resources/prezentations/spring_23_medialab/test_building_2/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}
void building_6_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("buildings_2").generator_description_blk_path, gen_info);
    Block &gen_mesh_parts = *gen_info.get_block("mesh_parts");
    MitsubaInterface::ModelInfo model_info;
    model_info.layout = dgen::ModelLayout(0, 3, 6, 8, 8);//default layout with pos, normals and tc

    for (int i=0;i<gen_mesh_parts.size();i++)
    {
      if (gen_mesh_parts.get_type(i) == Block::ValueType::STRING)
        model_info.parts.push_back({gen_mesh_parts.get_string(i), 
                                    "white.png", 
                                    MitsubaInterface::get_default_material()});
    }

    model_info.get_part("main_part")->texture_name = "wall2.png";
    model_info.get_part("interior")->texture_name = "concrete.png";
    model_info.get_part("windows")->material_name = "glass";
    model_info.get_part("wooden_parts")->texture_name = "wood6.png";
    model_info.get_part("metal_parts")->texture_name = "rusty_metal.png";
    model_info.get_part("roof")->texture_name = "roof1.png";
    std::vector<float> params = {2.000000, 4.000000, 0.040000, 0.375206, 0.336109, 2.000000, 9.000000, 0.330000, 6.000000, 1625.000000, 0.500000, 5.000000, 68.000000, 2.000000, 0.080000, 0.080000, 0.095000, 0.008677, 0.179252, 0.000000, 0.003194, 0.178314, 1.000000, 2.000000, 2.000000, 0.600000, 0.519676, 0.600000, 1.824628, 1.000000, 3.000000, 0.000000, 0.581000, 0.574231, 0.600000, 2.903548, 0.400000, 0.900000, 0.018000, 0.045000, 1.000000, 1.000000, 2.000000, 2.000000, 0.600000, 0.400000, 0.600000, 2.767931, 0.569821, 0.400000, 0.887747, 0.150000, 0.116000, 0.150000, 0.954397, 0.379385, 0.207553, 0.415835, -0.392796, 0.708893, -0.100000, 3.509702, 0.011455, 0.000000, 10.500000, 10.000000, 1.000000, 1.000000, 0.100000, 0.500000};
    dgen::DFModel res, res_det;
    dgen::dgen_test("buildings_2", params, res, false, dgen::ModelQuality(false, 0));
    dgen::dgen_test("buildings_2", params, res_det, false, dgen::ModelQuality(false, 3));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    visualizer::normalize_model(res.first);

    visualizer::transform(res_det.first, LiteMath::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    visualizer::normalize_model(res_det.first);
    render_normalized(mi, res_det, "../../saves/building_6_res/reconstructed_tex_raw.png",
                      "saves/building_6_res/monochrome_turntable", 1024, 256, 64, 3, 0, camera, nullptr, &model_info);
    render_normalized(mi, res_det, "../../saves/building_6_res/reconstructed_tex_raw.png",
                      "saves/building_6_res/monochrome_turntable_9", 1024, 256, 9, 3, 0, camera, nullptr, &model_info);

    render_normalized(mi, res, "../../saves/building_6_res/reconstructed_tex_raw.png",
                      "saves/building_6_res/mygen_turntable", 1024, 256, 64, 3, 0, camera);
    render_normalized(mi, res, "../../saves/building_6_res/reconstructed_tex_raw.png",
                      "saves/building_6_res/mygen_turntable_9", 1024, 256, 9, 3, 0, camera);
}
  /*
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("buildings_2").generator_description_blk_path, gen_info);
    Block &gen_mesh_parts = *gen_info.get_block("mesh_parts");

    std::vector<float> params = {1.000, 1.853, 0.040, 0.322, 0.005, 2.000, 5.000, 1.000, 5.000, 650.000, 1.000, 5.000, 341.000, 1.000, 0.080, 0.080, 0.100, 0.008, 0.410, 0.000, 0.024, 0.315, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 2.174, 1.000, 3.000, 0.000, 0.600, 0.486, 0.600, 1.484, 0.400, 0.500, 0.015, 0.050, 1.000, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 1.721, 0.064, 0.478, 0.737, 0.150, 0.100, 0.200, 0.416, 1.000, -0.160, -0.049, 0.954, 0.072, 0.523, 0.008, 0.000, 0.500, 10.000, 1.000, 100.000, 0.100, 0.300};
    dgen::DFModel res;
    dgen::dgen_test("buildings_2", params, res, false, dgen::ModelQuality(false, 0));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    
    auto res_bbox = visualizer::get_bbox(res.first);
    logerr("model bbox 3 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    visualizer::normalize_model(res.first);
    //for (int i=0;i<model.size();i++)
    //  logerr("%d %f",i, model[i]);
    res_bbox = visualizer::get_bbox(res.first);
    glm::vec3 sz1 = bbox.max_pos - bbox.min_pos;
    glm::vec3 sz2 = res_bbox.max_pos - res_bbox.min_pos;
    dgen::scale(res.first, sz1/sz2);
    logerr("model bbox 4 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);


    //render_normalized(mi, df_model, "../mitsuba_data/meshes/building/tex2.png",
    //                  "saves/building_2/reference", 256, 256, 64, 3, 0, camera);
    render_normalized(mi, res, "../mitsuba_data/meshes/building/tex2.png",
                      "saves/building_2/result", 256, 4, 64, 3, 0, camera);

    compare_utils::turntable_loss("/home/sammael/grade/saves/building_2/reference",
                                  "/home/sammael/grade/saves/building_2/result",
                                  64);

    compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-6/warp/turntable",
                                  "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-6/warp/reference/turntable",
                                  64);
    compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-12/warp/turntable",
                                  "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-12/warp/reference/turntable",
                                  64);
  */


void test_3diou_1()
{
    std::vector<float> params = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.074, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res;
    dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    visualizer::normalize_model(res.first);

    std::vector<float> params2 = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.07, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res2;
    dgen::dgen_test("dishes", params2, res2, false, dgen::ModelQuality(false, 2));
    visualizer::transform(res2.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    visualizer::normalize_model(res2.first);

    float iou = iou3d(res2.first, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/32);
    logerr("IoU %f", iou);
}

void cup_1_iou()
{
    std::vector<float> params = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.2, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res;
    dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    visualizer::normalize_model(res.first);

    auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/cup_model_1/cup_1.obj");
    visualizer::normalize_model(model);

    float iou = iou3d(model, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/64);
    logerr("Cup 1 IoU %f", iou);
}

void cup_4_iou()
{
    std::vector<float> params = {2.681, 3.268, 3.711, 3.896, 4.045, 4.068, 4.218, 4.540, 4.965, 0.906, 1.000, 0.043, 0.519, 0.101, 0.141, 0.183, 0.257, 0.283, 0.280, 0.274, 0.263, 0.252, 0.240, 0.223, 0.200, 0.171, 0.162, 0.155, 0.157, 0.166, 0.179, 0.207, 0.291, 1.091, 0.578, 1.950, 0.871, 0.805, 0.640, 0.674, 0.688, 0.720, 0.783, 0.810, 0.791, 0.582, 0.700, 0.620, 0.500, 0.552, 0.783, 1.271, 1.629, 0.082, 0.572, 0.402, 0.162, -0.090, 0.003, 0.590, 12.159, 666.696, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res;
    dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    visualizer::normalize_model(res.first);

    auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/cup_model_4/cup_4.obj");
    visualizer::normalize_model(model);

    float iou = iou3d(model, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/32);
    logerr("Cup 4 IoU %f", iou);
}

void building_2_iou()
{
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("buildings_2").generator_description_blk_path, gen_info);
    Block &gen_mesh_parts = *gen_info.get_block("mesh_parts");

    std::vector<float> params = {1.000, 1.853, 0.040, 0.322, 0.005, 2.000, 5.000, 1.000, 5.000, 650.000, 1.000, 5.000, 341.000, 1.000, 0.080, 0.080, 0.100, 0.008, 0.410, 0.000, 0.024, 0.315, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 2.174, 1.000, 3.000, 0.000, 0.600, 0.486, 0.600, 1.484, 0.400, 0.500, 0.015, 0.050, 1.000, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 1.721, 0.064, 0.478, 0.737, 0.150, 0.100, 0.200, 0.416, 1.000, -0.160, -0.049, 0.954, 0.072, 0.523, 0.008, 0.000, 0.500, 10.000, 1.000, 100.000, 0.100, 0.300};
    dgen::DFModel res;
    dgen::dgen_test("buildings_2", params, res, false, dgen::ModelQuality(false, 0));
    visualizer::transform(res.first, LiteMath::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    auto res_bbox = visualizer::get_bbox(res.first);
    visualizer::normalize_model(res.first);

    auto model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/test_building_2/original/original.obj");
    visualizer::normalize_model(model);

    float iou = iou3d(model, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/32);
    logerr("Building 2 IoU %f", iou);
}


void render_NGP_turntable(MitsubaInterface &mi, CameraSettings &camera, std::string mygen_name)
{
    auto model4 = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t4_simple.obj");
    auto model16 = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t16_simple.obj");
    auto model64 = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t64_simple.obj");
    // dgen::shift(model, glm::vec3(0,0,0.5));
    auto bbox = visualizer::get_bbox(model4);
    auto ref_model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/ref.obj");
    visualizer::transform(ref_model, LiteMath::rotate(glm::mat4(1.0f), -PI / 2, glm::vec3(0, 1, 0)));
    auto ref_bbox = visualizer::get_bbox(ref_model);
    logerr("model bbox 4 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    logerr("ref   bbox 4 (%f %f %f)(%f %f %f)", ref_bbox.min_pos.x, ref_bbox.min_pos.y, ref_bbox.min_pos.z,
           ref_bbox.max_pos.x, ref_bbox.max_pos.y, ref_bbox.max_pos.z);

    glm::vec3 sizes = ref_bbox.max_pos - ref_bbox.min_pos;
    float max_size = MAX(sizes.x, MAX(sizes.y, sizes.z));
    max_size = MAX(1e-6, max_size);
    visualizer::shift(ref_model, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(ref_model, glm::vec3(1 / max_size));
    visualizer::shift(model4, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(model4, glm::vec3(1 / max_size));
    visualizer::shift(model16, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(model16, glm::vec3(1 / max_size));
    visualizer::shift(model64, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(model64, glm::vec3(1 / max_size));
    // visualizer::normalize_model(model);
    // bbox = visualizer::get_bbox(model);

    dgen::DFModel df_model4 = {model4, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, df_model4, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP4_turntable", 1024, 64, 16, 3, 0, camera);
    dgen::DFModel df_model16 = {model16, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, df_model16, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP16_turntable", 1024, 64, 16, 3, 0, camera);
    dgen::DFModel df_model64 = {model64, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, df_model64, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP64_turntable", 1024, 64, 16, 3, 0, camera);
    dgen::DFModel ref_df_model = {ref_model, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, ref_df_model, "../../../grade_resources/prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref", 1024, 64, 16, 3, 0, camera);
}

void calc_NGP_3D_IoU(MitsubaInterface &mi, CameraSettings &camera, std::string mygen_name)
{
    auto model4 = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t4_simple.obj");
    auto model16 = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t16_simple.obj");
    auto model64 = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t64_simple.obj");
    // dgen::shift(model, glm::vec3(0,0,0.5));
    auto bbox = visualizer::get_bbox(model4);
    auto ref_model = dgen::load_obj("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/ref.obj");
    visualizer::transform(ref_model, LiteMath::rotate(glm::mat4(1.0f), -PI / 2, glm::vec3(0, 1, 0)));
    auto ref_bbox = visualizer::get_bbox(ref_model);
    logerr("model bbox 4 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    logerr("ref   bbox 4 (%f %f %f)(%f %f %f)", ref_bbox.min_pos.x, ref_bbox.min_pos.y, ref_bbox.min_pos.z,
           ref_bbox.max_pos.x, ref_bbox.max_pos.y, ref_bbox.max_pos.z);

    glm::vec3 sizes = ref_bbox.max_pos - ref_bbox.min_pos;
    float max_size = MAX(sizes.x, MAX(sizes.y, sizes.z));
    max_size = MAX(1e-6, max_size);
    visualizer::shift(ref_model, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(ref_model, glm::vec3(1 / max_size));
    visualizer::shift(model4, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(model4, glm::vec3(1 / max_size));
    visualizer::shift(model16, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(model16, glm::vec3(1 / max_size));
    visualizer::shift(model64, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    visualizer::scale(model64, glm::vec3(1 / max_size));

    float iou1 = iou3d(ref_model, model4, -1, -1, -1, 1, 1, 1, 1.0/16);
    float iou2 = iou3d(ref_model, model16, -1, -1, -1, 1, 1, 1, 1.0/16);
    float iou3 = iou3d(ref_model, model64, -1, -1, -1, 1, 1, 1, 1.0/16);

    logerr("%s NGP IoU %f %f %f", mygen_name.c_str(), iou1, iou2, iou3);
}

void compare_all_for_paper()
{
  logerr("mygen");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_1/reference_turntable",
                                "saves/cup_model_1/mygen_turntable",
                                64);
  logerr("image2mesh");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_1/reference_turntable",
                                "saves/cup1_img2mesh/turntable",
                                64);
  logerr("instant ngp");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_1/NGP64_turntable",
                                "../grade_resources/prezentations/spring_23_medialab/cup_model_1/NGP_turntable_ref",
                                16);
  logerr("diff SDF 1");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/cup_1/diffuse-1/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/cup_1/diffuse-1/warp/reference/turntable",
                                64);
  logerr("diff SDF 6");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/cup_1/diffuse-6/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/cup_1/diffuse-6/warp/reference/turntable",
                                64);

  logerr("cup_4 1 view");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable",
                                "saves/cup_model_4/mygen_turntable",
                                64);
  logerr("cup_4 2 view");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable",
                                "saves/cup_4_2cameras/mygen_turntable",
                                64);
  logerr("cup_4 4 view");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable",
                                "saves/cup_4_4cameras/mygen_turntable",
                                64);
  logerr("cup_4 8 view");
  compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/cup_model_4/reference_turntable",
                                "saves/cup_4_8cameras/mygen_turntable",
                                64);
}

void compare_sandbox(int argc, char **argv)
{
  MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");

  float fov_rad = 0.5;
  CameraSettings camera = MitsubaInterface::get_camera_from_scene_params({fov_rad});
  
  //compare_all_for_paper();
  
  logerr("mygen tree");
  compare_utils::turntable_loss("saves/tree_1/reference",
                                "saves/tree_1/mygen",
                                64);
  logerr("diff SDF 2");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_1/diffuse-2/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_1/diffuse-2/warp/reference/turntable",
                                64);
  logerr("diff SDF 6");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_1/diffuse-6/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_1/diffuse-6/warp/reference/turntable",
                                64);

  logerr("diff SDF 12");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_1/diffuse-12/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_1/diffuse-12/warp/reference/turntable",
                                64);

  logerr("mygen tree 2");
  compare_utils::turntable_loss("saves/tree_2/reference",
                                "saves/tree_2/mygen",
                                64);
  logerr("diff SDF 2");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_2/diffuse-2/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_2/diffuse-2/warp/reference/turntable",
                                64);
  logerr("diff SDF 6");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_2/diffuse-6/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_2/diffuse-6/warp/reference/turntable",
                                64);

  logerr("diff SDF 12");
  compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_2/diffuse-12/warp/turntable",
                                "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/tree_2/diffuse-12/warp/reference/turntable",
                                64);

  //cup_4_render_reference_turntable_2(mi, camera);
  //cup_1_render_img2mesh_turntable(mi, camera);
  //building_6_render_mygen_turntable(mi, camera);
  //cup_1_multicam_render_mygen_turntable(mi, camera);
  //std::string mygen_name = "cup_model_1_multicam";
  //compare_utils::turntable_loss("../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/reference_turntable",
  //                             "../grade_resources/prezentations/spring_23_medialab/" + mygen_name + "/mygen_turntable",
  //                              64);
  //building_2_render_mygen_turntable(mi, camera);

  //compare_and_print("test_building_2", "building_2");
  //compare_and_print("cup_model_4", "cup_4");
  //compare_and_print("cup_model_1", "cup_1");
  //building_2_render_mygen_turntable(mi, camera);
  //building_2_render_reference_turntable(mi, camera);
  //cup_4_render_reference_turntable(mi, camera);
  //cup_1_render_reference_turntable(mi, camera);
  //calc_NGP_3D_IoU(mi, camera, "cup_model_4"); 
  //calc_NGP_3D_IoU(mi, camera, "cup_model_1"); 
  //calc_NGP_3D_IoU(mi, camera, "test_building_2"); 
  //cup_1_iou();
  //cup_4_iou();
  //building_2_iou();
  //test_3diou_1();
  //test_3diou_2(); 
  //render_NGP_turntable(mi, camera, "cup_model_4"); 
  //render_NGP_turntable(mi, camera, "test_building_2");
}
/*
cup_model_4 NGP IoU 0.121806 0.202790 0.269309
cup_model_1 NGP IoU 0.299897 0.308010 0.456000
test_building_2 NGP IoU 0.134880 0.179714 0.153565
Cup 1 IoU 0.272983
Cup 4 IoU 0.422604
Cup 4 IoU 0.487518
test_building_2 DiffProcGen

Turntable Loss for 64 images
Textured PSNR = 18.67
With silhouette PSNR = 15.67
Silhouette IoU = 0.8856

test_building_2 Zero1to3

Turntable Loss for 9 images
Textured PSNR = 16.15
With silhouette PSNR = 12.70
Silhouette IoU = 0.6264

test_building_2 Instant-NGP 4

Turntable Loss for 16 images
Textured PSNR = 13.94
With silhouette PSNR = 12.72
Silhouette IoU = 0.7982

test_building_2 Instant-NGP 16

Turntable Loss for 16 images
Textured PSNR = 21.27
With silhouette PSNR = 19.94
Silhouette IoU = 0.9396

test_building_2 Instant-NGP 64

Turntable Loss for 16 images
Textured PSNR = 23.43
With silhouette PSNR = 19.14
Silhouette IoU = 0.9713

test_building_2 DiffSDF 1

Turntable Loss for 64 images
Textured PSNR = 11.86
With silhouette PSNR = 13.66
Silhouette IoU = 0.3215

test_building_2 DiffSDF 2

Turntable Loss for 64 images
Textured PSNR = 14.86
With silhouette PSNR = 13.93
Silhouette IoU = 0.7275

test_building_2 DiffSDF 6

Turntable Loss for 64 images
Textured PSNR = 23.05
With silhouette PSNR = 18.36
Silhouette IoU = 0.9666

test_building_2 DiffSDF 12

Turntable Loss for 64 images
Textured PSNR = 25.19
With silhouette PSNR = 20.20
Silhouette IoU = 0.9790

cup_model_4 DiffProcGen

Turntable Loss for 64 images
Textured PSNR = 22.60
With silhouette PSNR = 16.62
Silhouette IoU = 0.9669

cup_model_4 Zero1to3

Turntable Loss for 9 images
Textured PSNR = 13.87
With silhouette PSNR = 10.31
Silhouette IoU = 0.3332

cup_model_4 Instant-NGP 4

Turntable Loss for 16 images
Textured PSNR = 9.88
With silhouette PSNR = 11.94
Silhouette IoU = 0.5611

cup_model_4 Instant-NGP 16

Turntable Loss for 16 images
Textured PSNR = 16.69
With silhouette PSNR = 14.18
Silhouette IoU = 0.8779

cup_model_4 Instant-NGP 64

Turntable Loss for 16 images
Textured PSNR = 27.38
With silhouette PSNR = 21.60
Silhouette IoU = 0.9910

cup_model_4 DiffSDF 1

Turntable Loss for 64 images
Textured PSNR = 19.62
With silhouette PSNR = 14.37
Silhouette IoU = 0.3888

cup_model_4 DiffSDF 2

Turntable Loss for 64 images
Textured PSNR = 18.45
With silhouette PSNR = 13.90
Silhouette IoU = 0.4226

cup_model_4 DiffSDF 6

Turntable Loss for 64 images
Textured PSNR = 21.84
With silhouette PSNR = 17.15
Silhouette IoU = 0.5018

cup_model_4 DiffSDF 12

Turntable Loss for 64 images
Textured PSNR = 28.57
With silhouette PSNR = 23.54
Silhouette IoU = 0.6692

cup_model_1 DiffProcGen

Turntable Loss for 64 images
Textured PSNR = 24.67
With silhouette PSNR = 23.56
Silhouette IoU = 0.9404

cup_model_1 Zero1to3

Turntable Loss for 9 images
Textured PSNR = 13.70
With silhouette PSNR = 15.62
Silhouette IoU = 0.2251

cup_model_1 Instant-NGP 4

Turntable Loss for 16 images
Textured PSNR = 15.02
With silhouette PSNR = 13.08
Silhouette IoU = 0.7993

cup_model_1 Instant-NGP 16

Turntable Loss for 16 images
Textured PSNR = 18.04
With silhouette PSNR = 19.26
Silhouette IoU = 0.8577

cup_model_1 Instant-NGP 64

Turntable Loss for 16 images
Textured PSNR = 20.34
With silhouette PSNR = 23.31
Silhouette IoU = 0.9053

cup_model_1 DiffSDF 1

Turntable Loss for 64 images
Textured PSNR = 9.10
With silhouette PSNR = 16.40
Silhouette IoU = 0.5430

cup_model_1 DiffSDF 2

Turntable Loss for 64 images
Textured PSNR = 12.57
With silhouette PSNR = 16.99
Silhouette IoU = 0.7866

cup_model_1 DiffSDF 6

Turntable Loss for 64 images
Textured PSNR = 24.96
With silhouette PSNR = 23.80
Silhouette IoU = 0.9845

cup_model_1 DiffSDF 12

Turntable Loss for 64 images
Textured PSNR = 28.67
With silhouette PSNR = 26.27
Silhouette IoU = 0.9911
(base) sammael@sammael-520:~/grade$ 
*/