#include "sandbox.h"
#include <boost/algorithm/string.hpp>
#include "generation/scene_generation.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"
#include "tinyEngine/engine.h"
#include "tinyEngine/image.h"
#include "parameter_selection/impostor_similarity.h"
#include "parameter_selection/genetic_algorithm.h"
#include "tree_generators/GE_generator.h"
#include "parameter_selection/parameter_selection.h"
#include "tree_generators/all_generators.h"
#include "tree_generators/weber_penn_generator.h"
#include "parameter_selection/neural_selection.h"
#include "diff_generators/diff_geometry_generation.h"
#include "diff_generators/diff_optimization.h"
#include "diff_generators/mitsuba_python_interaction.h"
#include "graphics_utils/model_texture_creator.h"
#include "common_utils/optimization/optimization_benchmark.h"
#include "diff_generators/depth_extract_compare.h"
#include "graphics_utils/bilateral_filter.h"
#include "graphics_utils/resize_image.h"
#include "graphics_utils/unsharp_masking.h"
#include "graphics_utils/silhouette.h"
#include "graphics_utils/voxelization/voxelization.h"
#include "diff_generators/compare_utils.h"
#include <cppad/cppad.hpp>
#include <thread>
#include <chrono>
#include <time.h>
#include <csignal>
#include "diff_generators/obj_utils.h"
#include "diff_generators/simple_model_utils.h"
#include "diff_generators/iou3d.h"

void render_normalized(MitsubaInterface &mi, const dgen::DFModel &model, const std::string &texture_name, const std::string &folder_name,
                       int image_size, int spp, int rotations, float camera_dist, float camera_y, CameraSettings camera,
                       MitsubaInterface::ModelInfo *m_info = nullptr)
{
  MitsubaInterface::ModelInfo model_info = MitsubaInterface::ModelInfo::simple_mesh(texture_name, mi.get_default_material());
  if (m_info)
    model_info = *m_info;
  mi.init_scene_and_settings(MitsubaInterface::RenderSettings(image_size, image_size, spp, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_DEMO),
                             model_info);
  for (int i = 0; i < rotations; i++)
  {
    float phi = (2 * PI * i) / rotations;
    camera.target = glm::vec3(0, 0, 0);
    camera.up = glm::vec3(0, 1, 0);
    camera.origin = glm::vec3(camera_dist * sin(phi), camera_y, camera_dist * cos(phi));

    char path[1024];
    sprintf(path, "%s/frame-%04d.png", folder_name.c_str(), i);
    logerr("ss %s", path);
    std::vector<float> no_transform_scene_params = {0, 0, 0, 0, 0, 0, 100, 1000, 100, 100, 100, 0.01, camera.fov_rad};
    mi.render_model_to_file(model, path, camera, no_transform_scene_params);
  }
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
  compare_utils::turntable_loss("prezentations/spring_23_medialab/" + mygen_name + "/reference_turntable",
                                "prezentations/spring_23_medialab/" + mygen_name + "/mygen_turntable",
                                64);
  logerr("\n%s Zero1to3\n", mygen_name.c_str());
  compare_utils::turntable_loss("prezentations/spring_23_medialab/" + mygen_name + "/reference_turntable_9",
                                "prezentations/spring_23_medialab/" + mygen_name + "/zero1to3_turntable_9",
                                9);
  logerr("\n%s Instant-NGP 4\n", mygen_name.c_str());
  compare_utils::turntable_loss("prezentations/spring_23_medialab/" + mygen_name + "/NGP4_turntable",
                                "prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref",
                                16);
  logerr("\n%s Instant-NGP 16\n", mygen_name.c_str());
  compare_utils::turntable_loss("prezentations/spring_23_medialab/" + mygen_name + "/NGP16_turntable",
                                "prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref",
                                16);
  logerr("\n%s Instant-NGP 64\n", mygen_name.c_str());
  compare_utils::turntable_loss("prezentations/spring_23_medialab/" + mygen_name + "/NGP64_turntable",
                                "prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref",
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
  auto model = dgen::load_obj("prezentations/spring_23_medialab/cup_model_1/cup_1.obj");
  auto bbox = dgen::get_bbox(model);
  dgen::normalize_model(model);
  bbox = dgen::get_bbox(model);

  dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part", 0}}};
  render_normalized(mi, df_model, "../../prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                    "prezentations/spring_23_medialab/cup_model_1/reference_turntable", 1024, 64, 64, 3, 0, camera);
  render_normalized(mi, df_model, "../../prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                    "prezentations/spring_23_medialab/cup_model_1/reference_turntable_9", 256, 64, 9, 3, 0, camera);
}

void cup_1_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      std::vector<float> params = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.074, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
      dgen::DFModel res;
      dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
      dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
      
      auto bbox = dgen::get_bbox(res.first);
      dgen::normalize_model(res.first);
      bbox = dgen::get_bbox(res.first);

      render_normalized(mi, res, "../../prezentations/spring_23_medialab/cup_model_1/reconstructed_tex_complemented.png",
                        "prezentations/spring_23_medialab/cup_model_1/mygen_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, res, "../../prezentations/spring_23_medialab/cup_model_1/reconstructed_tex_complemented.png",
                        "prezentations/spring_23_medialab/cup_model_1/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}

void cup_4_render_reference_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      auto model = dgen::load_obj("prezentations/spring_23_medialab/cup_model_4/cup_4.obj");
      auto bbox = dgen::get_bbox(model);
      dgen::normalize_model(model);
      bbox = dgen::get_bbox(model);

      dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part",0}}};
      render_normalized(mi, df_model, "../../prezentations/spring_23_medialab/cup_model_4/cup_4_tex_inv.png",
                        "prezentations/spring_23_medialab/cup_model_4/reference_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, df_model, "../../prezentations/spring_23_medialab/cup_model_4/cup_4_tex_inv.png",
                        "prezentations/spring_23_medialab/cup_model_4/reference_turntable_9", 256, 64, 9, 3, 0, camera);
}

void cup_4_render_mygen_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      std::vector<float> params = {2.681, 3.268, 3.711, 3.896, 4.045, 4.068, 4.218, 4.540, 4.965, 0.906, 1.000, 0.043, 0.519, 0.101, 0.141, 0.183, 0.257, 0.283, 0.280, 0.274, 0.263, 0.252, 0.240, 0.223, 0.200, 0.171, 0.162, 0.155, 0.157, 0.166, 0.179, 0.207, 0.291, 1.091, 0.578, 1.950, 0.871, 0.805, 0.640, 0.674, 0.688, 0.720, 0.783, 0.810, 0.791, 0.582, 0.700, 0.620, 0.500, 0.552, 0.783, 1.271, 1.629, 0.082, 0.572, 0.402, 0.162, -0.090, 0.003, 0.590, 12.159, 666.696, 1.000, 79.108, 0.200, 0.250};
      dgen::DFModel res;
      dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
      dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
      
      auto bbox = dgen::get_bbox(res.first);
      dgen::normalize_model(res.first);
      bbox = dgen::get_bbox(res.first);

      render_normalized(mi, res, "../../prezentations/spring_23_medialab/cup_model_4/reconstructed_tex_complemented.png",
                        "prezentations/spring_23_medialab/cup_model_4/mygen_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, res, "../../prezentations/spring_23_medialab/cup_model_4/reconstructed_tex_complemented.png",
                        "prezentations/spring_23_medialab/cup_model_4/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}

void building_2_render_reference_turntable(MitsubaInterface &mi, CameraSettings &camera)
{
      auto model = dgen::load_obj("prezentations/spring_23_medialab/test_building_2/original/original.obj");
      auto bbox = dgen::get_bbox(model);
      logerr("model bbox 1 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
      dgen::normalize_model(model);
      bbox = dgen::get_bbox(model);
      logerr("model bbox 2 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);

      dgen::DFModel df_model = {model, dgen::PartOffsets{{"main_part",0}}};
      render_normalized(mi, df_model, "../../prezentations/spring_23_medialab/test_building_2/original/tex_inv.png",
                        "prezentations/spring_23_medialab/test_building_2/reference_turntable", 1024, 64, 64, 3, 0, camera);
      render_normalized(mi, df_model, "../../prezentations/spring_23_medialab/test_building_2/original/tex_inv.png",
                        "prezentations/spring_23_medialab/test_building_2/reference_turntable_9", 256, 64, 9, 3, 0, camera);
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
    dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    auto res_bbox = dgen::get_bbox(res.first);
    dgen::normalize_model(res.first);
    render_normalized(mi, res, "../mitsuba_data/meshes/building/tex2.png",
                      "prezentations/spring_23_medialab/test_building_2/mygen_turntable", 256, 64, 64, 3, 0, camera, &model_info);
    //render_normalized(mi, res, "../mitsuba_data/meshes/building/tex2.png",
    //                  "prezentations/spring_23_medialab/test_building_2/mygen_turntable_9", 256, 64, 9, 3, 0, camera);
}
  /*
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("buildings_2").generator_description_blk_path, gen_info);
    Block &gen_mesh_parts = *gen_info.get_block("mesh_parts");

    std::vector<float> params = {1.000, 1.853, 0.040, 0.322, 0.005, 2.000, 5.000, 1.000, 5.000, 650.000, 1.000, 5.000, 341.000, 1.000, 0.080, 0.080, 0.100, 0.008, 0.410, 0.000, 0.024, 0.315, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 2.174, 1.000, 3.000, 0.000, 0.600, 0.486, 0.600, 1.484, 0.400, 0.500, 0.015, 0.050, 1.000, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 1.721, 0.064, 0.478, 0.737, 0.150, 0.100, 0.200, 0.416, 1.000, -0.160, -0.049, 0.954, 0.072, 0.523, 0.008, 0.000, 0.500, 10.000, 1.000, 100.000, 0.100, 0.300};
    dgen::DFModel res;
    dgen::dgen_test("buildings_2", params, res, false, dgen::ModelQuality(false, 0));
    dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    
    auto res_bbox = dgen::get_bbox(res.first);
    logerr("model bbox 3 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    dgen::normalize_model(res.first);
    //for (int i=0;i<model.size();i++)
    //  logerr("%d %f",i, model[i]);
    res_bbox = dgen::get_bbox(res.first);
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
    dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    dgen::normalize_model(res.first);

    std::vector<float> params2 = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.07, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res2;
    dgen::dgen_test("dishes", params2, res2, false, dgen::ModelQuality(false, 2));
    dgen::transform(res2.first, glm::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    dgen::normalize_model(res2.first);

    float iou = iou3d(res2.first, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/32);
    logerr("IoU %f", iou);
}

void cup_1_iou()
{
    std::vector<float> params = {3.625, 3.821, 4.087, 4.098, 4.198, 4.241, 4.311, 4.343, 4.347, 0.748, 1.000, 0.2, 0.590, 0.214, 0.185, 0.226, 0.255, 0.278, 0.290, 0.293, 0.295, 0.290, 0.289, 0.282, 0.282, 0.284, 0.285, 0.289, 0.288, 0.295, 0.295, 0.275, 0.358, 1.155, 1.350, 0.969, 0.930, 0.799, 0.733, 0.710, 0.681, 0.672, 0.675, 0.672, 0.679, 0.687, 0.691, 0.717, 0.727, 0.768, 0.871, 1.537, 1.942, -0.152, 0.543, 0.239, -0.026, 2.957, 0.001, 15.396, 1.175, 667.281, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res;
    dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
    dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    dgen::normalize_model(res.first);

    auto model = dgen::load_obj("prezentations/spring_23_medialab/cup_model_1/cup_1.obj");
    dgen::normalize_model(model);

    float iou = iou3d(model, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/64);
    logerr("Cup 1 IoU %f", iou);
}

void cup_4_iou()
{
    std::vector<float> params = {2.681, 3.268, 3.711, 3.896, 4.045, 4.068, 4.218, 4.540, 4.965, 0.906, 1.000, 0.043, 0.519, 0.101, 0.141, 0.183, 0.257, 0.283, 0.280, 0.274, 0.263, 0.252, 0.240, 0.223, 0.200, 0.171, 0.162, 0.155, 0.157, 0.166, 0.179, 0.207, 0.291, 1.091, 0.578, 1.950, 0.871, 0.805, 0.640, 0.674, 0.688, 0.720, 0.783, 0.810, 0.791, 0.582, 0.700, 0.620, 0.500, 0.552, 0.783, 1.271, 1.629, 0.082, 0.572, 0.402, 0.162, -0.090, 0.003, 0.590, 12.159, 666.696, 1.000, 79.108, 0.200, 0.250};
    dgen::DFModel res;
    dgen::dgen_test("dishes", params, res, false, dgen::ModelQuality(false, 2));
    dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI, glm::vec3(0,1,0)));
    dgen::normalize_model(res.first);

    auto model = dgen::load_obj("prezentations/spring_23_medialab/cup_model_4/cup_4.obj");
    dgen::normalize_model(model);

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
    dgen::transform(res.first, glm::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,1,0)));
    auto res_bbox = dgen::get_bbox(res.first);
    dgen::normalize_model(res.first);

    auto model = dgen::load_obj("prezentations/spring_23_medialab/test_building_2/original/original.obj");
    dgen::normalize_model(model);

    float iou = iou3d(model, res.first, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 1.0/32);
    logerr("Building 2 IoU %f", iou);
}


void render_NGP_turntable(MitsubaInterface &mi, CameraSettings &camera, std::string mygen_name)
{
    auto model4 = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t4_simple.obj");
    auto model16 = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t16_simple.obj");
    auto model64 = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t64_simple.obj");
    // dgen::shift(model, glm::vec3(0,0,0.5));
    auto bbox = dgen::get_bbox(model4);
    auto ref_model = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/ref.obj");
    dgen::transform(ref_model, glm::rotate(glm::mat4(1.0f), -PI / 2, glm::vec3(0, 1, 0)));
    auto ref_bbox = dgen::get_bbox(ref_model);
    logerr("model bbox 4 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    logerr("ref   bbox 4 (%f %f %f)(%f %f %f)", ref_bbox.min_pos.x, ref_bbox.min_pos.y, ref_bbox.min_pos.z,
           ref_bbox.max_pos.x, ref_bbox.max_pos.y, ref_bbox.max_pos.z);

    glm::vec3 sizes = ref_bbox.max_pos - ref_bbox.min_pos;
    float max_size = MAX(sizes.x, MAX(sizes.y, sizes.z));
    max_size = MAX(1e-6, max_size);
    dgen::shift(ref_model, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(ref_model, glm::vec3(1 / max_size));
    dgen::shift(model4, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(model4, glm::vec3(1 / max_size));
    dgen::shift(model16, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(model16, glm::vec3(1 / max_size));
    dgen::shift(model64, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(model64, glm::vec3(1 / max_size));
    // dgen::normalize_model(model);
    // bbox = dgen::get_bbox(model);

    dgen::DFModel df_model4 = {model4, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, df_model4, "../../prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "prezentations/spring_23_medialab/" + mygen_name + "/NGP4_turntable", 1024, 64, 16, 3, 0, camera);
    dgen::DFModel df_model16 = {model16, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, df_model16, "../../prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "prezentations/spring_23_medialab/" + mygen_name + "/NGP16_turntable", 1024, 64, 16, 3, 0, camera);
    dgen::DFModel df_model64 = {model64, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, df_model64, "../../prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "prezentations/spring_23_medialab/" + mygen_name + "/NGP64_turntable", 1024, 64, 16, 3, 0, camera);
    dgen::DFModel ref_df_model = {ref_model, dgen::PartOffsets{{"main_part", 0}}};
    render_normalized(mi, ref_df_model, "../../prezentations/spring_23_medialab/cup_model_1/tex_inv.png",
                      "prezentations/spring_23_medialab/" + mygen_name + "/NGP_turntable_ref", 1024, 64, 16, 3, 0, camera);
}

void calc_NGP_3D_IoU(MitsubaInterface &mi, CameraSettings &camera, std::string mygen_name)
{
    auto model4 = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t4_simple.obj");
    auto model16 = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t16_simple.obj");
    auto model64 = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/t64_simple.obj");
    // dgen::shift(model, glm::vec3(0,0,0.5));
    auto bbox = dgen::get_bbox(model4);
    auto ref_model = dgen::load_obj("prezentations/spring_23_medialab/" + mygen_name + "/NGP_meshes/ref.obj");
    dgen::transform(ref_model, glm::rotate(glm::mat4(1.0f), -PI / 2, glm::vec3(0, 1, 0)));
    auto ref_bbox = dgen::get_bbox(ref_model);
    logerr("model bbox 4 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    logerr("ref   bbox 4 (%f %f %f)(%f %f %f)", ref_bbox.min_pos.x, ref_bbox.min_pos.y, ref_bbox.min_pos.z,
           ref_bbox.max_pos.x, ref_bbox.max_pos.y, ref_bbox.max_pos.z);

    glm::vec3 sizes = ref_bbox.max_pos - ref_bbox.min_pos;
    float max_size = MAX(sizes.x, MAX(sizes.y, sizes.z));
    max_size = MAX(1e-6, max_size);
    dgen::shift(ref_model, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(ref_model, glm::vec3(1 / max_size));
    dgen::shift(model4, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(model4, glm::vec3(1 / max_size));
    dgen::shift(model16, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(model16, glm::vec3(1 / max_size));
    dgen::shift(model64, -0.5f * (ref_bbox.max_pos + ref_bbox.min_pos));
    dgen::scale(model64, glm::vec3(1 / max_size));

    float iou1 = iou3d(ref_model, model4, -1, -1, -1, 1, 1, 1, 1.0/16);
    float iou2 = iou3d(ref_model, model16, -1, -1, -1, 1, 1, 1, 1.0/16);
    float iou3 = iou3d(ref_model, model64, -1, -1, -1, 1, 1, 1, 1.0/16);

    logerr("%s NGP IoU %f %f %f", mygen_name.c_str(), iou1, iou2, iou3);
}

void compare_sandbox(int argc, char **argv)
{
  MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");

  float fov_rad = 0.5;
  CameraSettings camera = MitsubaInterface::get_camera_from_scene_params({fov_rad});

  //building_2_render_mygen_turntable(mi, camera);

  //compare_and_print("test_building_2", "building_2");
  //compare_and_print("cup_model_4", "cup_4");
  //compare_and_print("cup_model_1", "cup_1");
  //building_2_render_mygen_turntable(mi, camera);
  building_2_render_reference_turntable(mi, camera);
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