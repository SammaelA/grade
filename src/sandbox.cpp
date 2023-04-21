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

void defaultSignalHandler(int signum)
{
  logerr("Interrupt signal received. Closing the program");
  exit(signum);  
}

void sandbox_main(int argc, char **argv, Scene *scene)
{
  // we don't init engine in sandbox, so need to init textures manager
  View view;
  view.lineWidth = 1.0f;
  view.init("Sandbox", 256, 256);
  engine::view = &view;

  Block textures_list;
  TextureManager textureManager = TextureManager("./resources/textures/", textures_list);
  engine::textureManager = &textureManager;
  signal(SIGINT, defaultSignalHandler);  
  signal(SIGTSTP, defaultSignalHandler); 

  CameraSettings camera;
  float h1 = 1.5;
  camera.fov_rad = 0.5;
  float h2 = h1 * tan((PI / 3) / 2) / tan(camera.fov_rad / 2);
  camera.origin = glm::vec3(0, 0.5, h2);
  camera.target = glm::vec3(0, 0.5, 0);
  camera.up = glm::vec3(0, 1, 0);

  std::vector<float> default_scene_params = {0.133, 0.543, 0.238, 0.088, -0.353, 0.023, 0.000, 0.500, 10.000, 1.000, 100.000, 0.1};
  if (argc >= 4 && std::string(argv[2]) == "-sil_test")
  {
    Texture t = textureManager.load_unnamed_tex(std::string(argv[3]));
    SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
    Texture tex = se.get_silhouette(t, 256, 256);
    textureManager.save_png(tex, "silhouette_test");
    engine::view->next_frame();
    logerr("Silhouette test completed. Saved to saves/silhouette_test.png");
    return;
  }
  else if (argc >= 4 && std::string(argv[2]) == "-opt_benchmark")
  {
    srand(time(NULL));
    std::string blk_name = std::string(argv[3]);
    Block b;
    load_block_from_file(blk_name, b);
    std::vector<std::string> reference_images;
    std::vector<std::string> reference_masks;
    if (argc == 4)
    {
      for (int i=0;i<b.size();i++)
      {
        if (b.get_name(i) == "reference")
        {
          reference_images.push_back(b.get_string(i));
          reference_masks.push_back("");
        }
      }
      if (reference_images.empty())
      {
        logerr("no reference images found in blk");
        return;
      }
    }
    else if (argc == 5)
    {
      reference_images.push_back(argv[4]);
      reference_masks.push_back("");
    }
    else if (argc == 6)
    {
      reference_images.push_back(argv[4]);
      reference_masks.push_back(argv[5]);
    }
    float av_loss = 0;
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    int count = MIN(b.get_int("count",1000), reference_images.size());
    for (int i=0;i<count;i++)
    {
      std::string &ref = reference_images[i];
      std::vector<std::string> split_res ={""};
      boost::algorithm::split(split_res, ref, boost::is_any_of("/"));
      b.set_string("saved_result_path", "saves/"+split_res.back()+"_result.png");
      b.set_string("saved_textured_path", "saves/"+split_res.back()+"_result_textured.png");

      Block ref_block;
      ref_block.set_string("textured", ref);
      ref_block.set_string("mask", reference_masks[i]);
      b.add_block("references", &ref_block);
      av_loss += dopt::image_based_optimization(b, mi);
    }
    av_loss /= count;
    debug("Benchmak finished. %d images tested\n", count);
    debug("Average loss: %.4f\n", av_loss);
    return;
  }
  else if (argc >= 4 && std::string(argv[2]) == "-opt")
  {
    std::vector<std::string> arg_names = {"program", "-sandbox", "-opt", "blk_name", "textured", "mask", "mask_1", "mask_2", "mask_3", "mask_4"};
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    std::string blk_name = std::string(argv[3]);
    Block b;
    Block ref_block;
    load_block_from_file(blk_name, b);
    for (int i = 4; i< MIN(arg_names.size(), argc); i++)
    {
      ref_block.set_string(arg_names[i], argv[i]);
    }
    b.add_block("references", &ref_block);
    dopt::image_based_optimization(b, mi);
  }
  else if (argc >=3 && std::string(argv[2]) == "-test_gen")
  {
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("dishes").generator_description_blk_path, gen_info);
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

    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    if (params.empty())
      params = {3.108, 3.521, 3.890, 4.132, 4.354, 4.457, 4.624, 4.696, 4.745, 1.139, 1.000, 0.041, 0.558, 
                0.137, 0.18, 0.25, 0.35, 0.409, 0.439, 0.465, 0.450, 0.413, 0.358, 0.315, 0.287, 0.264, 0.250, 0.244, 0.241, 0.247, 0.256, 0.240, 0.330, 
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    dgen::DFModel res;
    dgen::dgen_test("dishes", params, res);
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 2048, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_DEMO),
                               model_info);
    mi.render_model_to_file(res, "saves/test_result.png", camera, default_scene_params);
  }
  else if (argc >=3 && std::string(argv[2]) == "-test_gen_buildings")
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

    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    if (params.empty())
      params = {1.029, 
              1.926, 
              0.040, 
              0.100, 0.300, 
              3, 6, 1, 
              5, 650, 0.5, 
              2.000, 103.292, 0.749, 
              0.080, 0.080, 0.095, 
              0.001, 0.300, 0.000, 0.000, 0.241, 
              2.797, 1.443, 2.000, 0.790, 0.400, 0.400, 3.8, 
              1.761, 3.000, 0.000, 0.581, 0.667, 0.4, 1.6,
              0.400, 0.500, 0.018, 0.045, 1.000, 
              2.719, 2.000, 2.000, 0.600, 0.400, 0.400, 3.2, 
              0.550, 0.900, 0.906, 0.150, 
              0.116, 0.150, 
              0.251, 0.200, 
              0.482, 0.353, 2.745, 0.053, 2.831, -0.004, 0.000, 0.500, 10.000, 64.286, 100.000, 0.908, 0.650};
    params = {1.000, 1.853, 0.040, 0.322, 0.005, 2.000, 5.000, 1.000, 5.000, 650.000, 1.000, 5.000, 341.000, 1.000, 0.080, 0.080, 0.100, 0.008, 0.410, 0.000, 0.024, 0.315, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 2.174, 1.000, 3.000, 0.000, 0.600, 0.486, 0.600, 1.484, 0.400, 0.500, 0.015, 0.050, 1.000, 1.000, 2.000, 2.000, 0.600, 0.400, 0.600, 1.721, 0.064, 0.478, 0.737, 0.150, 0.100, 0.200, 0.416, 1.000, -0.160, -0.049, 0.954, 0.072, 0.523, 0.008, 0.000, 0.500, 10.000, 1.000, 100.000, 0.100, 0.300};
    dgen::DFModel res;
    dgen::dgen_test("buildings_2", params, res, false, dgen::ModelQuality(false, 3));
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(1024, 1024, 1024, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_DEMO),
                               model_info);
    std::vector<float> scene_params = {-0.160, 0.05, 0.954, 0.00, 0.523, 0.008, 0.000, 0.500, 10.000, 1.000, 100.000, 0.100, 0.300};
    mi.render_model_to_file(res, "saves/building_5.png", camera, scene_params);
  }
  else if (argc >=3 && std::string(argv[2]) == "-test_gen_buildings_multi")
  {
    Block gen_info;
    load_block_from_file(dgen::get_generator_by_name("buildings").generator_description_blk_path, gen_info);
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

    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    if (params.empty())
      params = {3, 3, 3, 3, 80,  3, 3, 3, 3, 25,   3, 3, 3, 3, 10,   1, 20};
    dgen::DFModel res;
    dgen::dgen_test("buildings", params, res);
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME),
                               model_info);
    mi.render_model_to_file(res, "saves/test_result.png", camera, default_scene_params);
    
    std::vector<float> no_transform_scene_params = {0,0,0,0,0,0, 100, 1000, 100, 100, 100, 0.01};
    mi.render_multicam_demo(MitsubaInterface::RenderSettings(512, 512, 32, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_CONST),
                            model_info, res, "saves/multicam_tex.png", no_transform_scene_params, camera, 3, 2);
    engine::view->next_frame();
  }
  else if (argc >= 3 && std::string(argv[2]) == "-opt_test")
  {
    opt::optimization_benchmark();
  }
  else if (argc >= 3 && std::string(argv[2]) == "-test_denoising")
  {
      Texture res_optimized = textureManager.load_unnamed_tex("saves/reconstructed_tex_raw.png");
      Texture mask_tex = textureManager.load_unnamed_tex("saves/reconstructed_mask.png");
      std::vector<ModelTex::tex_data> data = {{0, 0, 1, 0.375, 3, 1}, {0, 0.375, 1, 0.75, 3, 1}, {0, 0.75, 1, 1, 1, 10}};
    engine::view->next_frame();
      ModelTex mt;
      Texture comp = mt.symTexComplement(res_optimized, mask_tex, data);

      Texture res = BilateralFilter::perform(comp, 4, 0.5);
      Texture sharped = UnsharpMasking::perform(res, 1, 0.2);

      textureManager.save_png(res_optimized, "reconstructed_tex_raw_1");
      textureManager.save_png(sharped, "reconstructed_tex_complemented_1");
      textureManager.save_png(sharped, "reconstructed_tex_denoised_1");
    engine::view->next_frame();
  }
  else if (argc >= 3 && std::string(argv[2]) == "-voxelization_test")
  {
    voxelization::software_render_test_3d();
  }
  else if (argc >= 3 && std::string(argv[2]) == "-compare")
  {
    auto model = dgen::load_obj("/home/sammael/grade/resources/mitsuba_data/meshes/building/house.obj");
    auto bbox = dgen::get_bbox(model);
    logerr("model bbox 1 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    dgen::normalize_model(model);
    //for (int i=0;i<model.size();i++)
    //  logerr("%d %f",i, model[i]);
    bbox = dgen::get_bbox(model);
    logerr("model bbox 2 (%f %f %f)(%f %f %f)", bbox.min_pos.x, bbox.min_pos.y, bbox.min_pos.z, bbox.max_pos.x, bbox.max_pos.y, bbox.max_pos.z);
    dgen::save_obj("/home/sammael/grade/resources/mitsuba_data/meshes/building/house_norm.obj", model);
    //Differentiable-SDF-renderer project
    //To create this images 
    //python3.8 render_turntable.py building_2 --optconfig diffuse-6
    //Put constant emitter to building_2.xml <default name="emitter_scene" value="emitters/constant.xml"/>
    compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-6/warp/turntable",
                                  "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-6/warp/reference/turntable",
                                  64);
    compare_utils::turntable_loss("/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-12/warp/turntable",
                                  "/home/sammael/references_grade/differentiable-sdf-rendering/outputs/building_2/diffuse-12/warp/reference/turntable",
                                  64);
  }
  else
  {
    logerr("unknown sandbox command");
    logerr("./main -sandbox -h -- print help");
    logerr("./main -sandbox -opt -- optimization");
    logerr("./main -sandbox -sil_test <filename> -- silhouette test of file in argv[3]. Save to saves/silhouette_test.png");
    logerr("./main -sandbox -check_stability -- checks stability of diff procedural generator (dgen::create_cup)");
    logerr("./main -sandbox -opt_benchmark <blk> -- performs optimization with different references and settings (all set in blk)");
    logerr("./main -sandbox -test_gen <param> creates model with giver parameters and renders it with mitsuba");
    logerr("./main -sandbox -test_tex <param> tests texture reconstruction on a synthetic model");
    logerr("./main -sandbox -test_tex <param> tests texture reconstruction with mitsuba fine-tuning on a synthetic model");
    return;
  }
}