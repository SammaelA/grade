#include "sandbox.h"
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
#include <boost/algorithm/string.hpp>
#include <thread>
#include <chrono>
#include <time.h>
struct BS_Grid
{
  std::vector<std::pair<int, int>> bins;
  ParameterList base_set;
};

ParameterList grid_params(std::vector<std::pair<int, int>> bins, ParameterList &base_set)
{
  ParameterList params = base_set;
  int k = 0;
  for (auto &p : params.categorialParameters)
  {
    if (!p.second.fixed())
    {
      p.second.val = p.second.possible_values[bins[k].first];
      p.second.possible_values = {p.second.val};
    }
    k++;
  }
  for (auto &p : params.ordinalParameters)
  {
    if (!p.second.fixed())
    {
      p.second.val = p.second.min_val + bins[k].first;
      p.second.min_val = p.second.val;
      p.second.max_val = p.second.val;
    }
    k++;
  }
  for (auto &p : params.continuousParameters)
  {
    if (!p.second.fixed())
    {
      float step = (p.second.max_val - p.second.min_val) / bins[k].second;
      float min = p.second.min_val;
      p.second.val = min + (bins[k].first + urand()) * step;
      p.second.min_val = min + (bins[k].first) * step;
      p.second.max_val = min + (bins[k].first + 1) * step;
    }
    k++;
  }
  return params;
}

void fill_bs_grid_bins(BS_Grid &grid, ParameterList &params, int num_bins)
{
  for (auto &p : params.categorialParameters)
  {
    if (!p.second.fixed())
      grid.bins.push_back(std::pair<int, int>(0, p.second.possible_values.size()));
    else
      grid.bins.push_back(std::pair<int, int>(0, 0));
  }
  for (auto &p : params.ordinalParameters)
  {
    if (!p.second.fixed())
      grid.bins.push_back(std::pair<int, int>(0, p.second.max_val - p.second.min_val));
    else
      grid.bins.push_back(std::pair<int, int>(0, 0));
  }
  for (auto &p : params.continuousParameters)
  {
    if (!p.second.fixed())
      grid.bins.push_back(std::pair<int, int>(0, num_bins));
    else
      grid.bins.push_back(std::pair<int, int>(0, 0));
  }
}

void bruteforce_selection(const std::function<std::vector<float>(std::vector<ParameterList> &)> &f, int num_bins,
                          int detalization_count, int detalization_depth, int num_samples, float &best_val,
                          ParameterList &bestParams)
{
  std::vector<BS_Grid> progress_bars;
  progress_bars.emplace_back();
  progress_bars[0].base_set = bestParams;
  fill_bs_grid_bins(progress_bars[0], bestParams, num_bins);

  for (int layer = 0; layer < detalization_depth; layer++)
  {
    if (layer == detalization_depth - 1)
      detalization_count = 1;

    std::vector<BS_Grid> new_progress_bars;
    for (auto &grid : progress_bars)
    {
      std::vector<std::pair<float, std::vector<std::pair<int, int>>>> cur_best;
      int i = 0;
      while (i < grid.bins.size())
      {
        float sum_metric = 0;
        std::vector<ParameterList> params;
        for (int sample = 0; sample < num_samples; sample++)
        {
          params.push_back(grid_params(grid.bins, grid.base_set));
        }
        auto vec = f(params);
        for (auto &v : vec)
        {
          sum_metric += v;
          logerr("v = %f", v);
        }
        logerr("");
        sum_metric /= num_samples;
        if (cur_best.size() < detalization_count)
        {
          cur_best.push_back(std::pair<float, std::vector<std::pair<int, int>>>(sum_metric, grid.bins));
        }
        else
        {
          float min_val = 1e9;
          int min_pos = 0;
          for (int j = 0; j < cur_best.size(); j++)
          {
            if (cur_best[j].first < min_val)
            {
              min_val = cur_best[j].first;
              min_pos = j;
            }
          }
          if (min_val < sum_metric)
            cur_best[min_pos] = std::pair<float, std::vector<std::pair<int, int>>>(sum_metric, grid.bins);
        }

        if (grid.bins[i].first < grid.bins[i].second - 1)
        {
          grid.bins[i].first++;
        }
        else
        {
          i++;
          while (i < grid.bins.size() && grid.bins[i].second <= MAX(1, grid.bins[i].first + 1))
            i++;
          if (i < grid.bins.size())
          {
            grid.bins[i].first++;
            int t = -1;
            for (int j = 0; j < i; j++)
            {
              grid.bins[j].first = 0;
              if (grid.bins[j].second > 1 && t < 0)
                t = j;
            };
            if (t >= 0)
              i = t;
          }
        }
        /*
        for (int j=0;j<grid.bins.size();j++)
        {
            if (grid.bins[j].second > 1)
                debug("%d/%d ", grid.bins[j].first, grid.bins[j].second);
        }
        debug("  i = %d\n",i);
        */
      }
      if (layer < detalization_depth - 1)
      {
        for (auto &best : cur_best)
        {
          new_progress_bars.emplace_back();
          new_progress_bars.back().base_set = grid_params(best.second, grid.base_set);
          // new_progress_bars.back().base_set.print();
          fill_bs_grid_bins(new_progress_bars.back(), new_progress_bars.back().base_set, num_bins);
        }
      }
      else
      {
        best_val = cur_best[0].first;
        bestParams = grid_params(cur_best[0].second, grid.base_set);
      }
    }
    progress_bars = std::move(new_progress_bars);
  }
}

float dot_metric(Tree &single_tree, float dst_dot)
{
  double sum_dot = 0;
  int dot_cnt = 0;

  for (auto &bh : single_tree.branchHeaps)
  {
    for (auto &b : bh->branches)
    {
      glm::vec3 dir = normalize(b.segments.front().begin - b.segments.front().end);
      for (auto &j : b.joints)
      {
        for (auto *chb : j.childBranches)
        {
          glm::vec3 ch_dir = normalize(chb->segments.front().begin - chb->segments.front().end);
          float w = pow(10, 3 - b.level);
          sum_dot += w * SQR(0.5 - CLAMP(glm::dot(dir, ch_dir), 0, 1));
          dot_cnt += w;
        }
      }
    }
  }
  float dt = sum_dot / dot_cnt;
  float metric = 1 - dt;
  logerr("dot metric %f %f", dt, metric);
  return metric;
}

float f(float x)
{
  return 2 * (sin(PI * x + 0.1) * cos(PI * x + 0.1)) - 1e-8;
}

#include "graphics_utils/silhouette.h"

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

  CameraSettings camera;
  camera.origin = glm::vec3(0, 0.5, 1.5);
  camera.target = glm::vec3(0, 0.5, 0);
  camera.up = glm::vec3(0, 1, 0);

  if (argc >= 4 && std::string(argv[2]) == "-sil_test")
  {
    Texture t = engine::textureManager->load_unnamed_tex(std::string(argv[3]));
    SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
    Texture tex = se.get_silhouette(t, 256, 256);
    textureManager.save_png(tex, "silhouette_test");
    engine::view->next_frame();
    logerr("Silhouette test completed. Saved to saves/silhouette_test.png");
    return;
  }
  else if ((argc >= 3 && std::string(argv[2]) == "-check_stability"))
  {
    std::vector<float> reference_params{4 - 1.45, 4 - 1.0, 4 - 0.65, 4 - 0.45, 4 - 0.25, 4 - 0.18, 4 - 0.1, 4 - 0.05, 4,//spline point offsets
                                        0.4,// y_scale
                                        1, //has handle variant
                                        0.05, 0.35, 0.35};//handle params
    //dgen::check_stability(dgen::create_cup, reference_params, 4);
    return;
  }
  else if (argc >= 4 && std::string(argv[2]) == "-opt_benchmark")
  {
    srand(time(NULL));
    std::string blk_name = std::string(argv[3]);
    Block b;
    load_block_from_file(blk_name, b);
    std::vector<std::string> reference_images;
    for (int i=0;i<b.size();i++)
    {
      if (b.get_name(i) == "reference")
      {
        reference_images.push_back(b.get_string(i));
      }
    }
    if (reference_images.empty())
    {
      logerr("no reference images found in blk");
      return;
    }
    float av_loss = 0;
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    int count = MIN(b.get_int("count",1000), reference_images.size());
    for (int i=0;i<count;i++)
    {
      std::string &ref = reference_images[i];
      b.set_string("reference_path", ref);
      std::vector<std::string> split_res;
      boost::algorithm::split(split_res, ref, boost::is_any_of("/"));
      b.set_string("saved_result_path", "saves/"+split_res.back()+"_result.png");
      b.set_string("saved_textured_path", "saves/"+split_res.back()+"_result_textured.png");
      av_loss += dopt::image_based_optimization(b, mi);
    }
    av_loss /= count;
    debug("Benchmak finished. %d images tested\n", count);
    debug("Average loss: %.4f\n", av_loss);
    return;
  }
  else if (argc >=3 && std::string(argv[2]) == "-test_gen")
  {
    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    std::vector<float> res;
    dgen::dgen_test("dishes", params, res);
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
    mi.render_model_to_file(res, "saves/test_result.png", dgen::ModelLayout(), camera);
  }
  else if (argc >=3 && std::string(argv[2]) == "-test_gen_with_camera")
  {
    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    std::vector<float> res;
    dgen::dgen_test("dishes", params, res, true);
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
    mi.render_model_to_file(res, "saves/test_result.png", dgen::ModelLayout(), camera);
  }
  else if (argc >=3 && std::string(argv[2]) == "-test_tex")
  {
    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    std::vector<float> res;
    dgen::dgen_test("dishes", params, res);
    Model *m = new Model();
    visualizer::simple_mesh_to_model_332(res, m);
    m->update();

    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "porcelain_01.png"));
    mi.render_model_to_file(res, "saves/tex_colored.png", dgen::ModelLayout(), camera);
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE));
    mi.render_model_to_file(res, "saves/tex_sihouette.png", dgen::ModelLayout(), camera);

    engine::view->next_frame();
    Texture photo = textureManager.load_unnamed_tex("saves/tex_colored.png");
    Texture mask = textureManager.load_unnamed_tex("saves/tex_sihouette.png");
    ModelTex mt;
    Texture res_tex = mt.getTexbyUV(mask, *m, photo, 3, camera);
    textureManager.save_png(res_tex, "reconstructed_tex");

    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "../../saves/reconstructed_tex.png"));
    mi.render_model_to_file(res, "saves/tex_reconstructed.png", dgen::ModelLayout(), camera);

    std::vector<ModelTex::tex_data> data = {{0, 0, 1, 0.75, 4, -1}, {0, 0.75, 1, 1, 1, 1}};

    Texture comp = mt.symTexComplement(res_tex, data);
    textureManager.save_png(comp, "complement_tex");
    engine::view->next_frame();

    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "../../saves/complement_tex.png"));
  } 
  else if (argc >=3 && std::string(argv[2]) == "-test_tex_diff")
  {
    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    std::vector<float> res;
    dgen::dgen_test("dishes", params, res);
    Model *m = new Model();
    visualizer::simple_mesh_to_model_332(res, m);
    m->update();

    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "porcelain_01.png"));
    mi.render_model_to_file(res, "saves/tex_colored.png", dgen::ModelLayout(), camera);
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE));
    mi.render_model_to_file(res, "saves/tex_sihouette.png", dgen::ModelLayout(), camera);

    engine::view->next_frame();
    Texture photo = textureManager.load_unnamed_tex("saves/tex_colored.png");
    Texture mask = textureManager.load_unnamed_tex("saves/tex_sihouette.png");
    ModelTex mt;
    Texture res_tex = mt.getTexbyUV(mask, *m, photo, 3, camera);
    textureManager.save_png(res_tex, "reconstructed_tex");
    engine::view->next_frame();

    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "../../saves/reconstructed_tex.png"));
    mi.render_model_to_file(res, "saves/tex_reconstructed.png", dgen::ModelLayout(), camera);

    mi.init_optimization_with_tex("saves/tex_colored.png", "../../saves/reconstructed_tex.png", MitsubaInterface::LossFunction::LOSS_MSE, 1 << 16, 
                                  dgen::ModelLayout(0, 3, 6, 8, 8), 
                                  MitsubaInterface::RenderSettings(512, 512, 64, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST),
                                  true);

    for (int i=0;i<50;i++)
    {
      float loss = mi.render_and_compare(res, camera);
      logerr("%d loss = %f",i, loss);
    }
  } 
  else if (argc >= 3 && std::string(argv[2]) == "-opt_test")
  {
    opt::optimization_benchmark();
  }
  else if (argc >= 3 && std::string(argv[2]) == "-depth_test")
  {
    std::vector<float> params;
    for (int i=3;i<argc;i++)
    {
      params.push_back(std::stof(std::string(argv[i])));
    }
    std::vector<float> res;
    MitsubaInterface mi("scripts", "mitsuba_optimization_embedded");
    DepthLossCalculator dlc;
    
    dgen::dgen_test("dishes", params, res);
    Model *m = new Model();
    visualizer::simple_mesh_to_model_332(res, m);
    m->update();
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "porcelain_01.png"));
    mi.render_model_to_file(res, "saves/depth_test_tex_colored_1.png", dgen::ModelLayout(), camera);
    Texture d1 = dlc.get_depth(*m, camera, 256, 256);
    textureManager.save_png(d1, "depth_test_depth_1");

    params[0] += 1;
    params[1] += 1;
    res.clear();
    dgen::dgen_test("dishes", params, res);
    m = new Model();
    visualizer::simple_mesh_to_model_332(res, m);
    m->update();
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 256, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "porcelain_01.png"));
    mi.render_model_to_file(res, "saves/depth_test_tex_colored_2.png", dgen::ModelLayout(), camera);
    Texture d2 = dlc.get_depth(*m, camera, 256, 256);
    textureManager.save_png(d2, "depth_test_depth_2");
  
    float diff = dlc.get_loss(*m, d1, camera);
    logerr("depth diff %f", diff);

    engine::view->next_frame();
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
  return;
  //std::vector<float> model;
  //dgen::dgen_test(model);

  int quantiles[101];
  for (int i = 0; i < 101; i++)
    quantiles[i] = 0;
  int cnt_x = 100;
  int cnt_y = 10000;
  for (int j = 0; j < cnt_x; j++)
  {
    float x = 0.35;
    for (int i = 0; i < cnt_y; i++)
    {
      int pers = 100 * (0.5 * (x + 1));
      if (pers < 0 || pers >= 100)
      {
        // logerr("%.4f", x);
        pers = 100;
      }
      quantiles[pers]++;
      x = f(x);
    }
  }
  for (int i = 0; i < 100; i++)
    logerr("%.3f", (100 * (double)quantiles[i]) / (cnt_x * cnt_y));
  return;

  /*
  for (int i=0;i<2500;i++)
  {
      WeberPennParametersNative param;
      WeberPennGenerator::Tree t;
      t.init(param, true);
      t.make();
      t.clear();
  }
  return;
  */
  /*
   {
       metainfoManager->reload_all();
       TreeTypeData type = metainfoManager->get_tree_type("apple");
       scene.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(100, 100), 10);
       scene.heightmap->fill_const(0);

       ParameterList par;
       type.get_params()->write_parameter_list(par);

       Block b;
       load_block_from_file("weber_penn_gen_param_borders.blk", b);
       par.load_borders_from_blk(b);
       par.print();
       return;
       GroveGenerationData tree_ggd;
       tree_ggd.trees_count = 1;
       tree_ggd.types = {type};
       tree_ggd.name = "single_tree";
       tree_ggd.task = GenerationTask::IMPOSTORS;
       tree_ggd.impostor_generation_params.slices_n = 8;
       tree_ggd.impostor_generation_params.quality = 128;
       tree_ggd.impostor_generation_params.monochrome = true;
       tree_ggd.impostor_generation_params.normals_needed = false;
       tree_ggd.impostor_generation_params.leaf_opacity = 0.95;

       NeuralEstimator NE;
       NE.prepare_dataset(par, tree_ggd, "saves/NE_dataset", 256, 4, 512, 64);
       return;
   }
   */
  /*
      int cnt = 100;
      int _a;
      ParameterList par;
      type.params->write_parameter_list(par);
      std::vector<ParameterList> params;
      for (int i=0;i<cnt;i++)
      {
          params.push_back(par);
      }
      float imp_size = 128;
      GroveGenerationData tree_ggd;
      tree_ggd.trees_count = 1;
      tree_ggd.types = {type};
      tree_ggd.name = "single_tree";
      tree_ggd.task = GenerationTask::IMPOSTORS;
      tree_ggd.impostor_generation_params.slices_n = 8;
      tree_ggd.impostor_generation_params.quality = imp_size;
      tree_ggd.impostor_generation_params.monochrome = true;
      tree_ggd.impostor_generation_params.normals_needed = false;
      tree_ggd.impostor_generation_params.leaf_opacity = 0.33;

      ReferenceTree ref_tree;
      ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(cnt, 8, false);
      generate_for_par_selection(params, imp_sim, tree_ggd, scene.heightmap, ref_tree, _a, nullptr);
  */
  Block b, ref_info;

  metainfoManager->reload_all();
  load_block_from_file("parameter_selection_settings.blk", b);
  load_block_from_file("parameter_selection_reference.blk", ref_info);
  std::string add_str = "";
  if (argc == 3 && std::string(argv[2]) != "-render")
  {
    logerr("argv %s", argv[2]);
    add_str = argv[2];
  }
  if (add_str != "")
  {
    Block add_ref;
    load_block_from_string(add_str, add_ref);
    ref_info.add_detalization(add_ref);
  }
  ParameterSelector sel;
  auto res = sel.parameter_selection(ref_info, b, scene);
  /*
     metainfoManager->reload_all();
     scene.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(100, 100), 10);
     scene.heightmap->fill_const(0);
     TreeTypeData type = metainfoManager->get_tree_type("apple");
     float imp_size = 512;
     GroveGenerationData tree_ggd;
     tree_ggd.trees_count = 1;
     tree_ggd.types = {type};
     tree_ggd.name = "single_tree";
     tree_ggd.task = GenerationTask::IMPOSTORS;
     tree_ggd.impostor_generation_params.slices_n = 1;
     tree_ggd.impostor_generation_params.quality = imp_size;
     tree_ggd.impostor_generation_params.monochrome = true;
     tree_ggd.impostor_generation_params.normals_needed = false;
     tree_ggd.impostor_generation_params.leaf_opacity = 0.99;
     tree_ggd.impostor_generation_params.leaf_scale = 2.5;
     srand(2);
     ReferenceTree ref_tree;

     //create reference tree

         glm::vec3 pos = glm::vec3(0,0,0);
         glm::vec3 sz = type.get_params()->get_tree_max_size();
         LightVoxelsCube *ref_voxels = new LightVoxelsCube(pos + glm::vec3(0, sz.y - 10, 0),sz,
                                                           0.625f * type.get_params()->get_scale_factor(), 1.0f, 1, 2);

     int cnt = 1;
     std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
     for (int i=0;i<cnt*cnt;i++)
     {
         pos = glm::vec3(100*(i / cnt),0,100*(i % cnt));
         //logerr("GEN %d",i);
         AbstractTreeGenerator *gen = get_generator(type.generator_name);
         ref_voxels->fill(0);
         ref_voxels->relocate(pos + glm::vec3(0, sz.y - 10, 0));
         GrovePacker packer;
         Tree single_tree;
         tree_ggd.trees_count = 1;
         gen->plant_tree(pos, &(tree_ggd.types[0]));
         while (gen->iterate(*ref_voxels))
         {
         }
         gen->finalize_generation(&single_tree, *ref_voxels);
         packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
     }
     std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
     float time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
     logerr("took %.3f seconds, %.1f ms/tree", time/1000, time/(cnt*cnt));
     delete ref_voxels;
     engine::textureManager->save_png(scene.grove.impostors[1].atlas.tex(0), "original_atlas");

     Texture tex1 = engine::textureManager->load_unnamed_tex(image::base_img_path + "24_A_mono.png");
     Texture tex2 = engine::textureManager->load_unnamed_tex(image::base_img_path + "24_B_mono.png");
     //Texture tex3 = engine::textureManager->create_unnamed(tex1.get_W(), tex1.get_H());
     TextureAtlas tmp_atlas = TextureAtlas(2*tex1.get_W(), 2*tex1.get_H(), 1);
     tmp_atlas.set_grid(2*tex1.get_W(), 2*tex1.get_H());
     int id = tmp_atlas.add_tex();
     tmp_atlas.target_slice(id, 0);
     PostFx pixel_dist = PostFx("pixel_difference.fs");
     pixel_dist.use();
     pixel_dist.get_shader().texture("tex1", tex1);
     pixel_dist.get_shader().texture("tex2", tex2);
     pixel_dist.render();
     engine::textureManager->save_png(tmp_atlas.tex(0), "cmp");
 */
  /*
     LightVoxelsCube test = LightVoxelsCube(glm::vec3(0,0,0), glm::vec3(200,200,200),1.0f,1.0f,1,2);
     LightVoxelsCube ref = LightVoxelsCube(glm::vec3(0,0,0), glm::vec3(200,200,200),1.0f,1.0f,1,2);
     int cnt = 50000;
     srand(0);
     std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
     long r1, r2, r3;
     r1 = 0;
     r2 = 0;
     r3 = 0;
     for (int i=0;i<cnt;i++)
     {
         r1 = (3*r1 + 17) % 400;
         r2 = (5*r2 + 19) % 400;
         r3 = (7*r3 + 23) % 400;
         glm::vec3 pos = glm::vec3(r1 - 200, r2 - 200, r3 - 200);
         ref.set_occluder_pyramid2(pos, 1, 2, 7);
     }
     std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
     float time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
     logerr("Reference took %.3f seconds, %d op/s", time/1000, (int)(cnt*1000.0f/time));

     srand(0);
     t1 = std::chrono::steady_clock::now();
     r1 = 0;
     r2 = 0;
     r3 = 0;
     for (int i=0;i<cnt;i++)
     {
         r1 = (3*r1 + 17) % 400;
         r2 = (5*r2 + 19) % 400;
         r3 = (7*r3 + 23) % 400;
         glm::vec3 pos = glm::vec3(r1 - 200, r2 - 200, r3 - 200);
         test.set_occluder_pyramid_fast(pos, 1, 7);
     }
     t2 = std::chrono::steady_clock::now();
     time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
     logerr("Took %.3f seconds, %d op/s", time/1000, (int)(cnt*1000.0f/time));

     float *ref_data;
     glm::ivec3 ref_sizes;
     float *test_data;
     glm::ivec3 test_sizes;
     ref.get_data(&ref_data, ref_sizes);
     test.get_data(&test_data, test_sizes);
     if (ref_sizes != test_sizes)
     {
         logerr("AAA %d %d %d -- %d %d %d",ref_sizes.x, ref_sizes.y, ref_sizes.z, test_sizes.x, test_sizes.y, test_sizes.z);
     }
     else if (ref.get_size_cnt() == test.get_size_cnt())
     {
         int vox_cnt = ref.get_size_cnt();
         int wrong_voxels = 0;
         for (int i=0;i<vox_cnt;i++)
         {
             if (abs(ref_data[i] - test_data[i]) >= 0.001)
              wrong_voxels++;
         }
         logerr("Wrong voxels %d from %d", wrong_voxels,vox_cnt);
     }
  */
}