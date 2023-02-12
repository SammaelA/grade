#include "diff_optimization.h"
#include "differentiable_generators.h"
#include "mitsuba_python_interaction.h"
#include "common_utils/distribution.h"
#include <cppad/cppad.hpp>
#include "common_utils/utility.h"
#include <functional>
#include <chrono>
#include <algorithm>
#include "common_utils/blk.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/silhouette.h"
#include "save_utils/csv.h"
#include "graphics_utils/model_texture_creator.h"
#include "graphics_utils/modeling.h"
#include "common_utils/optimization/optimization.h"
#include "diff_generators/depth_extract_compare.h"
#include "graphics_utils/resize_image.h"

namespace dopt
{
  class UShortVecComparator
  {
  public:
    bool operator()(const std::vector<unsigned short> &v1, const std::vector<unsigned short> &v2) const
    {
      for (int i = 0; i < MIN(v1.size(), v2.size()); i++)
      {
        if (v1[i] < v2[i])
          return true;
        else if (v1[i] > v2[i])
          return false;
      }
      return false;
    }
  };

  class DiffFunctionEvaluator
  {
  public:
    ~DiffFunctionEvaluator()
    {
      for (auto f : functions)
      {
        if (f)
          delete f;
      }
    }
    void init(dgen::generator_func _model_creator, std::vector<unsigned short> &variant_positions)
    {
      model_creator = _model_creator;
      variant_params_positions = variant_positions;
    }
    std::vector<float> get(const std::vector<float> &params, dgen::ModelQuality mq = dgen::ModelQuality())
    {
      return functions[find_or_add(params, mq)]->Forward(0, params); 
    }
    std::vector<float> get_jac(const std::vector<float> &params, dgen::ModelQuality mq = dgen::ModelQuality())
    {
      return functions[find_or_add(params, mq)]->Jacobian(params); 
    }
    std::vector<float> get_transformed(const std::vector<float> &params, const std::vector<float> &camera_params, 
                                       dgen::ModelQuality mq = dgen::ModelQuality())
    {
      std::vector<float> f_model = functions[find_or_add(params, mq)]->Forward(0, params); 
      dgen::transform_by_scene_parameters(camera_params, f_model);  
      return f_model;
    }
  private:
    int find_or_add(const std::vector<float> &params, dgen::ModelQuality mq)
    {
      auto vs = get_variant_set(params);
      vs.push_back((unsigned short)mq.create_only_position);
      vs.push_back((unsigned short)mq.quality_level);
      auto it = variant_set_to_function_pos.find(vs);
      if (it == variant_set_to_function_pos.end())
      {
        debug("added new function {");
        for (auto &v : vs)
          debug("%d ", (int)v);
        debug("}\n");
        std::vector<dgen::dfloat> X(params.size());
        for (int i=0;i<params.size();i++)
          X[i] = params[i];
        std::vector<dgen::dfloat> Y;
        CppAD::Independent(X);
        model_creator(X, Y, mq);
        CppAD::ADFun<float> *f = new CppAD::ADFun<float>(X, Y); 
        functions.push_back(f);
        output_sizes.push_back(Y.size());
        int f_pos = functions.size()-1;
        variant_set_to_function_pos.emplace(vs, f_pos);

        return f_pos;
      }
      return it->second;
    }

    std::vector<unsigned short> get_variant_set(const std::vector<float> &params)
    {
      std::vector<unsigned short> vs;
      for (auto &pos : variant_params_positions)
      {
        vs.push_back((unsigned short)round(params[pos]));
      }
      return vs;
    }
    dgen::generator_func model_creator;
    std::map<std::vector<unsigned short>, int, UShortVecComparator> variant_set_to_function_pos;
    std::vector<CppAD::ADFun<float> *> functions;
    std::vector<int> output_sizes;//same size as functions vector
    std::vector<unsigned short> variant_params_positions;
  };

  struct OptimizationResult
  {
    std::vector<float> best_params;
    float best_err;
    int total_iters;
  };

  typedef std::vector<std::vector<std::vector<float>>> PresetsData;

  float parameters_error(const std::vector<float> &params, const std::vector<float> &ref_params,
                         const std::vector<float> &params_min, const std::vector<float> &params_max)
  {
    float err = 0;
    for (int i=0;i<params.size();i++)
    {
      err += abs(params[i] - ref_params[i]) / (params_max[i] - params_min[i]);
    }
    return err/params.size();
  }

  void load_presets_from_blk(Block &presets_blk, Block &gen_params, const std::vector<float> &init_params, 
                             PresetsData &presets_data /*output*/)
  {
    /*Presets are manually created valid sets of generator parmeters, representing different types of objects,
      that generator can create. They do not include scene parameters
      Presets blk has the following structure:
      block_1 {
        set_1_1:arr = {}
        ...
        set_1_n1:arr = {}
      }
      ...
      block_K {
        set_k_1:arr = {}
        ...
        set_k_nk:arr = {}
      }

      every possible concatenation {set_1_i1, set_2_i2, ..., set_k_ik} of these arrays (one from each block)
      is a valid preset that can be created by init_params_presets function from this data
    */
    presets_data.clear();
    for (int i = 0; i < presets_blk.size(); i++)
    {
      Block *b = presets_blk.get_block(i);
      if (b)
      {
        presets_data.emplace_back();
        auto &block =  presets_data.back();
        for (int j = 0; j < b->size(); j++)
        {
          if (b->get_type(j) == Block::ValueType::ARRAY)
          {
            block.emplace_back();
            b->get_arr(j, block.back());
          }
        }
      }
    }
  }

  float image_based_optimization(Block &settings_blk, MitsubaInterface &mi)
  {
    Block gen_params, scene_params, presets_blk;
    std::vector<float> params_min, params_max, init_params;
    std::vector<unsigned short> init_bins_count;
    std::vector<unsigned short> init_bins_positions;
    std::vector<unsigned short> variant_count;
    std::vector<unsigned short> variant_positions;
    PresetsData parameter_presets;
    dgen::GeneratorDescription generator = dgen::get_generator_by_name(settings_blk.get_string("procedural_generator"));
    load_block_from_file(generator.parameters_description_blk_path, gen_params);
    load_block_from_file(settings_blk.get_string("scene_description"), scene_params);
    load_block_from_file(generator.presets_blk_path, presets_blk);
    int gen_params_cnt = gen_params.size();
    int scene_params_cnt = scene_params.size();

    int verbose_level = settings_blk.get_int("verbose_level", 1);
    int ref_image_size = settings_blk.get_int("reference_image_size", 512);
    int sel_image_size = settings_blk.get_int("selection_image_size", 196);
    bool by_reference = settings_blk.get_bool("synthetic_reference", true);
    Block *textured_optimization = settings_blk.get_block("textured_optimization");
    std::string search_algorithm = settings_blk.get_string("search_algorithm", "simple_search");
    std::string save_stat_path = settings_blk.get_string("save_stat_path", "");
    std::string saved_result_path = settings_blk.get_string("saved_result_path", "saves/selected_final.png");
    std::string saved_textured_path = settings_blk.get_string("saved_textured_path", "saves/selected_textured.png");
    std::string saved_initial_path = settings_blk.get_string("saved_initial_path", "");
    std::string reference_path = settings_blk.get_string("reference_path", "");

    auto get_gen_params = [&](const std::vector<float> &params) -> std::vector<float>
    {
      std::vector<float> gp = std::vector<float>(params.begin(), params.begin() + gen_params_cnt);
      return gp;
    };
    auto get_camera_params = [&](const std::vector<float> &params) -> std::vector<float>
    {
      std::vector<float> gp = std::vector<float>(params.begin() + gen_params_cnt, params.end());
      return gp;
    };

    auto process_blk = [&](Block &blk){
      for (int i=0;i<blk.size();i++)
      {
        Block *pb = blk.get_block(i);
        if (!pb && pb->size()>0)
        {
          logerr("invalid parameter description\"%s\". It should be non-empty block", blk.get_name(i));
        }
        else
        {
          glm::vec2 min_max = pb->get_vec2("values", glm::vec2(1e9,-1e9));
          if (min_max.x > min_max.y)
            logerr("invalid parameter description\"%s\". It should have values:p2 with min and max values", blk.get_name(i));
          params_min.push_back(min_max.x);
          params_max.push_back(min_max.y);
          init_params.push_back(0.5*(min_max.x + min_max.y));

          int bins_cnt = pb->get_int("init_bins_count", 0);
          bool is_variant = pb->get_bool("is_variant", false);
          if (is_variant)
          {
            if (bins_cnt > 0)
            {
              logerr("invalid parameter description\"%s\". Variant parameter should not have explicit init_bins_count", blk.get_name(i));
            }
            int imin = round(min_max.x);
            int imax = round(min_max.y);
            if (abs((float)imin - min_max.x) > 1e-3 || abs((float)imax - min_max.y) > 1e-3)
            {
              logerr("invalid parameter description\"%s\". Variant parameter should have integer min max values", blk.get_name(i));
            }
            bins_cnt = imax - imin + 1;
            variant_count.push_back(bins_cnt);
            variant_positions.push_back(params_min.size()-1);
          }
          if (bins_cnt<0 || bins_cnt>512)
          {
            bins_cnt = 0;
            logerr("invalid parameter description\"%s\". Bin count should be in [0, 512] interval", blk.get_name(i));
          }

          if (bins_cnt>0)
          {
            init_bins_count.push_back((unsigned short)bins_cnt);
            init_bins_positions.push_back(params_min.size()-1);
          }
        }
      }
    };

    debug("Starting image-based optimization. Target function has %d parameters (%d for generator, %d for scene). %d SP %d var\n", 
          gen_params_cnt + scene_params_cnt, gen_params_cnt, scene_params_cnt, init_bins_count.size(), variant_count.size());

    int cameras_count = 1;
    if (by_reference)
    {
      Block *reference_cameras_blk = settings_blk.get_block("reference_cameras");
      if (reference_cameras_blk)
        cameras_count = reference_cameras_blk->size();
    }
    if (cameras_count != 1)
      logerr("ERROR: cameras_count != 1 is not supported properly.");
    process_blk(gen_params);
    for (int i=0;i<cameras_count;i++)
      process_blk(scene_params);//we duplicate scene parameters for each camera (currently scene == camera only)

    load_presets_from_blk(presets_blk, gen_params, init_params, parameter_presets);

    DiffFunctionEvaluator func;
    func.init(generator.generator, variant_positions);

    std::vector<Texture> reference_tex, reference_depth;

    CameraSettings camera;
    float h1 = 1.5;
    camera.fov_rad = 0.25;
    float h2 = h1 * tan((PI/3)/2) / tan(camera.fov_rad/2);
    camera.origin = glm::vec3(0, 0.5, h2);
    camera.target = glm::vec3(0, 0.5, 0);
    camera.up = glm::vec3(0, 1, 0);

    DepthLossCalculator dlc;

    if (by_reference)
    {
      std::vector<float> reference_params;

      settings_blk.get_arr("reference_params", reference_params);
      if (reference_params.size() != gen_params_cnt)
      {
        logerr("DOpt Error: reference_params has %d values, it should have %d", reference_params.size(), gen_params_cnt);
        return 1.0;
      }

      Block *reference_cameras_blk = settings_blk.get_block("reference_cameras");
      if (!reference_cameras_blk || reference_cameras_blk->size() == 0)
      {
        logerr("DOpt Error: reference_cameras block not found");
        return 1.0;
      }
      cameras_count = reference_cameras_blk->size();
      for (int i=0;i<reference_cameras_blk->size();i++)
      {
        std::vector<float> reference_camera_params;
        reference_cameras_blk->get_arr(i, reference_camera_params);
        std::vector<float> reference = func.get(reference_params, dgen::ModelQuality(false, 3));
        mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 512, MitsubaInterface::LLVM, 
                                                                     MitsubaInterface::TEXTURED_CONST, "texture_not_found.png"));
        mi.render_model_to_file(reference, "saves/reference.png", dgen::ModelLayout(0, 3, 6, 8, 8), camera, reference_camera_params);
        reference_tex.push_back(engine::textureManager->load_unnamed_tex("saves/reference.png"));

        Model *m = new Model();
        visualizer::simple_mesh_to_model_332(reference, m);
        m->update();
        reference_depth.push_back(dlc.get_depth(*m, camera, 128, 128));
        delete m;
      }
    }
    else
    {
      cameras_count = 1;
      reference_tex.push_back(engine::textureManager->load_unnamed_tex(reference_path));
      reference_depth.push_back(engine::textureManager->create_texture(128, 128));//TODO: estimate depth
    }

    CppAD::ADFun<float> f_reg;
    {
      std::vector<dgen::dfloat> X(init_params.size());
      for (int i = 0; i < init_params.size(); i++)
        X[i] = init_params[i];
      std::vector<dgen::dfloat> Y;
      CppAD::Independent(X);
      Y.resize(1);
      Y[0] = generator.params_regularizer(X);
      f_reg = CppAD::ADFun<float>(X, Y);
    }

    OptimizationResult opt_result{init_params, 1000, 0};

    Block *opt_settings = settings_blk.get_block(search_algorithm);
    if (!opt_settings)
    {
      logerr("Optimizer algorithm %s does not have settings block", search_algorithm.c_str());
      return 1.0;
    }
    else
    {
      int iters = 0;
      double total_time_ms = 0;

      int cnt = 0;
      std::vector<double> grad_stat;
      int model_quality = 0;

      opt::opt_func_with_grad F_silhouette = [&](std::vector<float> &params) -> std::pair<float,std::vector<float>>
      {
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        for (int i=0;i<params.size();i++)
          params[i] = CLAMP(params[i], params_min[i], params_max[i]);
        bool verbose = verbose_level > 1;
        if (verbose)
        {
          debug("params [");
          for (int j=0;j<params.size();j++)
          {
            debug("%.3f, ", params[j]);
          }
          debug("]\n");
        }
        std::vector<float> jac = func.get_jac(get_gen_params(params), dgen::ModelQuality(true, model_quality));
        std::vector<float> res = func.get(get_gen_params(params), dgen::ModelQuality(true, model_quality)); 
        std::vector<float> final_grad = std::vector<float>(params.size(), 0);

        float loss = mi.render_and_compare(res, camera, get_camera_params(params));
        mi.compute_final_grad(jac, gen_params_cnt, res.size()/FLOAT_PER_VERTEX, final_grad);
        float reg_q = 0.01;
        std::vector<float> reg_res = f_reg.Forward(0, params);
        std::vector<float> reg_jac = f_reg.Jacobian(params);
        //logerr("reg_res[0] = %f",reg_res[0]);
        loss += reg_q*MAX(0, reg_res[0]);
        for (int i=0;i<MIN(final_grad.size(), reg_jac.size());i++)
          final_grad[i] += reg_q*reg_jac[i];

        if (verbose)
        {
          debug("iter [%d] loss = %.3f\n", opt_result.total_iters, loss);
          debug("grad {");
          for (int j=0;j<final_grad.size();j++)
          {
            debug("%.3f ", final_grad[j]);
          }
          debug("}\n");
        }
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        total_time_ms += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        iters++;
        if (cnt == 0)
        {
          grad_stat = std::vector<double>(final_grad.size(), 0);
        }
        for (int j = 0; j < final_grad.size(); j++)
          grad_stat[j] += abs(final_grad[j]);
        cnt++;

        if (cnt % 100 == 0 && verbose)
        {
          debug("grad stat [");
          for (int j = 0; j < final_grad.size(); j++)
            debug("%.2f ", (float)(1000*grad_stat[j]/cnt));
          debug("]\n");
        }
        return std::pair<float,std::vector<float>>(loss, final_grad);
      };

      opt::init_params_func init_params_presets = [&parameter_presets]() -> std::vector<float>
      {
        std::vector<float> res;
        for (const auto &block : parameter_presets)
        {
          int set_n = urandi(0, block.size());
          for (float v : block[set_n])
            res.push_back(v);
        }
        return res;
      };

      constexpr int stages = 4;
      std::array<std::string, stages> optimizers = {search_algorithm, "adam", "adam", "adam"};
      std::array<int, stages> iterations = {0, 125, 100, 50};
      std::array<float, stages> lrs = {0, 0.01, 0.0075, 0.0075};
      std::array<int, stages> model_qualities = {0, 1, 1, 1};
      std::array<int, stages> image_sizes = {128, 256, 512, 1024};

      for (int stage = 0; stage < 4; stage++)
      {
        SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
        std::vector<std::string> reference_images_dir;
        std::vector<Texture> reference_mask;
        for (int i=0;i<cameras_count;i++)
        {
          reference_mask.push_back(se.get_silhouette(reference_tex[i], image_sizes[stage], image_sizes[stage]));
          reference_images_dir.push_back("saves/reference_" + std::to_string(i) + ".png");
          engine::textureManager->save_png_directly(reference_mask[i], reference_images_dir.back());
        }

        model_quality = model_qualities[stage];
        mi.init_optimization(reference_images_dir, MitsubaInterface::LOSS_MSE, 1 << 16, dgen::ModelLayout(0, 3, 3, 3, 8), 
                            MitsubaInterface::RenderSettings(image_sizes[stage], image_sizes[stage], 1, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE),
                            cameras_count, settings_blk.get_bool("save_intermediate_images", false));

        opt::Optimizer *opt = nullptr;
        if (optimizers[stage] == "adam")
          opt = new opt::Adam();
        else if (optimizers[stage] == "DE")
          opt = new opt::DifferentialEvolutionOptimizer();
        else if (optimizers[stage] == "memetic_classic")
          opt = new opt::MemeticClassic();
        else if (optimizers[stage] == "grid_search_adam")
        {
          opt_settings->add_arr("init_bins_count", init_bins_count);
          opt_settings->add_arr("init_bins_positions", init_bins_positions);
          opt = new opt::GridSearchAdam();
        }
        else
        {
          logerr("Unknown optimizer algorithm %s", optimizers[stage].c_str());
          return 1.0;
        }
        
        Block special_settings;
        special_settings.add_arr("initial_params", opt_result.best_params);
        special_settings.add_double("learning_rate", lrs[stage]);
        special_settings.add_int("iterations", iterations[stage]);
        special_settings.add_bool("verbose", true);
        Block &optimizer_settings = (stage == 0) ? *opt_settings : special_settings;

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        opt->optimize(F_silhouette, params_min, params_max, optimizer_settings, init_params_presets);
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

        double opt_time_ms = 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        opt_result.best_params = opt->get_best_result(&(opt_result.best_err));
        opt_result.total_iters = iters;
        
        debug("Stage %d Optimization stat\n", stage);
        debug("%.1f s total \n", 1e-3 * opt_time_ms);
        debug("%.1f s target function calc (%.1f ms/iter)\n", 1e-3 * total_time_ms, total_time_ms/iters);

        delete opt;
      }
    }

    std::vector<float> best_model = func.get(get_gen_params(opt_result.best_params), dgen::ModelQuality(false, 3));
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 512, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
    mi.render_model_to_file(best_model, saved_result_path, dgen::ModelLayout(), camera, get_camera_params(opt_result.best_params));
    if (saved_initial_path != "")
    {
      std::vector<float> initial_model = func.get(get_gen_params(init_params));
      mi.render_model_to_file(initial_model, saved_initial_path, dgen::ModelLayout(), camera, get_camera_params(init_params));
    }

    Texture mask_tex;
    ModelTex mt;
    {
      SilhouetteExtractor se = SilhouetteExtractor(1.0f, 0.075, 0.225);
      Texture reference_mask = se.get_silhouette(reference_tex[0], 256, 256);
      Model *m = new Model();
      std::vector<float> best_model_transformed = func.get_transformed(get_gen_params(opt_result.best_params),
                                                                       get_camera_params(opt_result.best_params),
                                                                       dgen::ModelQuality(false, 3));
      visualizer::simple_mesh_to_model_332(best_model_transformed, m);
      m->update();
      float rt_sz = MAX(reference_tex[0].get_W(), reference_tex[0].get_H());

      Texture res_tex = mt.getTexbyUV(reference_mask, *m, reference_tex[0], camera, mask_tex);
      engine::textureManager->save_png(res_tex, "reconstructed_tex");
      engine::textureManager->save_png(mask_tex, "reconstructed_mask");
    }
    
    if (textured_optimization && textured_optimization->get_bool("active", false))
    { 
      int cur_quality = 2;
      opt::opt_func_with_grad F_textured = [&](std::vector<float> &params) -> std::pair<float,std::vector<float>>
      {
        for (int i=0;i<params.size();i++)
          params[i] = CLAMP(params[i], params_min[i], params_max[i]);
        std::vector<float> jac = func.get_jac(get_gen_params(params), dgen::ModelQuality(true, 2));
        std::vector<float> res = func.get(get_gen_params(params), dgen::ModelQuality(true, 2)); 
        std::vector<float> final_grad = std::vector<float>(params.size(), 0);

        float loss = mi.render_and_compare(res, camera, get_camera_params(params));
        mi.compute_final_grad(jac, gen_params_cnt, res.size()/FLOAT_PER_VERTEX, final_grad);
        float reg_q = 0.01;
        std::vector<float> reg_res = f_reg.Forward(0, params);
        std::vector<float> reg_jac = f_reg.Jacobian(params);
        loss += reg_q*MAX(0, reg_res[0]);
        for (int i=0;i<MIN(final_grad.size(), reg_jac.size());i++)
          final_grad[i] += reg_q*reg_jac[i];

        return std::pair<float,std::vector<float>>(loss, final_grad);
      };

      opt::init_params_func init_params_null = []() -> std::vector<float>
      {
        logerr("init_params_null should never be called!!!");
        return std::vector<float>();
      };

      constexpr int stages = 3;
      std::array<int, stages> iterations = {100, 75, 50};
      std::array<float, stages> lrs = {0.0, 0.005, 0.005};
      std::array<int, stages> model_qualities = {1, 2, 2};
      std::array<int, stages> image_sizes = {128, 256, 512};
      std::array<int, stages> spps = {64, 256, 512};
      for (int stage=0;stage<stages;stage++)
      {
        cur_quality = model_qualities[stage];
        Texture reference_textured = ImageResizer::resize(reference_tex[0], image_sizes[stage], image_sizes[stage], ImageResizer::Type::CENTERED, glm::vec4(0,0,0,1));
        engine::textureManager->save_png(reference_textured, "reference_textured");
        mi.init_optimization_with_tex({"saves/reference_textured.png"}, "../../saves/reconstructed_tex.png", MitsubaInterface::LossFunction::LOSS_MSE,
                                        1 << 18, dgen::ModelLayout(0, 3, 6, 8, 8), 
                                        MitsubaInterface::RenderSettings(image_sizes[stage], image_sizes[stage], spps[stage], MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_CONST),
                                        1, true);
        Block adam_settings;
        adam_settings.add_arr("initial_params", opt_result.best_params);
        adam_settings.add_double("learning_rate", lrs[stage]);
        adam_settings.add_int("iterations", iterations[stage]);
        adam_settings.add_bool("verbose", true);
        opt::Optimizer *tex_opt = new opt::Adam();
        tex_opt->optimize(F_textured, params_min, params_max, adam_settings, init_params_null);
        opt_result.best_params = tex_opt->get_best_result(&(opt_result.best_err));

        delete tex_opt;
      }
      Texture res_optimized = engine::textureManager->load_unnamed_tex("saves/reconstructed_tex.png");
      std::vector<ModelTex::tex_data> data = {{0, 0, 1, 0.75, 3, -1}, {0, 0.75, 1, 1, 1, 1}};
      Texture comp = mt.symTexComplement(res_optimized, mask_tex, data);
      engine::textureManager->save_png(comp, "res_tex");
      sleep(1);

      std::vector<float> best_model_textured = func.get(get_gen_params(opt_result.best_params), dgen::ModelQuality(false, 3));
      mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 512, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "../../saves/reconstructed_tex.png"));
      mi.render_model_to_file(best_model_textured, saved_textured_path, dgen::ModelLayout(), camera, get_camera_params(opt_result.best_params));
    
      std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
      //debug("Optimization stat\n");
      //debug("%.1f s total \n", 1e-3 * all_time_ms);
      //debug("%.1f s target function calc (%.1f ms/iter)\n", 1e-3 * opt_time_ms, opt_time_ms/textured_optimization->get_int("iterations", 100));
    }
    debug("Model optimization finished. %d iterations total. Best result saved to \"%s\"\n", opt_result.total_iters, saved_result_path.c_str());
    debug("Best error: %f\n", opt_result.best_err);
    debug("Best params: [");
    for (int j = 0; j < opt_result.best_params.size(); j++)
    {
      debug("%.3f, ", opt_result.best_params[j]);
    }
    debug("]\n");
    mi.finish();

    return opt_result.best_err;
  }
}