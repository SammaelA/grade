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
#include "graphics_utils/image_expand.h"
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

  void load_presets_from_blk(Block &presets_blk, Block &gen_params, PresetsData &presets_data /*output*/)
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

  void process_parameters_blk(Block &blk, std::vector<float> &params_min, std::vector<float> &params_max,
                              std::vector<unsigned short> &variant_count, std::vector<unsigned short> &variant_positions)
  {
    /*
    parameters blk is a list of blocks with the following structure
    parameter_name
    {
      values:p2 = min_value, max_value
      [is_variant:b = true/false]
      [desc:s = "description"]
    }
    */
    for (int i = 0; i < blk.size(); i++)
    {
      Block *pb = blk.get_block(i);
      if (!pb && pb->size() > 0)
      {
        logerr("invalid parameter description\"%s\". It should be non-empty block", blk.get_name(i));
      }
      else
      {
        glm::vec2 min_max = pb->get_vec2("values", glm::vec2(1e9, -1e9));
        if (min_max.x > min_max.y)
          logerr("invalid parameter description\"%s\". It should have values:p2 with min and max values", blk.get_name(i));
        params_min.push_back(min_max.x);
        params_max.push_back(min_max.y);

        bool is_variant = pb->get_bool("is_variant", false);
        if (is_variant)
        {
          int imin = round(min_max.x);
          int imax = round(min_max.y);
          if (abs((float)imin - min_max.x) > 1e-3 || abs((float)imax - min_max.y) > 1e-3)
          {
            logerr("invalid parameter description\"%s\". Variant parameter should have integer min max values", blk.get_name(i));
          }
          variant_count.push_back(imax - imin + 1);
          variant_positions.push_back(params_min.size() - 1);
        }
      }
    }
  }

  opt::Optimizer *get_optimizer(const std::string &opt_name)
  {
    if (opt_name == "adam")
      return new opt::Adam();
    else if (opt_name == "DE")
      return new opt::DifferentialEvolutionOptimizer();
    else if (opt_name == "memetic_classic")
      return new opt::MemeticClassic();
    else if (opt_name == "grid_search_adam")
      return new opt::GridSearchAdam();
    else
    {
      logerr("Unknown optimizer algorithm %s", opt_name.c_str());
      return new opt::Adam();
    }
  }

  float image_based_optimization(Block &settings_blk, MitsubaInterface &mi)
  {
    Block gen_params, scene_params, presets_blk;
    std::vector<float> params_min, params_max;
    std::vector<unsigned short> variant_count, variant_positions;

    PresetsData parameter_presets;
    dgen::GeneratorDescription generator = dgen::get_generator_by_name(settings_blk.get_string("procedural_generator"));
    load_block_from_file(generator.parameters_description_blk_path, gen_params);
    load_block_from_file(settings_blk.get_string("scene_description"), scene_params);
    load_block_from_file(generator.presets_blk_path, presets_blk);
    int gen_params_cnt = gen_params.size();
    int scene_params_cnt = scene_params.size();

    int ref_image_size = settings_blk.get_int("reference_image_size", 512);
    int sel_image_size = settings_blk.get_int("selection_image_size", 196);
    bool by_reference = settings_blk.get_bool("synthetic_reference", true);
    Block *textured_optimization = settings_blk.get_block("textured_optimization");
    std::string search_algorithm = settings_blk.get_string("search_algorithm", "simple_search");
    std::string save_stat_path = settings_blk.get_string("save_stat_path", "");
    std::string saved_result_path = settings_blk.get_string("saved_result_path", "saves/selected_final.png");
    std::string saved_textured_path = settings_blk.get_string("saved_textured_path", "saves/selected_textured.png");
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

    process_parameters_blk(gen_params, params_min, params_max, variant_count, variant_positions);
    process_parameters_blk(scene_params, params_min, params_max, variant_count, variant_positions);
    load_presets_from_blk(presets_blk, gen_params, parameter_presets);

    debug("Starting image-based optimization. Target function has %d parameters (%d for generator, %d for scene). %d variant variables \n", 
          gen_params_cnt + scene_params_cnt, gen_params_cnt, scene_params_cnt, variant_count.size());

    DiffFunctionEvaluator func;
    func.init(generator.generator, variant_positions);

    CameraSettings camera = CameraSettings::get_default_mitsuba_preset();

    Texture reference_tex;
    Texture reference_mask;
    const int original_reference_size = 1024;

    {
      Texture reference_tex_raw;
      if (by_reference)
      {
        //create synthetic reference by rendering model with given parameters with mitsuba
        std::vector<float> reference_params;

        settings_blk.get_arr("reference_params", reference_params);
        if (reference_params.size() != gen_params_cnt)
        {
          logerr("DOpt Error: reference_params has %d values, it should have %d", reference_params.size(), gen_params_cnt);
          return 1.0;
        }
        std::vector<float> reference_scene_params;
        settings_blk.get_arr("reference_scene", reference_scene_params);
        std::vector<float> reference = func.get(reference_params, dgen::ModelQuality(false, 3));
        mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 512, MitsubaInterface::LLVM,
                                                                    MitsubaInterface::TEXTURED_CONST, "texture_not_found.png"));
        mi.render_model_to_file(reference, "saves/reference.png", dgen::ModelLayout(0, 3, 6, 8, 8), camera, reference_scene_params);
        reference_tex_raw = engine::textureManager->load_unnamed_tex("saves/reference.png");
      }
      else
      {
        reference_tex_raw = engine::textureManager->load_unnamed_tex(reference_path);
      }

      reference_tex = ImgExp::ImgExpanding(reference_tex_raw, original_reference_size);
      SilhouetteExtractor se = SilhouetteExtractor(0, 0.075, 0.225, 0.01);
      reference_mask = se.get_silhouette_simple(reference_tex, original_reference_size, original_reference_size);
      engine::textureManager->save_png(reference_mask, "ie_rsult2.png");
    }

    CppAD::ADFun<float> f_reg;
    {
      std::vector<dgen::dfloat> X(params_max.size());
      for (int i = 0; i < X.size(); i++)
        X[i] = i; //some random numbers, doesn't matter
      std::vector<dgen::dfloat> Y;
      CppAD::Independent(X);
      Y.resize(1);
      Y[0] = generator.params_regularizer(X);
      f_reg = CppAD::ADFun<float>(X, Y);
    }

    OptimizationResult opt_result{params_max, 1000, 0};

    Block *opt_settings = settings_blk.get_block(search_algorithm);
    if (!opt_settings)
    {
      logerr("Optimizer algorithm %s does not have settings block", search_algorithm.c_str());
      return 1.0;
    }

    int iters = 0;
    double total_time_ms = 0;
    int model_quality = 0;

    opt::opt_func_with_grad F_silhouette = [&](std::vector<float> &params) -> std::pair<float, std::vector<float>>
    {
      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

      for (int i = 0; i < params.size(); i++)
        params[i] = CLAMP(params[i], params_min[i], params_max[i]);

      std::vector<float> jac = func.get_jac(get_gen_params(params), dgen::ModelQuality(true, model_quality));
      std::vector<float> res = func.get(get_gen_params(params), dgen::ModelQuality(true, model_quality));
      std::vector<float> final_grad = std::vector<float>(params.size(), 0);

      float loss = mi.render_and_compare(res, camera, get_camera_params(params));

      mi.compute_final_grad(jac, gen_params_cnt, res.size() / FLOAT_PER_VERTEX, final_grad);

      float reg_q = 0.01;
      std::vector<float> reg_res = f_reg.Forward(0, params);
      std::vector<float> reg_jac = f_reg.Jacobian(params);
      loss += reg_q * MAX(0, reg_res[0]);
      for (int i = 0; i < MIN(final_grad.size(), reg_jac.size()); i++)
        final_grad[i] += reg_q * reg_jac[i];

      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      total_time_ms += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      iters++;

      return std::pair<float, std::vector<float>>(loss, final_grad);
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

    for (int stage = 0; stage < 2; stage++)
    {
      std::string reference_image_dir = "saves/reference.png";
      Texture reference_mask_resized = ImageResizer::resize(reference_mask, image_sizes[stage], image_sizes[stage], ImageResizer::Type::STRETCH);
      engine::textureManager->save_png_directly(reference_mask_resized, reference_image_dir);

      model_quality = model_qualities[stage];
      mi.init_optimization({reference_image_dir}, MitsubaInterface::LOSS_MSE, 1 << 16, dgen::ModelLayout(0, 3, 3, 3, 8),
                           MitsubaInterface::RenderSettings(image_sizes[stage], image_sizes[stage], 1, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE),
                           1, settings_blk.get_bool("save_intermediate_images", false));

      opt::Optimizer *opt = get_optimizer(optimizers[stage]);

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

      debug("Stage %d Model optimization stat\n", stage);
      debug("%.1f s total (%.1f ms/iter)\n", 1e-3 * opt_time_ms, opt_time_ms / iters);

      delete opt;
    }

    std::vector<float> best_model = func.get(get_gen_params(opt_result.best_params), dgen::ModelQuality(false, 3));
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 512, MitsubaInterface::LLVM, MitsubaInterface::MONOCHROME));
    mi.render_model_to_file(best_model, saved_result_path, dgen::ModelLayout(0, 3, 6, 8, 8), camera, get_camera_params(opt_result.best_params));

    Texture mask_tex;
    ModelTex mt;
    {
      Model *m = new Model();
      std::vector<float> best_model_transformed = func.get_transformed(get_gen_params(opt_result.best_params),
                                                                       get_camera_params(opt_result.best_params),
                                                                       dgen::ModelQuality(false, 3));
      visualizer::simple_mesh_to_model_332(best_model_transformed, m);
      m->update();
      float rt_sz = MAX(reference_tex.get_W(), reference_tex.get_H());

      Texture res_tex = mt.getTexbyUV(reference_mask, *m, reference_tex, camera, mask_tex);
      engine::textureManager->save_png(res_tex, "reconstructed_tex");
      engine::textureManager->save_png(mask_tex, "reconstructed_mask");
    }
    
    if (textured_optimization && textured_optimization->get_bool("active", false))
    { 
      iters = 0;
      total_time_ms = 0;

      constexpr int stages = 3;
      std::array<int, stages> iterations = {50, 50, 50};
      std::array<float, stages> lrs = {0.0, 0.0, 0.0}; //learning rate for model shape optimization. Temporary disabled
      std::array<int, stages> model_qualities = {1, 2, 2};
      std::array<int, stages> image_sizes = {128, 256, 512};
      std::array<int, stages> spps = {64, 256, 512};
      for (int stage = 0; stage < 2; stage++)
      {
        model_quality = model_qualities[stage];
        Texture reference_textured = ImageResizer::resize(reference_tex, image_sizes[stage], image_sizes[stage], ImageResizer::Type::CENTERED, glm::vec4(0,0,0,1));
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

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        tex_opt->optimize(F_silhouette, params_min, params_max, adam_settings, init_params_presets);
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

        double opt_time_ms = 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        opt_result.best_params = tex_opt->get_best_result(&(opt_result.best_err));
        opt_result.total_iters = iters;

        debug("Stage %d Texture optimization stat\n", stage);
        debug("%.1f s total (%.1f ms/iter)\n", 1e-3 * opt_time_ms, opt_time_ms / iters);

        delete tex_opt;
      }

      Texture res_optimized = engine::textureManager->load_unnamed_tex("saves/reconstructed_tex.png");
      std::vector<ModelTex::tex_data> data = {{0, 0, 1, 0.75, 3, -1}, {0, 0.75, 1, 1, 1, 1}};
      Texture comp = mt.symTexComplement(res_optimized, mask_tex, data);
      engine::textureManager->save_png(comp, "res_tex");
      
      sleep(1);

      std::vector<float> best_model_textured = func.get(get_gen_params(opt_result.best_params), dgen::ModelQuality(false, 3));
      mi.init_scene_and_settings(MitsubaInterface::RenderSettings(512, 512, 512, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST, "../../saves/reconstructed_tex.png"));
      mi.render_model_to_file(best_model_textured, saved_textured_path, dgen::ModelLayout(0, 3, 6, 8, 8), camera, get_camera_params(opt_result.best_params));
    }

    debug("Optimization finished. %d iterations total. Best result saved to \"%s\"\n", opt_result.total_iters, saved_result_path.c_str());
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