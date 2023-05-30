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
#include "graphics_utils/bilateral_filter.h"
#include "graphics_utils/resize_image.h"
#include "graphics_utils/unsharp_masking.h"
#include "graphics_utils/image_arithmetic.h"
#include "simple_model_utils.h"
#include <opencv2/opencv.hpp>
#include "compare.h"

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
      for (auto f : func_variants)
      {
        if (f.func)
          delete f.func;
      }
    }
    void init(int capacity, dgen::generator_func _model_creator, std::vector<unsigned short> &variant_positions)
    {
      model_creator = _model_creator;
      variant_params_positions = variant_positions;
      max_function_variants = capacity;
      func_variants.resize(capacity);
    }
    inline dgen::DFModel get(const std::vector<float> &params, dgen::ModelQuality mq = dgen::ModelQuality())
    {
      auto &FV = find_or_add(params, mq);
      std::vector<float> res = FV.func->Forward(0, params); 
      if (FV.model_size != res.size())
      {
        logerr("Model size(%d floats) is different from what is expected from function variant (%d floats)",
               res.size(), FV.model_size);
        logerr("This leads to incorrect work of autodiff. You should put \"is_variant:b = true\" in description of ANY generator parameter,\
                that can change the total size of model");
      }
      return {res, FV.part_offsets};
    }
    inline std::vector<float> get_jac(const std::vector<float> &params, dgen::ModelQuality mq = dgen::ModelQuality())
    {
      return find_or_add(params, mq).func->Jacobian(params); 
    }
    inline dgen::DFModel get_transformed(const std::vector<float> &params, const std::vector<float> &camera_params, 
                                       dgen::ModelQuality mq = dgen::ModelQuality())
    {
      auto &FV = find_or_add(params, mq);
      std::vector<float> f_model = FV.func->Forward(0, params); 
      dgen::transform_by_scene_parameters(camera_params, f_model);  
      return {f_model, FV.part_offsets};
    }
  private:
    struct FunctionVariant
    {
      CppAD::ADFun<float> *func = nullptr;
      dgen::ModelQuality model_quality;
      int model_size = 0;
      dgen::PartOffsets part_offsets; 
      int last_accessed = 0;
    };

    void print_FV_debug(const FunctionVariant &fv, const std::vector<unsigned short> &vs)
    {
      debug("Created new Function variant {");
      for (auto &v : vs)
        debug("%d ", (int)v);
      debug("}\n");
      debug("Expected size %d, part offsets {", fv.model_size);
      for (auto &p : fv.part_offsets)
        debug("%s : %d ", p.first.c_str(), p.second);
      debug("}\n");
    }

    FunctionVariant &find_or_add(const std::vector<float> &params, dgen::ModelQuality mq)
    {
      timestamp++;

      auto vs = get_variant_set(params);
      vs.push_back((unsigned short)mq.create_only_position);
      vs.push_back((unsigned short)mq.quality_level);
      auto it = variant_set_to_function_pos.find(vs);
      if (it == variant_set_to_function_pos.end())
      {
        //we need to create new function variant
        int func_pos = -1; 
        if (variant_set_to_function_pos.size() >= max_function_variants)
        {
          //find element with the smallest last_access and replace it with the new one
          int min_time = timestamp;
          int min_pos = 0;
          for (int i=0;i<func_variants.size();i++)
          {
            if (func_variants[i].last_accessed < min_time)
            {
              min_time = func_variants[i].last_accessed;
              min_pos = i;
            }
          }
          logerr("replacing variant %d, accessed %d (now %d)", min_pos, min_time, timestamp);
          func_pos = min_pos;
        }
        else
        {
          func_pos = next_free_function_pos;
          next_free_function_pos++;
        }
        std::vector<dgen::dfloat> X(params.size());
        for (int i=0;i<params.size();i++)
          X[i] = params[i];
        std::vector<dgen::dfloat> Y;
        CppAD::Independent(X);
        auto po = model_creator(X, Y, mq);
        CppAD::ADFun<float> *f = new CppAD::ADFun<float>(X, Y); 

        if (func_variants[func_pos].func)
          delete func_variants[func_pos].func;
        func_variants[func_pos] = FunctionVariant{f, mq, (int)Y.size(), po, timestamp};
        variant_set_to_function_pos.emplace(vs, func_pos);

        //print_FV_debug(func_variants[func_pos], vs);

        return func_variants[func_pos];
      }
      else
      {
        func_variants[it->second].last_accessed = timestamp;
        return func_variants[it->second];
      }
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
    std::vector<FunctionVariant> func_variants;
    std::vector<unsigned short> variant_params_positions;
    int timestamp = 0;
    int max_function_variants = 32;
    int next_free_function_pos = 0;
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

  void process_parameters_blk(Block &blk, std::vector<float> &params_min, std::vector<float> &params_max, std::vector<std::string> &params_names,
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
        logerr("invalid parameter description\"%s\". It should be non-empty block", blk.get_name(i).c_str());
      }
      else
      {
        glm::vec2 min_max = pb->get_vec2("values", glm::vec2(1e9, -1e9));
        if (min_max.x > min_max.y)
          logerr("invalid parameter description\"%s\". It should have values:p2 with min and max values", blk.get_name(i).c_str());
        params_min.push_back(min_max.x);
        params_max.push_back(min_max.y);
        params_names.push_back(blk.get_name(i));

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

  std::vector<float> get_params_mask(const std::vector<std::string> &all_names, const std::vector<std::string> &active_names)
  {
    std::vector<float> mask(all_names.size(), 0);

    for (const std::string &name : active_names)
    {
      for (int i=0;i<all_names.size();i++)
      {
        if (name == all_names[i])
        {
          mask[i] = 1;
          break;
        }
      }
    }

    return mask;
  }

  void inverse_mask(std::vector<float> &mask)
  {
    //expect mask to contain only 0/1 values
    for (float &v : mask)
      v = 1 - v;
  }

  std::vector<float> get_active_parameters_mask(Block &stage_blk, const std::vector<std::string> &all_names)
  {
    std::vector<std::string> pos_list, neg_list;
    stage_blk.get_arr("optimize_parameters", pos_list);
    stage_blk.get_arr("do_not_optimize_parameters", neg_list);
    if (pos_list.empty() && neg_list.empty())
    {
      return std::vector<float>(all_names.size(), 1);
    }
    else if (!pos_list.empty() && !neg_list.empty())
    {
      logerr("Error: only \"optimize_parameters\" or \"do_not_optimize_parameters\" list should be set, not both!");
      return std::vector<float>(all_names.size(), 1);
    }
    else if (!pos_list.empty())
    {
      std::vector<float> res = get_params_mask(all_names, pos_list);
      return res;
    }
    else //if (!neg_list.empty())
    {
      std::vector<float> res = get_params_mask(all_names, neg_list);
      inverse_mask(res);
      return res;
    }
  }

  MitsubaInterface::ModelInfo get_default_model_info_from_blk(Block &gen_mesh_parts)
  {
    MitsubaInterface::ModelInfo model_info;
    model_info.layout = dgen::ModelLayout(0, 3, 6, 8, 8);//default layout with pos, normals and tc

    for (int i=0;i<gen_mesh_parts.size();i++)
    {
      if (gen_mesh_parts.get_type(i) == Block::ValueType::STRING)
        model_info.parts.push_back({gen_mesh_parts.get_string(i), 
                                    "white.png", 
                                    MitsubaInterface::get_default_material()});
    }
    return model_info;
  }

  void save_opt_results(std::string save_dir, std::string prefix, DiffFunctionEvaluator &func, MitsubaInterface &mi,
                        MitsubaInterface::ModelInfo &default_model_info, dgen::GeneratorDescription &generator,
                        const CameraSettings &camera, int ref_image_size, std::vector<float> gen_params, 
                        std::vector<float> scene_params, bool save_turntable = false)
  {
    dgen::DFModel best_model = func.get(gen_params, dgen::ModelQuality(false, 3));
    MitsubaInterface::ModelInfo no_tex_model_info = default_model_info;
    if (generator.name == "buildings_2")
    {
      Block gen_info;
      load_block_from_file(generator.generator_description_blk_path, gen_info);
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
      no_tex_model_info = model_info;
    }
    mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 512, MitsubaInterface::LLVM, MitsubaInterface::TEXTURED_CONST),
                               no_tex_model_info);
    mi.render_model_to_file(best_model, save_dir+prefix+"_result.png",
                            camera, scene_params);
    
    Block params_blk;
    std::vector<float> full_params = gen_params;
    for (auto &p : scene_params)
      full_params.push_back(p);

    params_blk.add_arr("params", full_params);
    params_blk.add_arr("scene_params", scene_params);
    params_blk.add_arr("gen_params", gen_params);
    save_block_to_file("../"+save_dir+prefix+"_params.blk", params_blk);
  }

  float image_based_optimization(Block &settings_blk, MitsubaInterface &mi)
  {
    Block gen_description, scene_params, presets_blk;
    std::vector<float> params_min, params_max;
    std::vector<std::string> params_names;
    std::vector<unsigned short> variant_count, variant_positions;

    PresetsData parameter_presets;
    dgen::GeneratorDescription generator = dgen::get_generator_by_name(settings_blk.get_string("procedural_generator"));
    load_block_from_file(generator.generator_description_blk_path, gen_description);
    load_block_from_file(settings_blk.get_string("scene_description"), scene_params);
    load_block_from_file(generator.presets_blk_path, presets_blk);
    Block *_gen_params_p = gen_description.get_block("generator_parameters");
    if (!_gen_params_p)
    {
      logerr("ERROR: file %s do not have \"generator_parameters\" block", generator.generator_description_blk_path.c_str());
      return 1;
    }
    Block *_gen_mesh_parts_p = gen_description.get_block("mesh_parts");
    if (!_gen_params_p)
    {
      logerr("ERROR: file %s do not have \"mesh_parts\" block", generator.generator_description_blk_path.c_str());
      return 1;
    }
    Block &gen_params = *_gen_params_p;
    Block &gen_mesh_parts = *_gen_mesh_parts_p;
    int gen_params_cnt = gen_params.size();
    int scene_params_cnt = scene_params.size();

    int ref_image_size = settings_blk.get_int("reference_image_size", 512);
    bool by_reference = settings_blk.get_bool("synthetic_reference", true);

    auto get_gen_params = [&gen_params_cnt](const std::vector<float> &params) -> std::vector<float>
    {
      std::vector<float> gp = std::vector<float>(params.begin(), params.begin() + gen_params_cnt);
      return gp;
    };
    auto get_camera_params = [&gen_params_cnt](const std::vector<float> &params) -> std::vector<float>
    {
      std::vector<float> gp = std::vector<float>(params.begin() + gen_params_cnt, params.end());
      return gp;
    };

    process_parameters_blk(gen_params, params_min, params_max, params_names, variant_count, variant_positions);
    process_parameters_blk(scene_params, params_min, params_max, params_names, variant_count, variant_positions);
    load_presets_from_blk(presets_blk, gen_params, parameter_presets);

    bool save_results = true;
    bool save_progress = true;
    std::string save_dir = "./saves/"+settings_blk.get_string("save_folder", "opt_result")+"/";
    bool active_dir = prepare_directory(save_dir);
    save_results = active_dir;
    save_progress = active_dir;

    debug("Starting image-based optimization. Target function has %d parameters (%d for generator, %d for scene). %d variant variables \n", 
          gen_params_cnt + scene_params_cnt, gen_params_cnt, scene_params_cnt, variant_count.size());

    DiffFunctionEvaluator func;
    func.init(256, generator.generator, variant_positions);

    //there is no info about materials or textures, only the number of composite parts
    MitsubaInterface::ModelInfo default_model_info = get_default_model_info_from_blk(gen_mesh_parts);

    const int original_reference_size = settings_blk.get_int("original_reference_size", 1024);
    std::map<std::string, Texture> references;

    //if we have more than one camera, their positions must me described in either direct in
    //"cameras" block in settings blk file or indirect with "cameras_dir" string
    Block *cameras_blk = settings_blk.get_block("cameras");
    int cameras_count = cameras_blk ? cameras_blk->size() : 1;
    std::vector<CameraSettings> fixed_cameras;
    bool has_fixed_cameras = (cameras_blk != nullptr);
    if (has_fixed_cameras)
    {
      for (int i=0;i<cameras_count;i++)
      {
        Block *camera_blk = cameras_blk->get_block(i);
        if (camera_blk)
          fixed_cameras.push_back(dgen::load_camera_settings(*camera_blk));
        else
        {
          logerr("cameras_blk camera block %d corrupted!", i);
          fixed_cameras.push_back(CameraSettings());
        }
      }
    }
    else
      fixed_cameras.push_back(CameraSettings());

    for (int i=0;i<cameras_count;i++)
    {
      Block *references_blk = has_fixed_cameras ? cameras_blk->get_block(i) : settings_blk.get_block("references");
      if (!references_blk)
      {
        logerr("ERROR: camera %d references block not set", i);
        return 1;
      }
      Texture reference_tex_raw, reference_mask;
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
        dgen::DFModel reference = func.get(reference_params, dgen::ModelQuality(false, 1));

        MitsubaInterface::ModelInfo ref_mi = default_model_info;
        for (auto &p : ref_mi.parts)
          p.texture_name = "texture_not_found.png";

        mi.init_scene_and_settings(MitsubaInterface::RenderSettings(ref_image_size, ref_image_size, 512, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_CONST),
                                   ref_mi);
        CameraSettings cam = has_fixed_cameras ? fixed_cameras[i] : MitsubaInterface::get_camera_from_scene_params(reference_scene_params);
        mi.render_model_to_file(reference, "saves/reference.png", cam, reference_scene_params);
        reference_tex_raw = engine::textureManager->load_unnamed_tex("saves/reference.png");
      }
      else
      {
        reference_tex_raw = engine::textureManager->load_unnamed_tex(references_blk->get_string("textured"));
      }

      if (references_blk->get_string("mask") == "" || by_reference)
      {
        SilhouetteExtractor se = SilhouetteExtractor(MIN(reference_tex_raw.get_W(), reference_tex_raw.get_H())/256.0, 0, 0, 0.05);
        reference_mask = se.get_silhouette(reference_tex_raw, reference_tex_raw.get_W(), reference_tex_raw.get_H());
        
      }
      else
      {
        reference_mask = engine::textureManager->load_unnamed_tex(references_blk->get_string("mask"));
      }
      
      Texture reference_tex_masked = engine::textureManager->create_texture(reference_mask.get_W(),reference_mask.get_H());
      ImageArithmetics::mul(reference_tex_masked, reference_mask, reference_tex_raw, 1);
      
      references.emplace("textured_"+std::to_string(i), ImageResizer::resize(reference_tex_masked, original_reference_size, original_reference_size, ImageResizer::Type::CENTERED));
      references.emplace("mask_"+std::to_string(i), ImageResizer::resize(reference_mask, original_reference_size, original_reference_size, ImageResizer::Type::CENTERED));
      int optional_tex_cnt = 0;
      for (int i=0;i<references_blk->size();i++)
      {
        std::string tex_name = references_blk->get_name(i);
        if (references_blk->get_type(i) == Block::ValueType::STRING && tex_name != "textured" && tex_name != "mask")
        {
          Texture t = engine::textureManager->load_unnamed_tex(references_blk->get_string(i));
          t = ImageResizer::resize(t, original_reference_size, original_reference_size, ImageResizer::Type::CENTERED);
          references.emplace(tex_name+"_"+std::to_string(optional_tex_cnt), t);
          optional_tex_cnt++;
        }
      }
      if (save_results)
      {
        std::string ref_dir = save_dir + "references/";
        if (prepare_directory(ref_dir))
        {
          for (auto &p : references)
            engine::textureManager->save_png_directly(p.second, ref_dir+"reference_"+p.first+".png");
        }
      }
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

    int iters = 0;
    double total_time_ms = 0;
    int model_quality = 0;
    float regularization_alpha = settings_blk.get_double("regularization_alpha", 0.01);
    bool only_pos = true;

    opt::opt_func_with_grad F_to_optimize = [&](std::vector<float> &params) -> std::pair<float, std::vector<float>>
    {
      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

      for (int i = 0; i < params.size(); i++)
        params[i] = CLAMP(params[i], params_min[i], params_max[i]);
/*
      debug("params {");
      for (int i = 0; i < params.size(); i++)
      {
        debug("%.3f", params[i]);
        if (i != params.size()-1)
          debug(", ");
      }
      debug("}\n");
*/
      std::vector<float> jac = func.get_jac(get_gen_params(params), dgen::ModelQuality(only_pos, model_quality));
      dgen::DFModel res = func.get(get_gen_params(params), dgen::ModelQuality(only_pos, model_quality));
      std::vector<float> final_grad = std::vector<float>(params.size(), 0);

      std::vector<CameraSettings> cameras = has_fixed_cameras ? fixed_cameras : 
                                            std::vector<CameraSettings>{MitsubaInterface::get_camera_from_scene_params(get_camera_params(params))};
      float loss = mi.render_and_compare(res, cameras, get_camera_params(params));

      mi.compute_final_grad(jac, gen_params_cnt, res.first.size() / FLOAT_PER_VERTEX, final_grad);

      std::vector<float> reg_res = f_reg.Forward(0, params);
      std::vector<float> reg_jac = f_reg.Jacobian(params);
      loss += regularization_alpha * MAX(0, reg_res[0]);
      for (int i = 0; i < MIN(final_grad.size(), reg_jac.size()); i++)
        final_grad[i] += regularization_alpha * reg_jac[i];

      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      total_time_ms += 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      iters++;

      /*debug("grad [");
      for (int i = 0; i < final_grad.size(); i++)
      {
        debug("%.4f ", final_grad[i]);
      }
      debug("]\n");*/
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

    Block *silhouette_optimization_settings = settings_blk.get_block("silhouette_optimization_settings");
    if (!silhouette_optimization_settings)
    {
      logerr("Block \"silhouette_optimization_settings\" is required for optimization, but not found");
      return 1;
    }

    for (int stage = 0; stage < silhouette_optimization_settings->get_int("optimization_stages", 1); stage++)
    {
      Block *stage_blk = silhouette_optimization_settings->get_block("stage_"+std::to_string(stage));
      if (!stage_blk)
      {
        logerr("Block \"stage_%d\" is required for optimization, but not found", stage);
        return 1;
      }
      int im_sz = stage_blk->get_int("render_image_size", 128);
      std::vector<std::string> reference_image_dirs;
      for (int i=0;i<cameras_count;i++)
      {
        std::string reference_image_dir = "saves/reference_"+std::to_string(i)+".png";
        std::string reference_texture_name = stage_blk->get_string("reference_texture_name", "mask");
        reference_texture_name += "_"+std::to_string(i);
        Texture reference_mask_resized = ImageResizer::resize(references[reference_texture_name], im_sz, im_sz, ImageResizer::Type::STRETCH);
        engine::textureManager->save_png_directly(reference_mask_resized, reference_image_dir);
        reference_image_dirs.push_back(reference_image_dir);
      }

      model_quality = stage_blk->get_int("model_quality", 0);
      only_pos = false;
      mi.init_optimization(reference_image_dirs, MitsubaInterface::LOSS_MSE,
                           MitsubaInterface::RenderSettings(im_sz, im_sz, stage_blk->get_int("spp", 1), MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE),
                           default_model_info,
                           settings_blk.get_bool("save_intermediate_images", false));

      std::string search_algorithm = stage_blk->get_string("optimizer", "adam"); 
      std::vector<float> parameters_mask = get_active_parameters_mask(*stage_blk, params_names);
      opt::Optimizer *opt = get_optimizer(search_algorithm);
      Block *opt_settings = stage_blk->get_block("optimizer_settings");
      if (!opt_settings)
      {
        logerr("Optimizer algorithm %s does not have settings block", search_algorithm.c_str());
        return 1.0;
      }
      if (stage != 0)
        opt_settings->add_arr("initial_params", opt_result.best_params);
      opt_settings->add_arr("params_names", params_names);
      opt_settings->add_arr("derivatives_mult", parameters_mask);

      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
      opt->optimize(F_to_optimize, params_min, params_max, *opt_settings, init_params_presets);
      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

      double opt_time_ms = 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      opt_result.best_params = opt->get_best_result(&(opt_result.best_err));
      opt_result.total_iters = iters;

      float psnr = -opt_result.best_err;
      float mse = 1/pow(10, psnr/10);
      debug("Stage %d Model optimization stat\n", stage);
      debug("Error: PSNR = %.2f, MSE = %.5f\n", psnr, mse);
      debug("%.1f s total (%.1f ms/iter)\n", 1e-3 * opt_time_ms, opt_time_ms / iters);
      debug("Best params: [");
      for (int j = 0; j < opt_result.best_params.size(); j++)
      {
        debug("%.3f, ", opt_result.best_params[j]);
      }
      debug("]\n");

      iters = 0;
      delete opt;

      if (save_progress)
      {
        save_opt_results(save_dir, "silhouette_"+std::to_string(stage), func, mi, default_model_info, generator,
                        has_fixed_cameras ? fixed_cameras[0] : MitsubaInterface::get_camera_from_scene_params(get_camera_params(opt_result.best_params)),
                        512, 
                        get_gen_params(opt_result.best_params), 
                        get_camera_params(opt_result.best_params), 
                        false);
      }
    }

    Texture mask_tex, reconstructed_tex;
    ModelTex mt;
    {
      Model *m = new Model();
      dgen::DFModel best_model_transformed = func.get_transformed(get_gen_params(opt_result.best_params),
                                                            get_camera_params(opt_result.best_params),
                                                            dgen::ModelQuality(false, 3));
      visualizer::simple_mesh_to_model_332(best_model_transformed.first, m);
      m->update();
      for (int cam_id=0;cam_id<cameras_count;cam_id++)
      {
        //we use one (zero) camera to make initial texture and make mask as a maximum of all masks
        Texture mask_tex_cam, reconstructed_tex_cam, mask_tex_new;
        reconstructed_tex_cam = mt.getTexbyUV(references["mask_"+std::to_string(cam_id)], *m, references["textured_"+std::to_string(cam_id)],
                                          has_fixed_cameras ? fixed_cameras[cam_id] : MitsubaInterface::get_camera_from_scene_params(get_camera_params(opt_result.best_params)),
                                          mask_tex_cam);
        if (cam_id == 0)
        {
          reconstructed_tex = reconstructed_tex_cam;
          mask_tex = engine::textureManager->create_texture(reconstructed_tex.get_W(), reconstructed_tex.get_H());
        }
        mask_tex_new = engine::textureManager->create_texture(mask_tex.get_W(), mask_tex.get_H());
        ImageArithmetics::maximum(mask_tex_new, mask_tex, mask_tex_cam, 1, 1);
        mask_tex = mask_tex_new;
      }

      if (save_results)
      {
        engine::textureManager->save_png_directly(reconstructed_tex, save_dir+"reconstructed_tex.png");
        engine::textureManager->save_png_directly(mask_tex, save_dir+"reconstructed_mask.png");
      }
      engine::textureManager->save_png_directly(reconstructed_tex, "saves/reconstructed_tex.png");
    }

    Block *texture_optimization_settings = settings_blk.get_block("texture_optimization_settings");
    if (!texture_optimization_settings)
    {
      logerr("Block \"texture_optimization_settings\" is required for optimization, but not found");
      return 1;
    }
    if (texture_optimization_settings->get_int("optimization_stages", 0) > 0)
    { 
      //parameters related to light and materals cannot be optimized on silhouette stage
      //we should optimize them now along with texture
      //however, we don't want to change model shape on this stage as it may lead
      //to some "overgrowth" of model hidden with the black texture
      //That's why mask is used here
      std::vector<float> texture_only_parameters_mask = get_params_mask(params_names,
                                                                        {
                                                                          "light_translation_x",
                                                                          "light_translation_y",
                                                                          "light_translation_z",
                                                                          "light_size",
                                                                          "light_intensity",
                                                                          "ambient_light_intensity"
                                                                        });

      //choose material from the list of available ones
      std::vector<std::string> all_materials = mi.get_all_available_materials();
      std::string best_material = all_materials[0];
      std::string best_tex_name = "reconstructed_tex";

      //skip this stage if material is set manually in settings
      if (settings_blk.get_string("model_material") != "")
      {
        std::string mat_from_settings = settings_blk.get_string("model_material");
        bool valid = false;
        for (auto &material : all_materials)
        {
          if (material == mat_from_settings)
          {
            valid = true;
            break;
          }
        }
        if (valid)
          best_material = mat_from_settings;
        else
          logerr("Material \"%s\" from settings does not exist!", mat_from_settings.c_str());
        logerr("Taken material \"%s\" from settings.", best_material.c_str());
      }
      else
      {
        std::vector<float> init_params = opt_result.best_params;
        float best_val = 1e9;
      
        int mat_n = 0;
        for (auto &material : all_materials)
        {
          Block *stage_blk = texture_optimization_settings->get_block("stage_material_selection");
          if (!stage_blk)
          {
            logerr("Block \"stage_material_selection\" is required for optimization, but not found");
            return 1;
          }
          int im_sz = stage_blk->get_int("render_image_size", 128);

          model_quality = stage_blk->get_int("model_quality", 0);
          only_pos = false;
          //use only one camera to select material to make it faster
          Texture reference_textured = ImageResizer::resize(references["textured_0"], im_sz, im_sz, ImageResizer::Type::CENTERED, glm::vec4(0,0,0,1));

          dgen::DFModel cur_model = func.get(get_gen_params(opt_result.best_params), dgen::ModelQuality(false, model_quality));
          mi.init_scene_and_settings(MitsubaInterface::RenderSettings(im_sz, im_sz, 512, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE),
                                    default_model_info);
          mi.render_model_to_file(cur_model, "saves/ie_rsult3.png", 
                                  has_fixed_cameras ? fixed_cameras[0] : MitsubaInterface::get_camera_from_scene_params(get_camera_params(opt_result.best_params)),
                                  get_camera_params(opt_result.best_params));
          Texture mask_real = engine::textureManager->load_unnamed_tex("saves/ie_rsult3.png", 1);
          Texture mask_reference = ImageResizer::resize(references["mask_0"], im_sz, im_sz, ImageResizer::Type::CENTERED, glm::vec4(0,0,0,1));
          Texture mask_min = engine::textureManager->create_texture(im_sz, im_sz);
          ImageArithmetics::minimum(mask_min, mask_real, mask_reference, 1, 1);
          engine::textureManager->save_png(mask_min, "reference_textured_mask");

          Texture reference_textured_masked = engine::textureManager->create_texture(im_sz, im_sz);
          ImageArithmetics::mul(reference_textured_masked, reference_textured, mask_min, 1);
          engine::textureManager->save_png(reference_textured_masked, "reference_textured");
          std::string mat_ns = material;
          mat_ns.erase(std::remove_if(mat_ns.begin(), mat_ns.end(), isspace), mat_ns.end());
          std::string tex_name = "reconstructed_tex_" + mat_ns; 
          engine::textureManager->save_png(reconstructed_tex, tex_name);

          //TODO: select different textures and different materials for model parts
          MitsubaInterface::ModelInfo textured_model_info = default_model_info;
          textured_model_info.parts[0].material_name = material;
          textured_model_info.parts[0].texture_name = "../../saves/" + tex_name + ".png";

          mi.init_optimization_with_tex({"saves/reference_textured.png"}, MitsubaInterface::LossFunction::LOSS_MSE,
                                          MitsubaInterface::RenderSettings(im_sz, im_sz, stage_blk->get_int("spp", 128), MitsubaInterface::CUDA, 
                                          MitsubaInterface::TEXTURED_CONST),
                                          textured_model_info,
                                          stage_blk->get_double("texture_opt_lr", 0.25), settings_blk.get_bool("save_intermediate_images", false));
          std::string search_algorithm = stage_blk->get_string("optimizer", "adam"); 
          opt::Optimizer *tex_opt = get_optimizer(search_algorithm);
          Block *opt_settings = stage_blk->get_block("optimizer_settings");
          if (!opt_settings)
          {
            logerr("Optimizer algorithm %s does not have settings block", search_algorithm.c_str());
            return 1.0;
          }
          opt_settings->add_arr("initial_params", opt_result.best_params);
          opt_settings->add_arr("params_names", params_names);
          opt_settings->add_arr("derivatives_mult", texture_only_parameters_mask);

          std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
          tex_opt->optimize(F_to_optimize, params_min, params_max, *opt_settings, init_params_presets);
          std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

          double opt_time_ms = 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
          float err = 1e9;
          init_params = tex_opt->get_best_result(&err);

          if (err < best_val)
          {
            best_val = err;
            opt_result.best_err = err;
            opt_result.best_params = init_params;
            best_material = material;
            best_tex_name = tex_name;
          }

          float psnr = -err;
          float mse = 1/pow(10, psnr/10);
          debug("Stage %d/%d Material selection stat\n", mat_n + 1, all_materials.size());
          debug("Error: PSNR = %.2f, MSE = %.5f\n", psnr, mse);
          debug("%.1f s total (%.1f ms/iter)\n", 1e-3 * opt_time_ms, opt_time_ms / iters);

          iters = 0;
          total_time_ms = 0;
          mat_n++;
          delete tex_opt;
        }
        debug("Material \"%s\" selected, PSNR = %.2f\n", best_material.c_str(), -opt_result.best_err);
        sleep(3);
      }

      //when material is chosen, fine-tune the texture and scene parameters
      for (int stage = 0; stage < texture_optimization_settings->get_int("optimization_stages", 0); stage++)
      {
        Block *stage_blk = texture_optimization_settings->get_block("stage_"+std::to_string(stage));
        if (!stage_blk)
        {
          logerr("Block \"stage_%d\" is required for optimization, but not found", stage);
          return 1;
        }
        int im_sz = stage_blk->get_int("render_image_size", 128);

        model_quality = stage_blk->get_int("model_quality", 0);
        only_pos = false;
        std::vector<std::string> reference_textured_names;
        for (int i=0;i<cameras_count;i++)
        {
          std::string reference_textured_name = "reference_textured_"+std::to_string(i);
          Texture reference_textured = ImageResizer::resize(references["textured_"+std::to_string(i)], im_sz, im_sz, ImageResizer::Type::CENTERED, glm::vec4(0,0,0,1));

          dgen::DFModel cur_model = func.get(get_gen_params(opt_result.best_params), dgen::ModelQuality(false, model_quality));
          mi.init_scene_and_settings(MitsubaInterface::RenderSettings(im_sz, im_sz, 512, MitsubaInterface::LLVM, MitsubaInterface::SILHOUETTE),
                                    default_model_info);
          mi.render_model_to_file(cur_model, "saves/ie_rsult3.png", 
                                  has_fixed_cameras ? fixed_cameras[i] : MitsubaInterface::get_camera_from_scene_params(get_camera_params(opt_result.best_params)), 
                                  get_camera_params(opt_result.best_params));
          Texture mask_real = engine::textureManager->load_unnamed_tex("saves/ie_rsult3.png", 1);
          Texture mask_reference = ImageResizer::resize(references["mask_"+std::to_string(i)], im_sz, im_sz, ImageResizer::Type::CENTERED, glm::vec4(0,0,0,1));
          Texture mask_min = engine::textureManager->create_texture(im_sz, im_sz);
          ImageArithmetics::minimum(mask_min, mask_real, mask_reference, 1, 1);
          engine::textureManager->save_png(mask_min, reference_textured_name+"_mask");

          Texture reference_textured_masked = engine::textureManager->create_texture(im_sz, im_sz);
          ImageArithmetics::mul(reference_textured_masked, reference_textured, mask_min, 1);
          engine::textureManager->save_png(reference_textured_masked, reference_textured_name);
          reference_textured_names.push_back("saves/"+reference_textured_name+".png");
        }


        //TODO: select different textures and different materials for model parts
        MitsubaInterface::ModelInfo textured_model_info = default_model_info;
        textured_model_info.parts[0].material_name = best_material;
        textured_model_info.parts[0].texture_name = "../../saves/" + best_tex_name + ".png";
        mi.init_optimization_with_tex(reference_textured_names, MitsubaInterface::LossFunction::LOSS_MSE, 
                                      MitsubaInterface::RenderSettings(im_sz, im_sz, stage_blk->get_int("spp", 128), MitsubaInterface::CUDA, 
                                        MitsubaInterface::TEXTURED_CONST),
                                      textured_model_info,
                                      stage_blk->get_double("texture_opt_lr", 0.25), true);
        
        std::string search_algorithm = stage_blk->get_string("optimizer", "adam"); 
        opt::Optimizer *tex_opt = get_optimizer(search_algorithm);
        Block *opt_settings = stage_blk->get_block("optimizer_settings");
        if (!opt_settings)
        {
          logerr("Optimizer algorithm %s does not have settings block", search_algorithm.c_str());
          return 1.0;
        }
        opt_settings->add_arr("initial_params", opt_result.best_params);
        opt_settings->add_arr("params_names", params_names);
        opt_settings->add_arr("derivatives_mult", texture_only_parameters_mask);

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        tex_opt->optimize(F_to_optimize, params_min, params_max, *opt_settings, init_params_presets);
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

        double opt_time_ms = 1e-3 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        opt_result.best_params = tex_opt->get_best_result(&(opt_result.best_err));
        opt_result.total_iters = iters;

        float psnr = -opt_result.best_err;
        float mse = 1/pow(10, psnr/10);
        debug("Stage %d Texture optimization stat\n", stage);
        sleep(5);
        debug("Error: PSNR = %.2f, MSE = %.5f\n", psnr, mse);
        debug("%.1f s total (%.1f ms/iter)\n", 1e-3 * opt_time_ms, opt_time_ms / iters);
        debug("Best params: [");
        for (int j = 0; j < opt_result.best_params.size(); j++)
        {
          debug("%.3f, ", opt_result.best_params[j]);
        }
        debug("]\n");

        iters = 0;
        total_time_ms = 0;
        delete tex_opt;
      }

      Texture res_optimized = engine::textureManager->load_unnamed_tex("saves/" + best_tex_name + ".png");

      Texture res = BilateralFilter::perform(res_optimized, 2.5, 0.33);
      Texture comp, mask_complemented;
      if (cameras_count == 1)
      {
        std::vector<ModelTex::tex_data> data = {{0, 0, 1, 0.75, 3, 1}, {0, 0.75, 1, 1, 1, 4}};
        comp = mt.symTexComplement(res_optimized, mask_tex, data, &mask_complemented);
      }
      else
      {
        comp = res_optimized;
        mask_complemented = mask_tex;
      }

      engine::textureManager->save_png(res_optimized, "reconstructed_tex_raw");
      engine::textureManager->save_png(comp, "reconstructed_tex_complemented");
      engine::textureManager->save_png(mask_complemented, "reconstructed_mask_complemented");

      if (save_results)
      {
        engine::textureManager->save_png_directly(res_optimized, save_dir+"reconstructed_tex_raw.png");
        engine::textureManager->save_png_directly(comp, save_dir+"reconstructed_tex_complemented.png");
        engine::textureManager->save_png_directly(mask_complemented, save_dir+"reconstructed_mask_complemented.png");       

        cv::Mat image, mask, image_inpainted;
        image = cv::imread(save_dir+"reconstructed_tex_complemented.png");
        mask = cv::imread(save_dir+"reconstructed_mask_complemented.png", cv::ImreadModes::IMREAD_GRAYSCALE);
        for (int i=0;i<mask.size().height;i++)
          for (int j=0;j<mask.size().width;j++)
            mask.at<unsigned char>(i, j) = mask.at<unsigned char>(i, j) > 128 ? 0 : 1;

        cv::inpaint(image, mask, image_inpainted, 16, cv::INPAINT_TELEA);
        cv::imwrite(save_dir+"reconstructed_tex_complemented.png", image_inpainted);

        sleep(1);

        auto model_quality = dgen::ModelQuality(false, 3);
        if (generator.name == "buildings_2")
          model_quality = dgen::ModelQuality(false, 0);

        dgen::DFModel best_model_textured = func.get(get_gen_params(opt_result.best_params), model_quality);
        MitsubaInterface::ModelInfo textured_model_info = default_model_info;
        textured_model_info.parts[0].material_name = best_material;
        CameraSettings demo_cam = has_fixed_cameras ? fixed_cameras[0] : MitsubaInterface::get_camera_from_scene_params(get_camera_params(opt_result.best_params));
        
        textured_model_info.parts[0].texture_name = "../../"+save_dir+"reconstructed_tex_raw.png";
        mi.init_scene_and_settings(MitsubaInterface::RenderSettings(1024, 1024, 512, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_CONST),
                                  textured_model_info);
        mi.render_model_to_file(best_model_textured, save_dir+"selected_textured_raw.png", demo_cam, get_camera_params(opt_result.best_params));

        textured_model_info.parts[0].texture_name = "../../"+save_dir+"reconstructed_tex_complemented.png";
        mi.init_scene_and_settings(MitsubaInterface::RenderSettings(1024, 1024, 512, MitsubaInterface::CUDA, MitsubaInterface::TEXTURED_CONST),
                                  textured_model_info);
        mi.render_model_to_file(best_model_textured, save_dir+"selected_textured_complemented.png", demo_cam, get_camera_params(opt_result.best_params));

        if (save_progress)
        {
          float fov_rad = 0.5;
          CameraSettings turntable_cam = MitsubaInterface::get_camera_from_scene_params({fov_rad});
          render_mygen_cup(mi, turntable_cam, get_gen_params(opt_result.best_params), save_dir+"reconstructed_tex_complemented.png", save_dir);
        }
      }
    }
    debug("Optimization finished. %d iterations total. Best result saved to \"%s\"\n", opt_result.total_iters, save_dir.c_str());
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