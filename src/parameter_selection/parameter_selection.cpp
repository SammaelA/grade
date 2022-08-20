#include "parameter_selection.h"
#include "generation/generation_task.h"
#include "impostor_similarity.h"
#include "genetic_algorithm.h"
#include "generation/grove_packer.h"
#include "generation/grove_generator.h"
#include "tinyEngine/image.h"
#include "tree_generators/GE_generator.h"
#include "generation/metainfo_manager.h"
#include "tree_generators/all_generators.h"
#include "hydra_utils/hydra_scene_exporter.h"
#include "parameter_selection_utils.h"
#include "simulated_annealing.h"
#include "function_stat.h"
#include <thread>

bool debug_stat = false;

glm::vec4 tc_tr_mult(glm::vec4 tc_tr_1, glm::vec4 tc_tr_2)
{
    return glm::vec4(tc_tr_2.x + tc_tr_1.x * tc_tr_2.z,
                     tc_tr_2.y + tc_tr_1.y * tc_tr_2.w,
                     tc_tr_1.z * tc_tr_2.z,
                     tc_tr_1.w * tc_tr_2.w);
}

Texture load_reference(std::string name, int image_w, int image_h)
{
    Texture ref_raw = textureManager.load_unnamed_tex(image::base_img_path + name);
    Texture ref = textureManager.create_texture(image_w, image_h);
    PostFx ref_transform = PostFx("image_to_monochrome_impostor.fs");

    GLuint fbo;
    fbo = create_framebuffer();
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ref.texture, 0);
    glViewport(0, 0, image_w, image_h);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ref_transform.use();
    ref_transform.get_shader().texture("tex", ref_raw);
    ref_transform.get_shader().uniform("wood_color", glm::vec3(0.2, 0.2, 0.2));
    ref_transform.get_shader().uniform("leaves_color", glm::vec3(0, 0.15, 0));
    ref_transform.get_shader().uniform("background_color", glm::vec3(0, 0, 0));
    ref_transform.render();

    textureManager.delete_tex(ref_raw);
    delete_framebuffer(fbo);
    return ref;
}

void save_impostor_as_reference(ImpostorsData &imp, int tex_w, int tex_h, std::string name, TextureAtlas &ref_atl, glm::vec4 tc_tr)
{
    {
        TextureAtlas atl = TextureAtlas(imp.impostors.back().slices.size() * tex_w, tex_h, 1, 1);
        ref_atl = atl;
    }
    ref_atl.set_grid(tex_w, tex_h, false);
    PostFx copy = PostFx("copy_arr2.fs");
    for (int num_slice = 0; num_slice < imp.impostors.back().slices.size(); num_slice++)
    {
        glm::vec4 tex_transform = imp.atlas.tc_transform(imp.impostors.back().slices[num_slice].id);
        tex_transform = glm::vec4(tex_transform.z * tex_transform.x, tex_transform.w * tex_transform.y, tex_transform.x, tex_transform.y);
        tex_transform = tc_tr_mult(tc_tr, tex_transform);
        glm::vec3 tex = glm::vec3(0, 0, 0);
        imp.atlas.process_tc(imp.impostors.back().slices[num_slice].id, tex);
        int layer = tex.z;

        int tex_id = ref_atl.add_tex();
        ref_atl.target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex", imp.atlas.tex(0));
        copy.get_shader().uniform("tex_transform", tex_transform);
        copy.get_shader().uniform("layer", (float)layer);
        copy.render();
    }
    // textureManager.save_png(imp.atlas.tex(0), "ref_atlass_gauss_0");
    // textureManager.save_png(ref_atl.tex(0), "ref_atlass_gauss_1");
}

void vector_stat(std::vector<float> &vals, float *min_p, float *max_p, float *average_p, float *stddev_p)
{
    if (vals.empty())
        return;
    float mn = vals[0], mx = vals[0], sum = 0;
    for (float &v : vals)
    {
        mn = MIN(mn, v);
        mx = MAX(mx, v);
        sum += v;
    }
    sum /= vals.size();
    float dev = 0;
    for (float &v : vals)
    {
        dev += abs(v - sum);
    }
    dev /= vals.size();
    if (min_p)
        *min_p = mn;
    if (max_p)
        *max_p = mx;
    if (average_p)
        *average_p = sum;
    if (stddev_p)
        *stddev_p = dev;
}

float sel_quality(ParameterList &parList, ParameterList &referenceParList,
                  const std::function<std::vector<float>(std::vector<ParameterList> &)> &metric, int samples = 16)
{
    std::vector<float> res, ref;
    std::vector<ParameterList> parLists, referenceParLists;
    for (int i = 0; i < samples; i++)
    {
        parLists.push_back(parList);
        referenceParLists.push_back(referenceParList);
    }
    res = metric(parLists);
    ref = metric(referenceParLists);
    float mn1 = 0, mn2 = 0, mx1 = 0, mx2 = 0, av = 0, dev = 0;
    vector_stat(res, &mn1, &mx1, &av, &dev);
    debug("result stat [%.2f - %.2f] av=%.3f dev=%.4f\n", mn1, mx1, av, dev);
    vector_stat(ref, &mn2, &mx2, &av, &dev);
    debug("reference stat [%.2f - %.2f] av=%.3f dev=%.4f\n", mn2, mx2, av, dev);
    float q = MAX(MIN(mx1, mx2) - MAX(mn1, mn2), 1e-6) / MAX(MAX(mx1, mx2) - MIN(mn1, mn2), 1e-6);
    debug("quality: %.4f %f %f %f %f\n", q, MIN(mx1, mx2), MAX(mn1, mn2), MAX(mx1, mx2), MIN(mn1, mn2));

    return q;
}

void sel_quality(ParameterList &parList, const std::function<std::vector<float>(std::vector<ParameterList> &)> &metric,
                 int samples = 16)
{
    std::vector<float> res;
    std::vector<ParameterList> parLists;
    for (int i = 0; i < samples; i++)
    {
        parLists.push_back(parList);
    }
    res = metric(parLists);

    float mn1 = 0, mx1 = 0, av = 0, dev = 0;
    vector_stat(res, &mn1, &mx1, &av, &dev);
    debug("result stat [%.2f - %.2f] av=%.3f dev=%.4f\n", mn1, mx1, av, dev);
}

LightVoxelsCube *gen_voxels_for_selection(ReferenceTree &ref_tree)
{
    int sz_x = 25 * ceil(1.5 * ref_tree.info.BCyl_sizes.x / 25);
    int sz_y = 25 * ceil(1.5 * ref_tree.info.BCyl_sizes.y / 25);
    glm::vec3 ref_size = glm::vec3(sz_x, 0.5f * sz_y, sz_x);
    // logerr("ref size %d %d", sz_x, sz_y);
    return new LightVoxelsCube(glm::vec3(0, ref_size.y, 0), ref_size, 0.625f);
}

std::vector<float> generate_for_par_selection(std::vector<ParameterList> &params, ImpostorSimilarityCalc &imp_sim,
                                              ImpostorBaker::ImpostorGenerationParams &i_params, 
                                              std::vector<TreeTypeData> &types,
                                              ReferenceTree &ref_tree, int &cnt)
{
    textureManager.set_textures_tag(1);
    if (params.empty())
        return std::vector<float>();
    Tree *trees = new Tree[params.size()];
    GrovePacked tmp_g;
    for (int i = types.size(); i < params.size(); i++)
    {
        types.push_back(types[0]);
    }
    for (int i = 0; i < params.size(); i++)
    {
        types[i].get_params()->read_parameter_list(params[i]);
    }

    int num_threads = MIN(16, params.size());
    int step = ceil(params.size() / (float)num_threads);
    LightVoxelsCube **thr_voxels = new LightVoxelsCube *[num_threads];
    bool voxels_needed = false;
    for (auto &t : types)
    {
        auto *gen = get_generator(t.generator_name);
        if (gen->use_voxels_for_generation())
            voxels_needed = true;
        delete gen;
    }
    for (int i = 0; i < num_threads; i++)
    {
        if (voxels_needed)
            thr_voxels[i] = gen_voxels_for_selection(ref_tree);
        else
            thr_voxels[i] = new LightVoxelsCube(glm::vec3(0, 0, 0), glm::ivec3(1, 1, 1), 1.0f);
    }
    if (num_threads > 1)
    {
        std::vector<std::thread> threads;
        for (int i = 0; i < num_threads; i++)
        {
            int start_n = step * i;
            int stop_n = MIN(step * (i + 1), params.size());
            threads.push_back(std::thread(&ps_utils::gen_tree_task, start_n, stop_n, thr_voxels[i], &(types), trees));
        }

        for (auto &t : threads)
        {
            t.join();
        }
    }
    else
    {
        ps_utils::gen_tree_task(0, params.size(), thr_voxels[0], &(types), trees);
    }
    std::vector<float> res = std::vector<float>(params.size(), 1);
    int valid_cnt = 0;
    for (int i = 0; i < params.size(); i++)
    {
        if (GrovePacker::is_valid_tree(trees[i]))
        {
            valid_cnt++;
        }
        else
        {
            res[i] = 0;
        }
    }
    if (valid_cnt > 0)
    {
      GrovePacker packer;
      packer.init(Block(), types);
      packer.add_trees_to_grove(GrovePackingParams(GenerationTask::IMPOSTORS, i_params), tmp_g, trees, params.size(), false);
      cnt += params.size();
      std::vector<float> valid_res;
      imp_sim.calc_similarity(tmp_g, ref_tree, valid_res, trees, params.size(), debug_stat, true);
      int k = 0;
      for (auto &r : res)
      {
        if (r > 0)
        {
          r = valid_res[k];
          k++;
        }
      }
    }
    delete[] trees;
    for (int i = 0; i < num_threads; i++)
    {
        if (thr_voxels[i])
            delete thr_voxels[i];
    }
    delete[] thr_voxels;
    for (auto &imp : tmp_g.impostors)
    {
        imp.atlas.destroy();
    }
    textureManager.clear_unnamed_with_tag(1);
    return res;
}

void print_ref_tree_info(TreeCompareInfo &info)
{
    debug(" sz (%.3f %.3f) b_curv %.3f b_dens %.3f l_dens %.3f j_cnt %d t_th %.3f\n",
          info.BCyl_sizes.x, info.BCyl_sizes.y, info.branches_curvature,
          info.branches_density, info.leaves_density, info.joints_cnt,
          info.trunk_thickness);
}
void print_ref_tree_image_info(TreeImageInfo &info)
{
    debug(" cr_st %.3f l_sh %.3f b_sh %.3f\n", info.crown_start_level, info.crown_leaves_share, info.crown_branches_share);
}

void ParameterSelector::parameter_selection_internal(Block &selection_settings, Results &results, Scene &scene,
                                                     ReferenceTree &ref_tree, TreeTypeData *ref_type, bool save_result_image)
{
    debug("starting parameter selection for reference");
    print_ref_tree_info(ref_tree.info);
    print_ref_tree_image_info(ref_tree.image_info);
    std::string gen_name = selection_settings.get_string("generator_name", "simpliest_gen");
    float imp_size = selection_settings.get_int("impostor_size", 128);
    
    std::vector<TreeTypeData> tmp_types;
    ImpostorBaker::ImpostorGenerationParams impostor_generation_params;
    impostor_generation_params.slices_n = 1;
    impostor_generation_params.quality = imp_size;
    impostor_generation_params.monochrome = true;
    impostor_generation_params.normals_needed = false;
    impostor_generation_params.leaf_opacity = 0.75;
    impostor_generation_params.leaf_scale = 2.5;

    std::string ref_type_name = "";
    if (ref_type_name != "" && !ref_type)
    {
        tmp_types = {metainfoManager.get_tree_type(ref_type_name)};
    }
    else if (!ref_type || ref_type->generator_name != gen_name)
    {
        std::vector<TreeTypeData> types = metainfoManager.get_all_tree_types();
        bool found = false;
        for (auto &type : types)
        {
            if (type.generator_name == gen_name)
            {
                tmp_types = {type};
                found = true;
                break;
            }
        }
        if (!found)
        {
            logerr("no tree types found for generator %s. Probably it is not supported yet", gen_name.c_str());
            return;
        }
    }
    else
    {
        tmp_types = {*ref_type};
    }
    //tree_ggd.task = GenerationTask::IMPOSTORS;

    AbstractTreeGenerator::set_joints_limit(5000 * ceil(2 * ref_tree.info.joints_cnt / 5000.0f + 1));
    int imp_max_cnt = 32;
    std::string opt_name = "GA";
    my_opt::Optimizer *optimizer = nullptr;
    my_opt::MetaParameters *meta_parameters = nullptr;

    if (opt_name == "GA")
    {
        optimizer = new GeneticAlgorithm();
        GeneticAlgorithm::MetaParameters *mpp = new GeneticAlgorithm::MetaParameters();
        meta_parameters = mpp;
        GeneticAlgorithm::MetaParameters &mp = *mpp;
        mp.best_genoms_count = selection_settings.get_int("best_results_count", mp.best_genoms_count);
        mp.heaven_size = selection_settings.get_int("heaven_size", mp.heaven_size);
        mp.initial_population_size = selection_settings.get_int("initial_population_size", mp.initial_population_size);
        mp.max_population_size = selection_settings.get_int("max_population_size", mp.max_population_size);
        mp.elite = selection_settings.get_int("elite", mp.elite);
        mp.n_islands = selection_settings.get_int("n_islands", mp.n_islands);
        mp.migration_interval = selection_settings.get_int("migration_interval", mp.migration_interval);
        mp.clone_thr = selection_settings.get_double("clone_thr", mp.clone_thr);
        mp.migration_chance = selection_settings.get_double("migration_chance", mp.migration_chance);
        mp.evolution_stat = selection_settings.get_bool("evolution_stat", mp.evolution_stat);
        mp.debug_graph = selection_settings.get_bool("debug_graph", mp.debug_graph);
        mp.heaven_fine_tuning_count = selection_settings.get_int("heaven_fine_tuning_count", mp.heaven_fine_tuning_count);

        imp_max_cnt = MAX(MAX(32, mp.heaven_size * mp.heaven_recalc_n),
                          MAX(mp.initial_population_size, mp.max_population_size) +
                              (3 + 1) * (mp.heaven_fine_tuning_count) * mp.heaven_size);
    }
    else if (opt_name == "SA")
    {
        optimizer = new SimulatedAnnealing();
        SimulatedAnnealing::MetaParameters *mpp = new SimulatedAnnealing::MetaParameters();
    }

    ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(imp_max_cnt, impostor_generation_params.slices_n, false);
    ParameterList parList, bestParList;

    tmp_types[0].get_params()->write_parameter_list(parList);
    
    Block b;

    load_block_from_file(gen_name + "_param_borders.blk", b);
    parList.load_borders_from_blk(b);
    std::vector<ParameterList> initial_params;
    bool use_existing_presets = selection_settings.get_bool("use_existing_presets", true);
    if (use_existing_presets)
    {
        auto all_types = metainfoManager.get_all_tree_types();
        for (auto &t : all_types)
        {
            if (t.generator_name == tmp_types[0].generator_name)
            {
                initial_params.push_back(parList);
                t.get_params()->write_parameter_list(initial_params.back());
            }
        }
    }
    // parList.print();

    bestParList = parList;
    float best_metric = 0;
    int cnt = 0;

    std::vector<std::pair<float, ParameterList>> best_pars;
    auto func = [&](std::vector<std::vector<float>> &params) -> std::vector<float>
    {
        std::vector<ParameterList> gen_params(params.size(), parList);
        for (int i = 0; i < params.size(); i++)
        {
            gen_params[i].from_simple_list(params[i], true, true);
            // logerr("params size %d. Res size %d %d %d", params[i].size(), gen_params[i].categorialParameters.size(), gen_params[i].ordinalParameters.size(),
            //        gen_params[i].continuousParameters.size());
        }
        return generate_for_par_selection(gen_params, imp_sim, impostor_generation_params, tmp_types, ref_tree, cnt);
    };

    my_opt::ExitConditions ex_c;
    ex_c.function_calculated = selection_settings.get_int("function_calculated", ex_c.function_calculated);
    ex_c.function_reached = selection_settings.get_double("function_reached", ex_c.function_reached);
    ex_c.generations = selection_settings.get_int("generations", ex_c.generations);
    ex_c.time_elapsed_seconds = selection_settings.get_double("time_elapsed_seconds", ex_c.time_elapsed_seconds);

    my_opt::OptFunction optF;
    optF.f = func;
    optF.name = gen_name;
    optF.version = selection_settings.get_int("version", 1);

    std::vector<float> parList_f;
    parList.to_simple_list(parList_f, true, true);
    bool test_functions = false;
    bool selection = true;
    if (test_functions)
    {
        ex_c.function_calculated = 50000;
        ex_c.function_reached = 1e9;
        ((GeneticAlgorithm::MetaParameters *)(meta_parameters))->use_genes_importance = false;
        ((GeneticAlgorithm::MetaParameters *)(meta_parameters))->use_function_stat = false;
        my_opt::get_test_function("alpine100", optF, parList_f);
        int tries = 16;
        int cnt = 10;
        std::vector<glm::vec4> stat;
        for (int i = 0; i < cnt; i++)
        {
            float res = 0;
            float mn = 1;
            float mx = 0;
            float val = 25 * (i + 1);
            ((GeneticAlgorithm::MetaParameters *)(meta_parameters))->initial_population_size = val;
            for (int ii = 0; ii < tries; ii++)
            {
                optimizer = new GeneticAlgorithm();
                //((GeneticAlgorithm::MetaParameters *)(meta_parameters))->tree_GA_iters = 10 + 10*(i/2);
                std::vector<std::pair<float, std::vector<float>>> best_pars_f;
                std::vector<std::vector<float>> initial_params_f;
                for (auto &p : initial_params)
                {
                    initial_params_f.emplace_back();
                    p.to_simple_list(initial_params_f.back(), true, true);
                }
                optimizer->perform(parList_f, meta_parameters, ex_c, optF, best_pars_f, initial_params_f);
                best_pars.clear();
                for (auto &p : best_pars_f)
                {
                    best_pars.emplace_back();
                    best_pars.back().first = p.first;
                    best_pars.back().second = parList;
                    best_pars.back().second.from_simple_list(p.second, true, true);
                }

                bestParList = best_pars[0].second;
                best_metric = best_pars[0].first;
                res += best_metric;
                if (best_metric < mn)
                    mn = best_metric;
                if (best_metric > mx)
                    mx = best_metric;
                delete optimizer;
            }
            logerr("av_res = %f", res / tries);
            stat.push_back(glm::vec4(val, res / tries, mn, mx));
        }
        for (auto &v : stat)
        {
            logerr("val: %.3f, av: %.4f, [%.4f, %.4f]", v.x, v.y, v.z, v.w);
        }
        return;
    }
    else if (selection)
    {
        std::vector<std::pair<float, std::vector<float>>> best_pars_f;
        std::vector<std::vector<float>> initial_params_f;
        for (auto &p : initial_params)
        {
            initial_params_f.emplace_back();
            p.to_simple_list(initial_params_f.back(), true, true);
        }
        bool add_stat = selection_settings.get_bool("add_stat", false);
        bool replace_stat = selection_settings.get_bool("replace_stat", false);
        int stat_count = selection_settings.get_int("stat_count", 51200);
        if (add_stat || replace_stat)
        {
            my_opt::FunctionStat f_stat;
            my_opt::get_function_stat(parList_f, optF, f_stat, stat_count, replace_stat);
            f_stat.print();
        }
        optimizer->perform(parList_f, meta_parameters, ex_c, optF, best_pars_f, initial_params_f);
        for (auto &p : best_pars_f)
        {
            best_pars.emplace_back();
            best_pars.back().first = p.first;
            best_pars.back().second = parList;
            best_pars.back().second.from_simple_list(p.second, true, true);
        }

        bestParList = best_pars[0].second;
        best_metric = best_pars[0].first;
        delete optimizer;

        debug_stat = true;
        auto func_q = [&](std::vector<ParameterList> &params) -> std::vector<float>
        {
            return generate_for_par_selection(params, imp_sim, impostor_generation_params, tmp_types, ref_tree, cnt);
        };
        if (ref_type && ref_type->generator_name == gen_name)
        {
            ParameterList referenceParList;
            ref_type->get_params()->write_parameter_list(referenceParList);
            sel_quality(bestParList, referenceParList, func_q, 32);
        }
        else
        {
            sel_quality(bestParList, func_q, 32);
        }
        debug_stat = false;
    }
    else
    {
        best_pars = {std::pair<float, ParameterList>(0, bestParList)};
        AbstractTreeGenerator::set_joints_limit(1000000);
    }
    // create preapred tree
    {
        LightVoxelsCube *res_voxels = gen_voxels_for_selection(ref_tree);
        Tree *trees = new Tree[best_pars.size()];
        auto type = tmp_types[0];
        tmp_types.clear();
        for (int i = 0; i < best_pars.size(); i++)
        {
            tmp_types.emplace_back(type);
        }

        for (int i = 0; i < best_pars.size(); i++)
        {
            AbstractTreeGenerator *gen = get_generator(type.generator_name);
            glm::vec3 pos = glm::vec3(200 * (1 + i / 5), 0, 200 * (i % 5));
            res_voxels->fill(0);
            res_voxels->relocate(glm::vec3(0, res_voxels->get_center().y, 0) + pos);

            tmp_types[i].get_params()->read_parameter_list(best_pars[i].second);
            gen->plant_tree(pos, &(tmp_types[i]));
            while (gen->iterate(*res_voxels))
            {
            }
            gen->finalize_generation(trees + i, *res_voxels);
        }

      {
        GrovePacker packer;
        packer.init(Block(), tmp_types);
        auto igp = ImpostorBaker::ImpostorGenerationParams();
        if (save_result_image)
        {
            igp.slices_n = 1;
            // igp.need_top_view = false;
            igp.normals_needed = false;
            igp.quality = 512;
        }
        packer.add_trees_to_grove(GrovePackingParams(GenerationTask::IMPOSTORS | GenerationTask::MODELS, igp),
                                  scene.grove, trees, best_pars.size(), false);
      }
        int i = 0;
        for (auto &imp : scene.grove.impostors[1].impostors)
        {
            TreeCompareInfo info;
            ImpostorSimilarityCalc::get_tree_compare_info(imp, trees[i], info);
            print_ref_tree_info(info);
            i++;
        }
        delete[] trees;
        delete res_voxels;
    }

    AbstractTreeGenerator::set_joints_limit(1000000);
    results.best_candidates = tmp_types;
    results.best_res = best_metric;
}

void prepare_to_transform_reference_image(Texture &t, glm::vec3 background_color,
                                          glm::vec4 &tc_transform, float &width_height, float &tr_thick_height)
{
    unsigned char *data = nullptr;
    int w, h;
    if (!t.is_valid())
    {
        logerr("invalid reference texture");
        data = nullptr;
        w = 0;
        h = 0;
    }
    else if (t.type == GL_TEXTURE_2D)
    {
        w = t.get_W();
        h = t.get_H();
        data = new unsigned char[4 * w * h];

        glBindTexture(t.type, t.texture);

        glGetTexImage(t.type,
                      0,
                      GL_RGBA,
                      GL_UNSIGNED_BYTE,
                      data);
        glBindTexture(t.type, 0);
    }
    else
    {
        logerr("invalid reference texture format");
        data = nullptr;
        w = 0;
        h = 0;
    }
    if (data)
    {
        std::vector<glm::ivec2> min_max;
        int min_x = w;
        int max_x = -1;
        int min_y = h;
        int max_y = -1;
        for (int y = 0; y < h; y++)
        {
            min_max.push_back(glm::ivec2(-1, -1));
            for (int x = 0; x < w; x++)
            {
                int pos = 4 * ((h - y - 1) * w + x);
                glm::vec3 color = (1.0f / 255) * glm::vec3(data[pos], data[pos + 1], data[pos + 2]);
                if (length(color - background_color) > 0.05)
                {
                    min_max.back().x = x;
                    break;
                }
            }
            if (min_max.back().x >= 0)
            {
                for (int x = w - 1; x >= 0; x--)
                {
                    int pos = 4 * ((h - y - 1) * w + x);
                    glm::vec3 color = (1.0f / 255) * glm::vec3(data[pos], data[pos + 1], data[pos + 2]);
                    if (length(color - background_color) > 0.05)
                    {
                        min_max.back().y = x;
                        break;
                    }
                }
                min_x = MIN(min_x, min_max.back().x);
                max_x = MAX(max_x, min_max.back().y);
                min_y = MIN(min_y, y);
                max_y = MAX(max_y, y);
            }
            // logerr("[%d %d]", min_max.back().x, min_max.back().y);
        }
        // logerr("borders [%d %d] - [%d %d]",min_x,min_y,max_x,max_y);
        float sym_line = -1;
        for (auto &v : min_max)
        {
            if (v.x >= 0)
            {
                sym_line = (v.y + v.x) / 2.0;
                logerr("sym line %f %d %d %d %d", sym_line, v.y, v.x, max_x, min_x);
                break;
            }
        }

        tc_transform.x = (float)(min_x) / w;
        tc_transform.y = (float)(min_y) / h;
        tc_transform.z = (float)(max_x - min_x) / w;
        tc_transform.w = (float)(max_y - min_y) / h;

        float br = ReferenceTree::border_size;
        glm::vec4 tc_tr_2 = glm::vec4(-br, -br, 1, 1) / (1 - 2 * br);
        tc_transform = tc_tr_mult(tc_transform, tc_tr_2);

        width_height = 0.5 * (float)(max_x - min_x) / (max_y - min_y);

        float base_th = 0;
        float sum_th = 0;
        int len = 0;
        for (auto &v : min_max)
        {
            if (v.x >= 0)
            {
                float th = (float)(v.y - v.x) / h;
                if (base_th == 0)
                {
                    if (th < 0.25)
                    {
                        // logerr("base th [%d %d]", v.x, v.y);
                        base_th = th;
                        sum_th = th;
                        len++;
                    }
                    else
                    {
                        tr_thick_height = -1;
                        debug("cannot find visible trunk on reference image\n");
                        break;
                    }
                }
                else if (th < 1.5 * base_th)
                {
                    // logerr("th [%d %d]", v.x, v.y);
                    sum_th += th;
                    len++;
                }
                else
                {
                    break;
                }
            }
        }
        if (len > 0)
            tr_thick_height = 0.5 * sum_th / len;
        // logerr("trunk th %f", tr_thick_height);
        delete[] data;
    }
}

void get_reference_info(Block &reference_info, ReferenceTree &ref_tree, TreeCompareInfo &explicit_info, int reference_images_cnt,
                        float original_tex_aspect_ratio, float original_tex_tr_thickness)
{
    Block &ri = reference_info;
    explicit_info.BCyl_sizes = glm::vec2(0.5 * ri.get_double("width", -1), ri.get_double("height", -1));
    explicit_info.BCyl_sizes /= (1 - 2 * ReferenceTree::border_size);
    explicit_info.branches_curvature = ri.get_double("branches_curvature", -1);
    explicit_info.branches_density = ri.get_double("branches_density", -1);
    explicit_info.joints_cnt = ri.get_int("joints_cnt", -1);
    explicit_info.leaves_density = ri.get_double("leaves_density", -1);
    explicit_info.trunk_thickness = ri.get_double("trunk_thickness", -1);
    explicit_info.tropism = ri.get_double("tropism", -100);
    if (explicit_info.BCyl_sizes.x < 0)
    {
        auto status_str = ri.get_string("width", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.width_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.width_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree width is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.width_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_IMAGE" && reference_images_cnt > 0)
        {
            ref_tree.width_status = TCIFeatureStatus::FROM_IMAGE;
            if (explicit_info.BCyl_sizes.y > 0)
                ref_tree.info.BCyl_sizes.x = original_tex_aspect_ratio * explicit_info.BCyl_sizes.y;
            else
            {
                auto h_status_str = ri.get_string("height", "");
                if (h_status_str != "FROM_IMAGE")
                {
                    logerr("if width status is FROM_IMAGE height status should be either FROM_IMAGE or explicit");
                }
                else
                {
                    float max_height = ri.get_double("max_height", 200);
                    ref_tree.info.BCyl_sizes.y = max_height / 1.5;
                    ref_tree.info.BCyl_sizes.x = original_tex_aspect_ratio * ref_tree.info.BCyl_sizes.y;
                }
            }
            logerr("%f %f %f", ref_tree.info.BCyl_sizes.x, original_tex_aspect_ratio, explicit_info.BCyl_sizes.y);
        }
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.width_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for width", status_str.c_str());
    }
    else
    {
        ref_tree.info.BCyl_sizes.x = explicit_info.BCyl_sizes.x;
    }

    if (explicit_info.BCyl_sizes.y < 0)
    {
        ref_tree.height_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("height", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.height_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.height_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree height is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.height_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_IMAGE" && reference_images_cnt > 0)
        {
            ref_tree.height_status = TCIFeatureStatus::FROM_IMAGE;
        }
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.height_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for height", status_str.c_str());
    }
    else
    {
        ref_tree.info.BCyl_sizes.y = explicit_info.BCyl_sizes.y;
    }

    if (explicit_info.branches_density < 0)
    {
        ref_tree.branches_density_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("branches_density", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.branches_density_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.branches_density_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree branches_density is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.branches_density_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.branches_density_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for branches_density", status_str.c_str());
    }
    else
    {
        ref_tree.info.branches_density = explicit_info.branches_density;
    }

    if (explicit_info.leaves_density < 0)
    {
        ref_tree.leaves_density_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("leaves_density", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.leaves_density_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.leaves_density_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree leaves_density is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.leaves_density_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.leaves_density_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for leaves_density", status_str.c_str());
    }
    else
    {
        ref_tree.info.leaves_density = explicit_info.leaves_density;
    }

    if (explicit_info.tropism < -1)
    {
        ref_tree.tropism_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("tropism", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.tropism_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.tropism_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree tropism is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.tropism_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.tropism_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for tropism", status_str.c_str());
    }
    else
    {
        ref_tree.info.leaves_density = explicit_info.leaves_density;
    }

    if (explicit_info.branches_curvature < 0)
    {
        ref_tree.branches_curvature_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("branches_curvature", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.branches_curvature_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.branches_curvature_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree branches_curvature is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.branches_curvature_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.branches_curvature_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for branches_curvature", status_str.c_str());
    }
    else
    {
        ref_tree.info.branches_curvature = explicit_info.branches_curvature;
    }

    if (explicit_info.trunk_thickness < 0)
    {
        ref_tree.trunk_thickness_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("trunk_thickness", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.trunk_thickness_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.trunk_thickness_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree trunk_thickness is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
            ref_tree.trunk_thickness_status = TCIFeatureStatus::DONT_CARE;
        else if (status_str == "FROM_IMAGE" && reference_images_cnt > 0)
        {
            ref_tree.trunk_thickness_status = TCIFeatureStatus::FROM_IMAGE;
            ref_tree.info.trunk_thickness = original_tex_tr_thickness * ref_tree.info.BCyl_sizes.y;
        }
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.trunk_thickness_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for trunk_thickness", status_str.c_str());
    }
    else
    {
        ref_tree.info.trunk_thickness = explicit_info.trunk_thickness;
    }

    if (explicit_info.joints_cnt < 0)
    {
        ref_tree.joints_cnt_status = TCIFeatureStatus::DONT_CARE;
        auto status_str = ri.get_string("joints_cnt", "");
        if (status_str == "")
        {
            if (reference_images_cnt == 0)
            {
                ref_tree.joints_cnt_status = TCIFeatureStatus::FROM_TYPE;
            }
            else
            {
                ref_tree.joints_cnt_status = TCIFeatureStatus::DONT_CARE;
                debug("expected tree joints_cnt is not set. DONT_CARE mode will be used by default\n");
            }
        }
        else if (status_str == "DONT_CARE")
        {
            ref_tree.joints_cnt_status = TCIFeatureStatus::DONT_CARE;
            int max_joints = ri.get_int("max_joints", 25000);
            debug("joints_cnt status is set as DONT_CARE, but there is still a limit max_joints = %d\n", max_joints);
            ref_tree.info.joints_cnt = 0.5 * max_joints;
        }
        else if (status_str == "FROM_TYPE" && reference_images_cnt == 0)
            ref_tree.joints_cnt_status = TCIFeatureStatus::FROM_TYPE;
        else
            logerr("Parameter selection error: parameter status %s is not supported for joints_cnt", status_str.c_str());
    }
    else
    {
        ref_tree.info.joints_cnt = explicit_info.joints_cnt;
    }
}

ParameterSelector::Results ParameterSelector::parameter_selection(Block &reference_info, Block &selection_settings,
                                                                  Scene *demo_scene)
{
    Scene inner_scene;
    Scene &scene = demo_scene ? *demo_scene : inner_scene;
    scene.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(100, 100), 10);
    scene.heightmap->fill_const(0);
    float imp_size = selection_settings.get_int("impostor_size", 128);
    float original_tex_aspect_ratio = 1;
    float original_tex_tr_thickness = 0;
    ReferenceTree ref_tree;
    TreeTypeData reference_ttd;

    bool save_selection_stat = true;
    int sel_stat_w = 512;
    int sel_stat_h = 512;

    TextureAtlas sel_stat_atl = TextureAtlas(4 * sel_stat_w, sel_stat_h, 1, 1);
    sel_stat_atl.set_grid(sel_stat_w, sel_stat_h, false);
    std::string reference_name = "ERROR";
    float reference_aspect_ratio = 1;

    int reference_images_cnt = 0;
    for (int i = 0; i < reference_info.size(); i++)
    {
        if (reference_info.get_name(i) == "reference_image" && reference_info.get_type(i) == Block::ValueType::BLOCK)
            reference_images_cnt++;
    }
    if (reference_images_cnt == 0)
    {
        std::string type_name = reference_info.get_string("reference_type", "");
        if (type_name == "")
        {
            logerr("Parameter selection error: No images or plant type selected as a reference");
            return Results();
        }
        int type_id = metainfoManager.get_tree_type_id_by_name(type_name);
        if (type_id < 0)
        {
            logerr("Parameter selection error: reference type %s is invalid", type_name.c_str());
            return Results();
        }
        else
        {
            // existed type is a reference for selection
            reference_ttd = metainfoManager.get_tree_type(type_id);
            ParameterList referenceParList;
            reference_ttd.get_params()->write_parameter_list(referenceParList);
            // referenceParList.print();
            reference_ttd.get_params()->read_parameter_list(referenceParList);
            ImpostorBaker::ImpostorGenerationParams igp;
            igp.quality = imp_size;
            igp.monochrome = true;
            igp.normals_needed = false;
            igp.leaf_opacity = 0.75;
            igp.leaf_scale = 2.5;

            LightVoxelsCube *ref_voxels = new LightVoxelsCube(glm::vec3(0, 0, 0), 2.0f * reference_ttd.get_params()->get_tree_max_size(),
                                                              0.625f * reference_ttd.get_params()->get_scale_factor());
            logerr("AAAA %f %f %f", reference_ttd.get_params()->get_tree_max_size().x, reference_ttd.get_params()->get_tree_max_size().y,
                   reference_ttd.get_params()->get_tree_max_size().z);
            ReferenceTree ref_tree_init;
            ref_tree_init.info.BCyl_sizes = glm::vec2(reference_ttd.get_params()->get_tree_max_size().x, reference_ttd.get_params()->get_tree_max_size().y);
            ref_tree_init.info.joints_cnt = 100000;
            // create reference tree
            for (int i = 0; i < 10; i++)
            {
                AbstractTreeGenerator::set_joints_limit(5000 * ceil(2 * ref_tree_init.info.joints_cnt / 5000.0f + 1));
                ref_voxels = gen_voxels_for_selection(ref_tree_init);
                GrovePacked grove;
                AbstractTreeGenerator *gen = get_generator(reference_ttd.generator_name);
                ref_voxels->fill(0);
                Tree single_tree;
                gen->plant_tree(glm::vec3(0, 0, 0), &(reference_ttd));
                while (gen->iterate(*ref_voxels))
                {
                }
                gen->finalize_generation(&single_tree, *ref_voxels);
                {
                  GrovePacker packer;
                  packer.init(Block(), {reference_ttd});
                  packer.add_trees_to_grove(GrovePackingParams(GenerationTask::IMPOSTORS | GenerationTask::MODELS, igp), 
                                            grove, &single_tree, 1, false);
                }
                // save_impostor_as_reference(scene.grove.impostors[1], imp_size, imp_size, "imp_ref", ref_tree.atlas);
                // ref_atlas_transform(ref_tree.atlas);
                ImpostorSimilarityCalc::get_tree_compare_info(grove.impostors[1].impostors.back(), single_tree, ref_tree_init.info);
                delete ref_voxels;
            }
            ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(1, igp.slices_n);
            ref_voxels = gen_voxels_for_selection(ref_tree_init);
            {
                AbstractTreeGenerator::set_joints_limit(5000 * ceil(2 * ref_tree_init.info.joints_cnt / 5000.0f + 1));
                ref_voxels = gen_voxels_for_selection(ref_tree_init);
                GrovePacked grove;
                AbstractTreeGenerator *gen = get_generator(reference_ttd.generator_name);
                ref_voxels->fill(0);
                Tree single_tree;
                gen->plant_tree(glm::vec3(0, 0, 0), &(reference_ttd));
                while (gen->iterate(*ref_voxels))
                {
                }
                gen->finalize_generation(&single_tree, *ref_voxels);
                {
                  GrovePacker packer;
                  packer.init(Block(), {reference_ttd});
                  packer.add_trees_to_grove(GrovePackingParams(GenerationTask::IMPOSTORS | GenerationTask::MODELS, igp), 
                                            grove, &single_tree, 1, false);
                }
                logerr("%d", single_tree.branchHeaps[0]->branches.size());
                logerr("impostors %d", grove.impostors[1].impostors.size());
                save_impostor_as_reference(grove.impostors[1], imp_size, imp_size, "imp_ref", ref_tree.atlas, glm::vec4(0, 0, 1, 1));
                textureManager.save_png(ref_tree.atlas.tex(0), "reference_atlas_01");
                imp_sim.get_reference_tree_image_info(ref_tree, 1);
                ImpostorSimilarityCalc::get_tree_compare_info(grove.impostors[1].impostors.back(), single_tree, ref_tree.info);
                ref_tree.info.BCyl_sizes *= glm::vec2(ref_tree.image_info.tc_transform.z, ref_tree.image_info.tc_transform.w);

                float br = ReferenceTree::border_size;
                glm::vec4 tc_tr_2 = glm::vec4(-br, -br, 1, 1) / (1 - 2 * br);
                glm::vec4 tc_tr = tc_tr_mult(ref_tree.image_info.tc_transform, tc_tr_2);
                save_impostor_as_reference(grove.impostors[1], imp_size, imp_size, "imp_ref", ref_tree.atlas, tc_tr);
                textureManager.save_png(ref_tree.atlas.tex(0), "reference_atlas_02");
            }
            delete ref_voxels;
        }
    }
    else
    {
        // a set of images is a reference for selection
        {
            TextureAtlas atl = TextureAtlas(reference_images_cnt * imp_size, imp_size, 1, 1);
            ref_tree.atlas = atl;
        }
        ref_tree.atlas.set_grid(imp_size, imp_size, false);
        PostFx ref_transform = PostFx("image_to_monochrome_impostor.fs");
        original_tex_aspect_ratio = 0;
        ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(1, reference_images_cnt);

        bool ref_stat_saved = false;
        for (int i = 0; i < reference_info.size(); i++)
        {
            if (reference_info.get_name(i) == "reference_image" && reference_info.get_type(i) == Block::ValueType::BLOCK)
            {
                Block *ref_image_blk = reference_info.get_block(i);
                if (!ref_image_blk)
                    continue;
                std::string name = ref_image_blk->get_string("image", "");
                if (name == "")
                    continue;
                reference_name = name;
                Texture ref_raw = textureManager.load_unnamed_tex(image::base_img_path + name);
                reference_aspect_ratio = ref_raw.get_H() / (float)ref_raw.get_W();
                glm::vec4 ref_tc_transform;
                float aspect_ratio = -1;
                float tr_thickness = -1;
                prepare_to_transform_reference_image(ref_raw, ref_image_blk->get_vec3("background_color", glm::vec3(0, 0, 0)),
                                                     ref_tc_transform, aspect_ratio, tr_thickness);
                original_tex_aspect_ratio += aspect_ratio;
                original_tex_tr_thickness += tr_thickness;
                int slice_id = ref_tree.atlas.add_tex();
                ref_tree.atlas.target_slice(slice_id, 0);
                ref_transform.use();
                ref_transform.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
                ref_transform.get_shader().texture("tex", ref_raw);
                ref_transform.get_shader().uniform("wood_color",
                                                   ref_image_blk->get_vec3("wood_color", glm::vec3(0.2, 0.2, 0.2)));
                ref_transform.get_shader().uniform("leaves_color",
                                                   ref_image_blk->get_vec3("leaves_color", glm::vec3(0, 0.15, 0)));
                ref_transform.get_shader().uniform("background_color",
                                                   ref_image_blk->get_vec3("background_color", glm::vec3(0, 0, 0)));
                ref_transform.get_shader().uniform("ref_tc_transform", ref_tc_transform);
                ref_transform.render();
                glMemoryBarrier(GL_ALL_BARRIER_BITS);

                if (save_selection_stat && !ref_stat_saved)
                {
                    PostFx copy("copy.fs");
                    int tex_id = sel_stat_atl.add_tex();
                    sel_stat_atl.target_slice(tex_id, 0);
                    copy.use();
                    copy.get_shader().texture("tex", ref_raw);
                    copy.get_shader().uniform("tex_transform", ref_tc_transform);
                    copy.render();
                    glMemoryBarrier(GL_ALL_BARRIER_BITS);

                    PostFx copy_arr("copy_arr2.fs");
                    tex_id = sel_stat_atl.add_tex();
                    sel_stat_atl.target_slice(tex_id, 0);
                    copy_arr.use();
                    copy_arr.get_shader().texture("tex", ref_tree.atlas.tex(0));
                    copy_arr.get_shader().uniform("tex_transform", glm::vec4(0, 1, 1, -1));
                    copy_arr.render();
                    glMemoryBarrier(GL_ALL_BARRIER_BITS);
                    ref_stat_saved = true;
                }

                textureManager.delete_tex(ref_raw);
            }
        }
        original_tex_aspect_ratio /= reference_images_cnt;
        original_tex_tr_thickness /= reference_images_cnt;
        textureManager.save_png(ref_tree.atlas.tex(0), "reference_atlas_0");
        imp_sim.get_reference_tree_image_info(ref_tree);
    }
    textureManager.save_png(ref_tree.atlas.tex(0), "reference_atlas_1");
    if (save_selection_stat)
    {
        PostFx copy_arr("copy_arr2.fs");
        int tex_id = sel_stat_atl.add_tex();
        sel_stat_atl.target_slice(tex_id, 0);
        copy_arr.use();
        copy_arr.get_shader().texture("tex", ref_tree.atlas.tex(0));
        copy_arr.get_shader().uniform("tex_transform", glm::vec4(0, 1, 1, -1));
        copy_arr.render();
        glMemoryBarrier(GL_ALL_BARRIER_BITS);
    }
    TreeCompareInfo explicit_info;
    get_reference_info(reference_info, ref_tree, explicit_info, reference_images_cnt, original_tex_aspect_ratio, original_tex_tr_thickness);

    bool save_result_image = selection_settings.get_bool("save_result_image", false);
    Results res;
    parameter_selection_internal(selection_settings, res, scene, ref_tree, reference_images_cnt > 0 ? nullptr : &reference_ttd,
                                 save_result_image);

    std::string type_name = selection_settings.get_string("save_best_result", "");
    if (selection_settings.get_bool("save_best_result_auto", false))
    {
        int version = selection_settings.get_int("version", -1);
        type_name = reference_name + "_selected_params_" + std::to_string(version);
    }
    if (type_name != "")
    {
      metainfoManager.add_tree_type(res.best_candidates[0], type_name);
      metainfoManager.save_all();
    }
    
    if (save_result_image)
    {
        int tex_w = 512;
        int tex_h = 512;
        int start_imp = 1;
        auto &imp = scene.grove.impostors[1];
        int id = imp.impostors.front().slices[0].id;
        PostFx copy = PostFx("copy_arr2.fs");

        TextureAtlas atlas_copy = TextureAtlas(imp.atlas.get_sizes().x, imp.atlas.get_sizes().y, 1, 1);
        atlas_copy.set_grid(imp.atlas.get_sizes().x, imp.atlas.get_sizes().y, false);
        int tex_id = atlas_copy.add_tex();
        atlas_copy.target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex", imp.atlas.tex(0));
        copy.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        copy.get_shader().uniform("layer", (float)0);
        copy.render();
        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        // save result imposter
        std::map<int, TreeImageInfo> image_infos;
        ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(imp.impostors.size(), 1);
        imp_sim.get_tree_image_info(atlas_copy, image_infos, false);
        auto it = image_infos.find(id);
        glm::vec4 tex_transform = imp.atlas.tc_transform(id);
        tex_transform = glm::vec4(tex_transform.z * tex_transform.x, tex_transform.w * tex_transform.y, tex_transform.x, tex_transform.y);
        if (it != image_infos.end())
        {
            float br = 0.09;
            glm::vec4 tc_tr_2 = glm::vec4(-br, -br, 1, 1) / (1 - 2 * br);
            glm::vec4 tr_tc_3 = tc_tr_mult(it->second.tc_transform, tc_tr_2);
            tex_transform = tc_tr_mult(tr_tc_3, tex_transform);
        }
        else
        {
            logerr("asasff");
        }

        tex_transform = glm::vec4(tex_transform.x, tex_transform.y + tex_transform.w, tex_transform.z, -tex_transform.w);
        glm::vec3 tex = glm::vec3(0, 0, 0);
        imp.atlas.process_tc(id, tex);
        int layer = tex.z;
        tex_id = sel_stat_atl.add_tex();
        sel_stat_atl.target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex", imp.atlas.tex(0));
        copy.get_shader().uniform("tex_transform", tex_transform);
        copy.get_shader().uniform("layer", (float)layer);
        copy.render();
        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        logerr("dshdjsds %f", reference_aspect_ratio);
        TextureAtlas atlas_asp = TextureAtlas(sel_stat_atl.get_sizes().x,
                                              reference_aspect_ratio * sel_stat_atl.get_sizes().y, 1, 1);
        atlas_asp.set_grid(sel_stat_atl.get_sizes().x, reference_aspect_ratio * sel_stat_atl.get_sizes().y, false);
        tex_id = atlas_asp.add_tex();
        atlas_asp.target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex", sel_stat_atl.tex(0));
        copy.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        copy.get_shader().uniform("layer", (float)0);
        copy.render();
        glMemoryBarrier(GL_ALL_BARRIER_BITS);

        int iters = selection_settings.get_int("function_calculated", -1);
        int version = selection_settings.get_int("version", -1);
        std::string gen_name = selection_settings.get_string("generator_name", "unknown_gen");
        std::string file_name = "selection/" + gen_name + "/" + reference_name + "_v" + std::to_string(version) + "_i" + std::to_string(iters) + "_q0." + std::to_string((int)(100 * res.best_res));
        textureManager.save_png(atlas_asp.tex(0), file_name);

        // hydra scene
        {
            Block export_settings;
            float sz = 0.5 * MAX(ref_tree.info.BCyl_sizes.y, ref_tree.info.BCyl_sizes.x);
            logerr("%f %f %f", sz, ref_tree.info.BCyl_sizes.x, ref_tree.info.BCyl_sizes.y);
            export_settings.add_vec3("camera_look_at", glm::vec3(0, 53, 0));
            export_settings.add_vec3("camera_pos", glm::vec3(-18.6, 53, 8));
            export_settings.add_bool("need_terrain", true);
            export_settings.add_bool("white_terrain", true);
            export_settings.add_string("demo_copy_dir", "saves/" + file_name + "_res");
            hydra::export_scene("param_selection_scene", scene, export_settings);
        }
    }

    return res;
}