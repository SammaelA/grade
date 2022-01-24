#include "parameter_selection.h"
#include "generation/generation_settings.h"
#include "impostor_similarity.h"
#include "genetic_algorithm.h"
#include "generation/grove_packer.h"
#include "generation/grove_generator.h"
#include "tinyEngine/image.h"
#include "tree_generators/GE_generator.h"
#include "generation/metainfo_manager.h"
#include <thread>

Texture load_reference(std::string name, int image_w, int image_h)
{
    Texture ref_raw = textureManager.load_unnamed_tex(image::base_img_path + name);
    Texture ref = textureManager.create_unnamed(image_w, image_h);
    PostFx ref_transform = PostFx("image_to_monochrome_impostor.fs");

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
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
    glDeleteFramebuffers(1, &fbo);
    return ref;
}

void save_impostor_as_reference(ImpostorsData &imp, int tex_w, int tex_h, std::string name, TextureAtlas &ref_atl)
{
    {
        TextureAtlas atl = TextureAtlas(8 * tex_w, tex_h, 1, 1);
        ref_atl = atl;
    }
    ref_atl.set_grid(tex_w, tex_h, false);
    PostFx copy = PostFx("copy_arr2.fs");
    for (int num_slice = 0; num_slice < 8; num_slice++)
    {
        glm::vec4 tex_transform = imp.atlas.tc_transform(imp.impostors.back().slices[num_slice].id);
        tex_transform = glm::vec4(tex_transform.z * tex_transform.x, tex_transform.w * tex_transform.y, tex_transform.x, tex_transform.y);
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
    //textureManager.save_png(imp.atlas.tex(0), "ref_atlass_gauss_0");
    //textureManager.save_png(ref_atl.tex(0), "ref_atlass_gauss_1");
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

void ref_atlas_transform(TextureAtlas &atl)
{
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glm::ivec4 sizes = atl.get_sizes();
    glm::ivec2 slice_size = atl.get_slice_size();
    TextureAtlas atl_tmp = TextureAtlas(sizes.x, sizes.y, atl.layers_count(), 1);
    atl_tmp.set_grid(slice_size.x, slice_size.y, false);

    PostFx gauss = PostFx("gaussian_blur_atlas.fs");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl_tmp.target(l, 0);
        gauss.use();
        gauss.get_shader().texture("tex", atl.tex(0));
        gauss.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        gauss.get_shader().uniform("layer", (float)l);
        gauss.get_shader().uniform("pass", 0);
        gauss.get_shader().uniform("tex_size_inv", glm::vec2(1.0f / sizes.x, 1.0f / sizes.y));
        gauss.get_shader().uniform("slice_size", slice_size);
        gauss.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl_tmp.tex(0), "atlass_gauss_0");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl.target(l, 0);
        gauss.use();
        gauss.get_shader().texture("tex", atl_tmp.tex(0));
        gauss.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        gauss.get_shader().uniform("layer", l);
        gauss.get_shader().uniform("pass", 1);
        gauss.get_shader().uniform("tex_size_inv", glm::vec2(1.0f / sizes.x, 1.0f / sizes.y));
        gauss.get_shader().uniform("slice_size", slice_size);
        gauss.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl.tex(0), "atlass_gauss_1");

    PostFx sil_fill = PostFx("silhouette_fill.fs");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl_tmp.target(l, 0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl.tex(0));
        sil_fill.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        sil_fill.get_shader().uniform("layer", (float)l);
        sil_fill.get_shader().uniform("radius", 8);
        sil_fill.get_shader().uniform("dir_threshold", 4);
        sil_fill.get_shader().uniform("tex_size_inv", glm::vec2(1.0f / sizes.x, 1.0f / sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", slice_size);
        sil_fill.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl.target(l, 0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl_tmp.tex(0));
        sil_fill.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        sil_fill.get_shader().uniform("layer", (float)l);
        sil_fill.get_shader().uniform("radius", 4);
        sil_fill.get_shader().uniform("dir_threshold", 6);
        sil_fill.get_shader().uniform("tex_size_inv", glm::vec2(1.0f / sizes.x, 1.0f / sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", slice_size);
        sil_fill.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl_tmp.target(l, 0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl.tex(0));
        sil_fill.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        sil_fill.get_shader().uniform("layer", (float)l);
        sil_fill.get_shader().uniform("radius", 4);
        sil_fill.get_shader().uniform("dir_threshold", 6);
        sil_fill.get_shader().uniform("tex_size_inv", glm::vec2(1.0f / sizes.x, 1.0f / sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", slice_size);
        sil_fill.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl_tmp.tex(0), "atlass_gauss_2");

    PostFx sil_sharp = PostFx("silhouette_sharpen.fs");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl.target(l, 0);
        sil_sharp.use();
        sil_sharp.get_shader().texture("tex", atl_tmp.tex(0));
        sil_sharp.get_shader().uniform("tex_transform", glm::vec4(0, 0, 1, 1));
        sil_sharp.get_shader().uniform("layer", (float)l);
        sil_sharp.get_shader().uniform("threshold", 0.05f);
        sil_sharp.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl.tex(0), "atlass_gauss_3");
}

void gen_tree(LightVoxelsCube &voxels, TreeTypeData *type, Tree *tree)
{
    AbstractTreeGenerator *gen = GroveGenerator::get_generator(type->generator_name);
    voxels.fill(0);
    gen->plant_tree(glm::vec3(0, 0, 0), type);
    while (gen->iterate(voxels))
    {
    }
    gen->finalize_generation(tree, voxels);
}
void gen_tree_task(int start_n, int stop_n, LightVoxelsCube *vox, std::vector<TreeTypeData> *types, Tree *trees)
{
    for (int i = start_n; i < stop_n; i++)
    {
        //logerr("gen tree %d in %d %d", i, start_n, stop_n);
        gen_tree(*vox, &((*types)[i]), trees + i);
    }
}

LightVoxelsCube *gen_voxels_for_selection(ReferenceTree &ref_tree)
{
    glm::vec3 ref_size = glm::vec3(ref_tree.info.BCyl_sizes.x, ref_tree.info.BCyl_sizes.y, ref_tree.info.BCyl_sizes.x);
    return new LightVoxelsCube(glm::vec3(0, 0, 0), 2.0f * ref_size, 0.625f);
}

std::vector<float> generate_for_par_selection(std::vector<ParameterList> &params, ImpostorSimilarityCalc &imp_sim,
                                              GroveGenerationData &tree_ggd, Heightmap *flat_hmap,
                                              ReferenceTree &ref_tree, int &cnt)
{
    if (params.empty())
        return std::vector<float>();
    GrovePacker packer;
    Tree *trees = new Tree[params.size()];
    GrovePacked tmp_g;
    for (int i = tree_ggd.types.size(); i < params.size(); i++)
    {
        tree_ggd.types.push_back(tree_ggd.types[0]);
    }
    tree_ggd.trees_count = params.size();
    for (int i = 0; i < params.size(); i++)
    {
        tree_ggd.types[i].params->read_parameter_list(params[i]);
    }

    int num_threads = MIN(8, params.size());
    int step = ceil(params.size() / (float)num_threads);
    LightVoxelsCube **thr_voxels = new LightVoxelsCube *[num_threads];
    for (int i = 0; i < num_threads; i++)
    {
        thr_voxels[i] = gen_voxels_for_selection(ref_tree);
    }
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++)
    {
        int start_n = step * i;
        int stop_n = MIN(step * (i + 1), params.size());
        threads.push_back(std::thread(&gen_tree_task, start_n, stop_n, thr_voxels[i], &(tree_ggd.types), trees));
    }

    for (auto &t : threads)
    {
        t.join();
    }
    packer.add_trees_to_grove(tree_ggd, tmp_g, trees, flat_hmap, false);
    ref_atlas_transform(tmp_g.impostors[1].atlas);
    //textureManager.save_png(tmp_g.impostors[1].atlas.tex(0),"imp"+std::to_string(cnt));
    //logerr("generate %d",cnt);
    cnt += params.size();
    std::vector<float> res;
    imp_sim.calc_similarity(tmp_g, ref_tree, res, trees);
    delete[] trees;
    for (int i = 0; i < num_threads; i++)
    {
        delete thr_voxels[i];
    }
    delete[] thr_voxels;
    return res;
}
void ParameterSelector::parameter_selection_internal(Block &selection_settings, Results &results, Scene &scene,
                                                     ReferenceTree &ref_tree, TreeTypeData *ref_type)
{
    debug("starting parameter selection for reference (%.3f %.3f) %.3f %.3f %d %.3f %.3f\n", 
          ref_tree.info.BCyl_sizes.x, ref_tree.info.BCyl_sizes.y, ref_tree.info.branches_curvature,
          ref_tree.info.branches_density, ref_tree.info.joints_cnt, ref_tree.info.leaves_density,
          ref_tree.info.trunk_thickness);
    std::string gen_name = selection_settings.get_string("generator_name", "simpliest");
    float imp_size = selection_settings.get_int("impostor_size", 128);
    GroveGenerationData tree_ggd;
    {
        std::vector<TreeTypeData> types = metainfoManager.get_all_tree_types();
        bool found = false;
        for (auto &type : types)
        {
            if (type.generator_name == gen_name)
            {
                tree_ggd.types = {type};
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

    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    tree_ggd.impostor_generation_params.slices_n = 8;
    tree_ggd.impostor_generation_params.quality = imp_size;
    tree_ggd.impostor_generation_params.monochrome = true;
    tree_ggd.impostor_generation_params.normals_needed = false;
    tree_ggd.impostor_generation_params.leaf_opacity = 1.0;

    GETreeGenerator::set_joints_limit(1.667 * ref_tree.info.joints_cnt);

    GeneticAlgorithm::MetaParameters mp;
    mp.best_genoms_count = selection_settings.get_int("best_results_count", mp.best_genoms_count);
    mp.initial_population_size = selection_settings.get_int("initial_population_size", mp.initial_population_size);
    mp.max_population_size = selection_settings.get_int("max_population_size", mp.max_population_size);
    mp.elite = selection_settings.get_int("elite", mp.elite);
    mp.n_islands = selection_settings.get_int("n_islands", mp.n_islands);
    mp.migration_interval = selection_settings.get_int("migration_interval", mp.migration_interval);
    mp.clone_thr = selection_settings.get_double("clone_thr", mp.clone_thr);
    mp.migration_chance = selection_settings.get_double("migration_chance", mp.migration_chance);
    mp.evolution_stat = selection_settings.get_bool("evolution_stat", mp.evolution_stat);
    mp.debug_graph = selection_settings.get_bool("debug_graph", mp.debug_graph);

    GeneticAlgorithm::ExitConditions ex_c;
    Block *ex_bl = selection_settings.get_block("exit_conditions");
    if (ex_bl)
    {
        ex_c.function_calculated = ex_bl->get_int("function_calculated", ex_c.function_calculated);
        ex_c.function_reached = ex_bl->get_double("function_reached", ex_c.function_reached);
        ex_c.generations = ex_bl->get_int("generations", ex_c.generations);
        ex_c.time_elapsed_seconds = ex_bl->get_double("time_elapsed_seconds", ex_c.time_elapsed_seconds);
    }
    int imp_max_cnt = MAX(50, MAX(mp.initial_population_size, mp.max_population_size));
    ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(imp_max_cnt, 8, false);
    ParameterList parList, bestParList;

    tree_ggd.types[0].params->write_parameter_list(parList);
    BlkManager man;
    Block b;

    man.load_block_from_file(gen_name + "_param_borders.blk", b);
    parList.load_borders_from_blk(b);
    parList.print();

    bestParList = parList;
    float best_metric = 0;
    int cnt = 0;

    bestParList = parList;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::vector<std::pair<float, ParameterList>> best_pars;
    auto func = [&](std::vector<ParameterList> &params) -> std::vector<float> {
        return generate_for_par_selection(params, imp_sim, tree_ggd, scene.heightmap, ref_tree, cnt);
    };

    GeneticAlgorithm GA;
    GA.perform(parList, mp, ex_c, func, best_pars);
    bestParList = best_pars[0].second;
    best_metric = best_pars[0].first;
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    bestParList.print();
    debug("best metric %f took %.2f seconds and %d tries to find\n", best_metric, time / 1000, cnt);

    if (ref_type && ref_type->generator_name == gen_name)
    {
        ParameterList referenceParList;
        ref_type->params->write_parameter_list(referenceParList);
        sel_quality(bestParList, referenceParList, func, 32);
    }

    //create preapred tree
    {
        LightVoxelsCube *res_voxels = gen_voxels_for_selection(ref_tree);
        tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
        GrovePacker packer;
        Tree *trees = new Tree[best_pars.size()];
        int s = sizeof(Tree);
        tree_ggd.trees_count = best_pars.size();

        auto type = tree_ggd.types[0];
        tree_ggd.types.clear();
        for (int i = 0; i < best_pars.size(); i++)
        {
            tree_ggd.types.emplace_back(type);
        }

        for (int i = 0; i < best_pars.size(); i++)
        {
            AbstractTreeGenerator *gen = GroveGenerator::get_generator(type.generator_name);
            glm::vec3 pos = glm::vec3(100 * (1 + i / 5), 0, 100 * (i % 5));
            res_voxels->fill(0);
            res_voxels->relocate(pos);

            tree_ggd.types[i].params->read_parameter_list(best_pars[i].second);
            gen->plant_tree(pos, &(tree_ggd.types[i]));
            while (gen->iterate(*res_voxels))
            {
            }
            gen->finalize_generation(trees + i, *res_voxels);
        }
        packer.add_trees_to_grove(tree_ggd, scene.grove, trees, scene.heightmap, false);
        delete[] trees;
        delete res_voxels;
    }

    GETreeGenerator::set_joints_limit(1000000);
    results.best_candidates = tree_ggd.types;
}

ParameterSelector::Results ParameterSelector::parameter_selection(TreeTypeData reference_tree_type,
                                                                  Block &selection_settings, Scene *demo_scene)
{

    Scene inner_scene;
    Scene &scene = demo_scene ? *demo_scene : inner_scene;
    scene.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(100, 100), 10);
    scene.heightmap->fill_const(0);

    float imp_size = selection_settings.get_int("impostor_size", 128);
    GroveGenerationData tree_ggd;
    tree_ggd.trees_count = 1;
    tree_ggd.types = {reference_tree_type};
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    tree_ggd.impostor_generation_params.slices_n = 8;
    tree_ggd.impostor_generation_params.quality = imp_size;
    tree_ggd.impostor_generation_params.monochrome = true;
    tree_ggd.impostor_generation_params.normals_needed = false;
    tree_ggd.impostor_generation_params.leaf_opacity = 0.33;

    ReferenceTree ref_tree;
    LightVoxelsCube *ref_voxels = new LightVoxelsCube(glm::vec3(0, 0, 0), 2.0f * reference_tree_type.params->get_tree_max_size(),
                                                      0.625f * reference_tree_type.params->get_scale_factor());
    //create reference tree
    {
        AbstractTreeGenerator *gen = GroveGenerator::get_generator(reference_tree_type.generator_name);
        ref_voxels->fill(0);
        tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
        GrovePacker packer;
        Tree single_tree;
        tree_ggd.trees_count = 1;
        gen->plant_tree(glm::vec3(0, 0, 0), &(tree_ggd.types[0]));
        while (gen->iterate(*ref_voxels))
        {
        }
        gen->finalize_generation(&single_tree, *ref_voxels);
        packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
        save_impostor_as_reference(scene.grove.impostors[1], imp_size, imp_size, "imp_ref", ref_tree.atlas);
        ref_atlas_transform(ref_tree.atlas);
        ref_tree.tex = textureManager.load_unnamed_tex(image::base_img_path + "imp_ref.png");
        ImpostorSimilarityCalc::get_tree_compare_info(scene.grove.impostors[1].impostors.back(), single_tree, ref_tree.info);
    }
    delete ref_voxels;
    Results res;
    parameter_selection_internal(selection_settings, res, scene, ref_tree, &reference_tree_type);
    return res;
}

ParameterSelector::Results ParameterSelector::parameter_selection(Block &reference_info, Block &selection_settings,
                                                                  Scene *demo_scene)
{
    int reference_images_cnt = 0;
    for (int i=0;i<reference_info.size();i++)
    {
        if (reference_info.get_name(i) == "reference_image" && reference_info.get_type(i) == Block::ValueType::BLOCK)
            reference_images_cnt++;
    }
    if (reference_images_cnt == 0)
    {
        logerr("No images selected as a reference for parameter selection");
        return Results();
    }
    
    Scene inner_scene;
    Scene &scene = demo_scene ? *demo_scene : inner_scene;
    scene.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(100, 100), 10);
    scene.heightmap->fill_const(0);
    float imp_size = selection_settings.get_int("impostor_size", 128);
    ReferenceTree ref_tree;
    ref_tree.info.BCyl_sizes = reference_info.get_vec2("BCyl_sizes",glm::vec2(0,0));
    ref_tree.info.branches_curvature = reference_info.get_double("branches_curvature",0);
    ref_tree.info.branches_density = reference_info.get_double("branches_density",0);
    ref_tree.info.joints_cnt =reference_info.get_int("joints_cnt",0);
    ref_tree.info.leaves_density = reference_info.get_double("leaves_density",0);
    ref_tree.info.trunk_thickness = reference_info.get_double("trunk_thickness",0);

    {
        TextureAtlas atl = TextureAtlas(reference_images_cnt * imp_size, imp_size, 1, 1);
        ref_tree.atlas = atl;
    }
    ref_tree.atlas.set_grid(imp_size, imp_size, false);
    PostFx ref_transform = PostFx("image_to_monochrome_impostor.fs");
    for (int i=0;i<reference_info.size();i++)
    {
        if (reference_info.get_name(i) == "reference_image" && reference_info.get_type(i) == Block::ValueType::BLOCK)
        {
            Block *ref_image_blk = reference_info.get_block(i);
            if (!ref_image_blk)
                continue;
            std::string name = ref_image_blk->get_string("image","");
            if (name == "")
                continue;
            Texture ref_raw = textureManager.load_unnamed_tex(image::base_img_path + name);
            int slice_id = ref_tree.atlas.add_tex();
            ref_tree.atlas.target_slice(slice_id, 0);
            ref_transform.use();
            ref_transform.get_shader().texture("tex", ref_raw);
            ref_transform.get_shader().uniform("wood_color",
                                               ref_image_blk->get_vec3("wood_color", glm::vec3(0.2, 0.2, 0.2)));
            ref_transform.get_shader().uniform("leaves_color",
                                               ref_image_blk->get_vec3("leaves_color", glm::vec3(0, 0.15, 0)));
            ref_transform.get_shader().uniform("background_color",
                                               ref_image_blk->get_vec3("background_color", glm::vec3(0, 0, 0)));
            ref_transform.render();

            textureManager.delete_tex(ref_raw);
        }
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    ref_atlas_transform(ref_tree.atlas);
    textureManager.save_png(ref_tree.atlas.tex(0),"reference_atlas");
    Results res;
    parameter_selection_internal(selection_settings, res, scene, ref_tree, nullptr);
    return res;
}