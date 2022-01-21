#include "sandbox.h"
#include "generation/scene_generator.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"
#include "graphics_utils/texture_manager.h"
#include "tinyEngine/image.h"
#include "parameter_selection/impostor_similarity.h"
#include "parameter_selection/genetic_algorithm.h"
#include "tree_generators/GE_generator.h"
#include <thread>
#include <chrono>
struct BS_Grid
{
    std::vector<std::pair<int, int>> bins;
    ParameterList base_set;
};

ParameterList grid_params(std::vector<std::pair<int, int>> bins, ParameterList &base_set)
{
    ParameterList params = base_set;
                        int k=0;
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
                            float step = (p.second.max_val - p.second.min_val)/bins[k].second;
                            float min = p.second.min_val;
                            p.second.val = min + (bins[k].first + urand())*step; 
                            p.second.min_val = min + (bins[k].first)*step;
                            p.second.max_val = min + (bins[k].first + 1)*step;
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

    for (int layer = 0;layer<detalization_depth;layer++)
    {
        if (layer == detalization_depth - 1)
            detalization_count = 1;

        std::vector<BS_Grid> new_progress_bars;
        for (auto &grid : progress_bars)
        {
            std::vector<std::pair<float, std::vector<std::pair<int, int>>>> cur_best;
            int i=0;
            while (i < grid.bins.size())
            {
                float sum_metric = 0;
                std::vector<ParameterList> params; 
                for (int sample = 0;sample<num_samples;sample++)
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
                    for (int j = 0;j<cur_best.size();j++)
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
                        for (int j=0;j<i;j++)
                        {
                            grid.bins[j].first = 0;
                            if (grid.bins[j].second > 1 && t<0)
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
                    //new_progress_bars.back().base_set.print();
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
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ref_transform.use();
    ref_transform.get_shader().texture("tex", ref_raw);
    ref_transform.get_shader().uniform("wood_color", glm::vec3(0.2,0.2,0.2));
    ref_transform.get_shader().uniform("leaves_color", glm::vec3(0,0.15,0));
    ref_transform.get_shader().uniform("background_color", glm::vec3(0,0,0));
    ref_transform.render();

    textureManager.delete_tex(ref_raw);
    glDeleteFramebuffers(1, &fbo);
    return ref;
}

void save_impostor_as_reference(ImpostorsData &imp, int tex_w, int tex_h, std::string name, TextureAtlas &ref_atl)
{
    {
        TextureAtlas atl = TextureAtlas(8*tex_w, tex_h, 1, 1);
        ref_atl = atl;
    }
    ref_atl.set_grid(tex_w,tex_h,false);
    PostFx copy = PostFx("copy_arr2.fs");
    for (int num_slice = 0; num_slice < 8; num_slice++)
    {
        glm::vec4 tex_transform = imp.atlas.tc_transform(imp.impostors.back().slices[num_slice].id);
        tex_transform = glm::vec4(tex_transform.z*tex_transform.x, tex_transform.w*tex_transform.y, tex_transform.x, tex_transform.y);
        glm::vec3 tex = glm::vec3(0,0,0);
        imp.atlas.process_tc(imp.impostors.back().slices[num_slice].id, tex);
        int layer = tex.z;
        
        int tex_id = ref_atl.add_tex();
        ref_atl.target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex", imp.atlas.tex(0));
        copy.get_shader().uniform("tex_transform", tex_transform);
        copy.get_shader().uniform("layer",(float)layer);
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
    for (int i=0;i<samples;i++)
    {
        parLists.push_back(parList);
        referenceParLists.push_back(referenceParList);
    }
    res = metric(parLists);
    ref = metric(referenceParLists);
    float mn1=0, mn2=0, mx1=0, mx2=0, av=0, dev=0;
    vector_stat(res, &mn1, &mx1, &av, &dev);
    debug("result stat [%.2f - %.2f] av=%.3f dev=%.4f\n", mn1, mx1, av, dev);
    vector_stat(ref, &mn2, &mx2, &av, &dev);
    debug("reference stat [%.2f - %.2f] av=%.3f dev=%.4f\n", mn2, mx2, av, dev);
    float q = MAX(MIN(mx1, mx2) - MAX(mn1, mn2),1e-6)/MAX(MAX(mx1, mx2) - MIN(mn1, mn2),1e-6);
    debug("quality: %.4f %f %f %f %f\n", q, MIN(mx1, mx2), MAX(mn1, mn2),MAX(mx1, mx2), MIN(mn1, mn2));

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
    for (int l=0;l<atl.layers_count();l++)
    {
        atl_tmp.target(l,0);
        gauss.use();
        gauss.get_shader().texture("tex", atl.tex(0));
        gauss.get_shader().uniform("tex_transform", glm::vec4(0,0,1,1));
        gauss.get_shader().uniform("layer",(float)l);
        gauss.get_shader().uniform("pass", 0);
        gauss.get_shader().uniform("tex_size_inv", glm::vec2(1.0f/sizes.x, 1.0f/sizes.y));
        gauss.get_shader().uniform("slice_size", slice_size);
        gauss.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl_tmp.tex(0), "atlass_gauss_0");
    for (int l=0;l<atl.layers_count();l++)
    {
        atl.target(l,0);
        gauss.use();
        gauss.get_shader().texture("tex", atl_tmp.tex(0));
        gauss.get_shader().uniform("tex_transform", glm::vec4(0,0,1,1));
        gauss.get_shader().uniform("layer",l);
        gauss.get_shader().uniform("pass", 1);
        gauss.get_shader().uniform("tex_size_inv", glm::vec2(1.0f/sizes.x, 1.0f/sizes.y));
        gauss.get_shader().uniform("slice_size", slice_size);
        gauss.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl.tex(0), "atlass_gauss_1");

    PostFx sil_fill = PostFx("silhouette_fill.fs");
    for (int l=0;l<atl.layers_count();l++)
    {
        atl_tmp.target(l,0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl.tex(0));
        sil_fill.get_shader().uniform("tex_transform", glm::vec4(0,0,1,1));
        sil_fill.get_shader().uniform("layer",(float)l);
        sil_fill.get_shader().uniform("radius", 8);
        sil_fill.get_shader().uniform("dir_threshold", 4);
        sil_fill.get_shader().uniform("tex_size_inv", glm::vec2(1.0f/sizes.x, 1.0f/sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", slice_size);
        sil_fill.render();
    } 
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    for (int l=0;l<atl.layers_count();l++)
    {
        atl.target(l,0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl_tmp.tex(0));
        sil_fill.get_shader().uniform("tex_transform", glm::vec4(0,0,1,1));
        sil_fill.get_shader().uniform("layer",(float)l);
        sil_fill.get_shader().uniform("radius", 4);
        sil_fill.get_shader().uniform("dir_threshold", 6);
        sil_fill.get_shader().uniform("tex_size_inv", glm::vec2(1.0f/sizes.x, 1.0f/sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", slice_size);
        sil_fill.render();
    } 
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    for (int l=0;l<atl.layers_count();l++)
    {
        atl_tmp.target(l,0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl.tex(0));
        sil_fill.get_shader().uniform("tex_transform", glm::vec4(0,0,1,1));
        sil_fill.get_shader().uniform("layer",(float)l);
        sil_fill.get_shader().uniform("radius", 4);
        sil_fill.get_shader().uniform("dir_threshold", 6);
        sil_fill.get_shader().uniform("tex_size_inv", glm::vec2(1.0f/sizes.x, 1.0f/sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", slice_size);
        sil_fill.render();
    } 
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //textureManager.save_png(atl_tmp.tex(0), "atlass_gauss_2");

    PostFx sil_sharp = PostFx("silhouette_sharpen.fs");
    for (int l=0;l<atl.layers_count();l++)
    {
        atl.target(l,0);
        sil_sharp.use();
        sil_sharp.get_shader().texture("tex", atl_tmp.tex(0));
        sil_sharp.get_shader().uniform("tex_transform", glm::vec4(0,0,1,1));
        sil_sharp.get_shader().uniform("layer",(float)l);
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
    gen->plant_tree(glm::vec3(0,0,0),type);
    while (gen->iterate(voxels)){}
    gen->finalize_generation(tree,voxels);
}
void gen_tree_task(int start_n, int stop_n, LightVoxelsCube *vox, std::vector<TreeTypeData> *types, Tree *trees)
{
    for (int i=start_n;i<stop_n;i++)
    {
        //logerr("gen tree %d in %d %d", i, start_n, stop_n);
        gen_tree(*vox, &((*types)[i]), trees + i);
    }
}

void sandbox_main(int argc, char **argv, Scene &scene)
{
    metainfoManager.reload_all();

    scene.heightmap = new Heightmap(glm::vec3(0,0,0),glm::vec2(100,100),10);
    scene.heightmap->fill_const(0);

    float imp_size = 128;
    GroveGenerationData tree_ggd;
    tree_ggd.trees_count = 1;
    //TreeTypeData type = metainfoManager.get_tree_type("simpliest_tree_default");
    TreeTypeData type = metainfoManager.get_tree_type("small_oak");
    tree_ggd.types = {type};
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    tree_ggd.impostor_generation_params.slices_n = 8;
    tree_ggd.impostor_generation_params.quality = imp_size;
    tree_ggd.impostor_generation_params.monochrome = true;
    tree_ggd.impostor_generation_params.normals_needed = false;
    tree_ggd.impostor_generation_params.leaf_opacity = 0.33;

    ReferenceTree ref_tree;
    GeneticAlgorithm::MetaParameters mp;
    int imp_max_cnt = MAX(mp.initial_population_size, mp.max_population_size);
    ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(imp_max_cnt, 8, false);
    LightVoxelsCube voxels = LightVoxelsCube(glm::vec3(0,0,0), 2.0f*type.params->get_tree_max_size(), 0.625f*type.params->get_scale_factor());
    ParameterList parList, bestParList, referenceParList;

    tree_ggd.types[0].params->write_parameter_list(parList);
    BlkManager man;
    Block b;

    //man.load_block_from_file("simpliest_gen_param_borders.blk", b);
    man.load_block_from_file("ge_gen_param_borders.blk", b);
    parList.load_borders_from_blk(b);
    parList.print();

    bestParList = parList;
    referenceParList = parList;
    float best_metric = 0;
    int cnt = 0;

    auto func = [&](std::vector<ParameterList> &params) -> std::vector<float>
    {
        if (params.empty())
            return std::vector<float>();
        GrovePacker packer;
        Tree *trees = new Tree[params.size()];
        GrovePacked tmp_g;
        for (int i=tree_ggd.types.size();i<params.size();i++)
        {
            tree_ggd.types.push_back(tree_ggd.types[0]);
        }
        tree_ggd.trees_count = params.size();
        for (int i=0;i<params.size();i++)
        {
            tree_ggd.types[i].params->read_parameter_list(params[i]);
        }

        int num_threads = MIN(8, params.size());
        int step = ceil(params.size()/(float)num_threads);
        LightVoxelsCube **thr_voxels = new LightVoxelsCube*[num_threads];
        for (int i=0;i<num_threads;i++)
        {
            thr_voxels[i] = new LightVoxelsCube(glm::vec3(0,0,0), 2.0f*type.params->get_tree_max_size(), 0.625f*type.params->get_scale_factor());
        }
        std::vector<std::thread> threads;
        for (int i=0;i<num_threads;i++)
        {
            int start_n = step*i;
            int stop_n = MIN(step*(i+1), params.size());
            threads.push_back(std::thread(&gen_tree_task, start_n, stop_n, thr_voxels[i],&(tree_ggd.types), trees));
        }

        for (auto &t : threads)
        {
            t.join();
        }
        packer.add_trees_to_grove(tree_ggd, tmp_g, trees, scene.heightmap, false);
        ref_atlas_transform(tmp_g.impostors[1].atlas);
        //textureManager.save_png(tmp_g.impostors[1].atlas.tex(0),"imp"+std::to_string(cnt));
        //logerr("generate %d",cnt);
        cnt+=params.size();
        std::vector<float> res;
        imp_sim.calc_similarity(tmp_g, ref_tree, res, trees);
        delete[] trees;
        for (int i=0;i<num_threads;i++)
        {
            delete thr_voxels[i];
        }
        delete[] thr_voxels;
        return res;
    };

    //create reference tree
    {
        AbstractTreeGenerator *gen = GroveGenerator::get_generator(type.generator_name);
        voxels.fill(0);
        tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
        GrovePacker packer;
        Tree single_tree;
        tree_ggd.trees_count = 1;
        tree_ggd.types[0].params->read_parameter_list(referenceParList);
        gen->plant_tree(glm::vec3(0,0,0),&(tree_ggd.types[0]));
        while (gen->iterate(voxels))
        {
            
        }
        gen->finalize_generation(&single_tree,voxels);
        packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
        save_impostor_as_reference(scene.grove.impostors[1], imp_size, imp_size, "imp_ref", ref_tree.atlas);
        ref_atlas_transform(ref_tree.atlas);
        ref_tree.tex = textureManager.load_unnamed_tex(image::base_img_path + "imp_ref.png");
        ImpostorSimilarityCalc::get_tree_compare_info(scene.grove.impostors[1].impostors.back(), single_tree, ref_tree.info);
        GETreeGenerator::set_joints_limit(1.667*ref_tree.info.joints_cnt);
    }


    bestParList = parList;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::vector<std::pair<float,ParameterList>> best_pars;
    if (true)
    {
        //bruteforce_selection(func, 4, 1, 1, 4, best_metric, bestParList);
        GeneticAlgorithm GA;
        GA.perform(parList, GeneticAlgorithm::MetaParameters(), GeneticAlgorithm::ExitConditions(), func, best_pars);
        bestParList = best_pars[0].second;
        best_metric = best_pars[0].first;
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        float time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        bestParList.print();
        debug("best metric %f took %.2f seconds and %d tries to find\n",best_metric, time/1000, cnt);
        //calculate_selection_quality
        sel_quality(bestParList, referenceParList, func, 32);
    }
    else
    {
       best_pars.push_back(std::pair<float,ParameterList>(0, bestParList));
       best_pars.push_back(std::pair<float,ParameterList>(0, bestParList));
    }
    //create preapred tree
    {
        tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
        GrovePacker packer;
        Tree *trees = new Tree[best_pars.size()];
        int s = sizeof(Tree);
        tree_ggd.trees_count = best_pars.size();

        auto type = tree_ggd.types[0];
        tree_ggd.types.clear();
        for (int i=0;i<best_pars.size();i++)
        {
            tree_ggd.types.emplace_back(type);
        }
        //tree_ggd.types = std::vector<TreeTypeData>(best_pars.size(), type);
        for (int i=0;i<best_pars.size();i++)
        {
            AbstractTreeGenerator *gen = GroveGenerator::get_generator(type.generator_name);
            glm::vec3 pos = glm::vec3(100*(1 + i/5),0,100*(i%5));
            voxels.fill(0);
            voxels.relocate(pos);

            tree_ggd.types[i].params->read_parameter_list(best_pars[i].second);
            gen->plant_tree(pos,&(tree_ggd.types[i]));
            while (gen->iterate(voxels))
            {
                
            }
            gen->finalize_generation(trees+i,voxels);
        }
        packer.add_trees_to_grove(tree_ggd, scene.grove, trees, scene.heightmap, false);
        delete[] trees;
        //TextureAtlas atl;
        //save_impostor_as_reference(scene.grove.impostors[1], imp_size, imp_size, "imp_res", atl);
    }
}