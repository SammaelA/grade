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

void sandbox_main(int argc, char **argv, Scene *scene)
{
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
        logerr("argv %s",argv[2]);
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