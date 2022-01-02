#include "sandbox.h"
#include "generation/scene_generator.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"
#include "graphics_utils/texture_manager.h"
#include "tinyEngine/image.h"
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

void bruteforce_selection(std::function<float(ParameterList &)> &f, int num_bins, int detalization_count, 
                          int detalization_depth, int num_samples, float &best_val, ParameterList &bestParams)
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
                for (int sample = 0;sample<num_samples;sample++)
                {
                    ParameterList params = grid_params(grid.bins, grid.base_set);
                    sum_metric += f(params);
                }
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
                    new_progress_bars.back().base_set.print();
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

Texture load_reference(std::string name)
{
    Texture ref_raw = textureManager.load_unnamed_tex(image::base_img_path + name);
    Texture ref = textureManager.create_unnamed(ref_raw.get_W(), ref_raw.get_H());
    PostFx ref_transform = PostFx("image_to_monochrome_impostor.fs");
    
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ref.texture, 0);
    glViewport(0, 0, ref_raw.get_W(), ref_raw.get_H());
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ref_transform.use();
    ref_transform.get_shader().texture("tex", ref_raw);
    ref_transform.get_shader().uniform("wood_color", glm::vec3(0.2,0.2,0.2));
    ref_transform.get_shader().uniform("leaves_color", glm::vec3(0,0.15,0));
    ref_transform.get_shader().uniform("background_color", glm::vec3(0,0,0));
    ref_transform.render();

    textureManager.delete_tex(ref_raw);
    textureManager.save_png(ref, "transformed_"+name);
    return ref;
}

void sandbox_main(int argc, char **argv, Scene &scene)
{
    metainfoManager.reload_all();

    scene.heightmap = new Heightmap(glm::vec3(0,0,0),glm::vec2(100,100),10);
    scene.heightmap->fill_const(0);

    Texture reference = load_reference("reference_tree_test.png");
    return;

    GroveGenerationData tree_ggd;
    tree_ggd.trees_count = 1;
    TreeTypeData type = metainfoManager.get_tree_type("simpliest_tree_default");
    AbstractTreeGenerator *gen = GroveGenerator::get_generator(type.generator_name);
    tree_ggd.types = {type};
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    LightVoxelsCube voxels = LightVoxelsCube(glm::vec3(0,0,0),type.params->get_tree_max_size(), type.params->get_scale_factor());

    ParameterList parList;
    ParameterList bestParList;
    tree_ggd.types[0].params->write_parameter_list(parList);
    parList.continuousParameters.at("branch_angle_0").min_val = 0;
    parList.continuousParameters.at("branch_angle_0").max_val = PI/2;
    parList.continuousParameters.at("branch_angle_1").min_val = 0;
    parList.continuousParameters.at("branch_angle_1").max_val = PI/2;
    parList.continuousParameters.at("branch_angle_2").min_val = 0;
    parList.continuousParameters.at("branch_angle_2").max_val = PI/2;
    //parList.continuousParameters.at("branch_angle_3").min_val = 0;
    //parList.continuousParameters.at("branch_angle_3").max_val = PI/2;
    parList.print();
    bestParList = parList;
    float best_metric = 0;
    int cnt = 0;
    std::function<float(ParameterList &)> func = [&](ParameterList &params) -> float
    {
        GrovePacker packer;
        Tree single_tree;
        GrovePacked tmp_g;
        tree_ggd.types[0].params->read_parameter_list(params);
        gen->plant_tree(glm::vec3(0,0,0),&(tree_ggd.types[0]));
        while (gen->iterate(voxels))
        {
            
        }
        gen->finalize_generation(&single_tree,voxels);
        packer.add_trees_to_grove(tree_ggd, tmp_g, &single_tree, scene.heightmap, false);
        //textureManager.save_png(tmp_g.impostors[1].atlas.tex(0),"imp0");
        logerr("generate %d",cnt);
        cnt++;
        return dot_metric(single_tree, 0.5);
    };
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    /*
    for (int i=0;i<max_iters;i++)
    {
        GrovePacker packer;
        Tree single_tree;
        GrovePacked tmp_g;
        for (auto &p : parList.continuousParameters)
        {
            if (!p.second.fixed())
                p.second.val = urand(p.second.min_val, p.second.max_val);
        }
        tree_ggd.types[0].params->read_parameter_list(parList);
        gen->plant_tree(glm::vec3(0,0,0),&(tree_ggd.types[0]));
        while (gen->iterate(voxels))
        {
            
        }
        gen->finalize_generation(&single_tree,voxels);
        packer.add_trees_to_grove(tree_ggd, tmp_g, &single_tree, scene.heightmap, false);
        //textureManager.save_png(tmp_g.impostors[1].atlas.tex(0),"imp0");

        double sum_dot = 0;
        int dot_cnt = 0;
        float dst_dot = 0.5;
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
                        sum_dot += SQR(0.5 - glm::dot(dir, ch_dir));
                        dot_cnt++;
                    }
                }
            }
        }
        float dt = sum_dot/dot_cnt;
        float metric = 1 - dt;
        logerr("%d dot metric %f %f",i, dt, metric);
        if (metric > best_metric)
        {
            best_metric = metric;
            bestParList = parList;
        }
    }
    */
    bestParList = parList;
    bruteforce_selection(func, 4, 2, 2, 1, best_metric, bestParList);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    //bestParList.print();
    logerr("best metric %f took %.2f seconds to find",best_metric, time/1000);
    tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
    GrovePacker packer;
    Tree single_tree;
        tree_ggd.types[0].params->read_parameter_list(bestParList);
        gen->plant_tree(glm::vec3(0,0,0),&(tree_ggd.types[0]));
        while (gen->iterate(voxels))
        {
            
        }
        gen->finalize_generation(&single_tree,voxels);
        packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
}