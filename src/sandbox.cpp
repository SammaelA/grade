#include "sandbox.h"
#include "generation/scene_generator.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"
#include "graphics_utils/texture_manager.h"
#include "tinyEngine/image.h"
#include "parameter_selection/impostor_similarity.h"
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

void save_impostor_as_reference(ImpostorsData &imp, int tex_w, int tex_h, int num_slice, std::string name)
{
    Texture res(textureManager.create_unnamed(tex_w, tex_h));
    PostFx copy = PostFx("copy_arr2.fs");
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res.texture, 0);
    glViewport(0, 0, tex_w, tex_h);
    glClearColor(1,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glm::vec4 tex_transform = imp.atlas.tc_transform(imp.impostors.back().slices[num_slice].id);
    tex_transform = glm::vec4(tex_transform.z*tex_transform.x, tex_transform.w*tex_transform.y, tex_transform.x, tex_transform.y);
    //tex_transform.y += tex_transform.w;
    //tex_transform.w *= -1;
    logerr("tc %f %f %f %f", tex_transform.x, tex_transform.y, tex_transform.z, tex_transform.w);
    glm::vec3 tex = glm::vec3(0,0,0);
    imp.atlas.process_tc(imp.impostors.back().slices[num_slice].id, tex);
    int layer = tex.z;

    copy.use();
    copy.get_shader().texture("tex", imp.atlas.tex(0));
    copy.get_shader().uniform("tex_transform", tex_transform);
    copy.get_shader().uniform("layer",layer);
    copy.render();
    
    glBindTexture(imp.atlas.tex(0).type, 0);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    textureManager.save_png_directly(res, image::base_img_path + name + ".png");
    textureManager.delete_tex(res);
    glDeleteFramebuffers(1, &fbo);
}

void sandbox_main(int argc, char **argv, Scene &scene)
{
    metainfoManager.reload_all();

    scene.heightmap = new Heightmap(glm::vec3(0,0,0),glm::vec2(100,100),10);
    scene.heightmap->fill_const(0);

    GroveGenerationData tree_ggd;
    tree_ggd.trees_count = 1;
    TreeTypeData type = metainfoManager.get_tree_type("simpliest_tree_default");
    tree_ggd.types = {type};
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    tree_ggd.impostor_generation_params.slices_n = 8;
    tree_ggd.impostor_generation_params.quality = 128;
    tree_ggd.impostor_generation_params.monochrome = true;
    tree_ggd.impostor_generation_params.normals_needed = false;
    tree_ggd.impostor_generation_params.leaf_opacity = 1;

    ReferenceTree ref_tree;
    AbstractTreeGenerator *gen = GroveGenerator::get_generator(type.generator_name);
    ImpostorSimilarityCalc imp_sim = ImpostorSimilarityCalc(8, 8, false);
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
    
    /*
    parList.continuousParameters.at("branch_angle_0").min_val = 0.52;
    parList.continuousParameters.at("branch_angle_0").max_val = 0.52;
    parList.continuousParameters.at("branch_angle_1").min_val = 0.52;
    parList.continuousParameters.at("branch_angle_1").max_val = 0.52;
    parList.continuousParameters.at("branch_angle_2").min_val = 0.79;
    parList.continuousParameters.at("branch_angle_2").max_val = 0.79;*/
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
        std::vector<float> res;
        imp_sim.calc_similarity(tmp_g, ref_tree, res);
        return res[0];
        return dot_metric(single_tree, 0.5);
    };

    //create reference tree
    {
        tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
        GrovePacker packer;
        Tree single_tree;
        tree_ggd.types[0].params->read_parameter_list(parList);
        gen->plant_tree(glm::vec3(0,0,0),&(tree_ggd.types[0]));
        while (gen->iterate(voxels))
        {
            
        }
        gen->finalize_generation(&single_tree,voxels);
        packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
        save_impostor_as_reference(scene.grove.impostors[1], 256, 256, 0, "imp_ref");
        ref_tree.tex = textureManager.load_unnamed_tex(image::base_img_path + "imp_ref.png");
        ImpostorSimilarityCalc::get_tree_compare_info(scene.grove.impostors[1].impostors.back(), ref_tree.info);
    }


    bestParList = parList;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    bruteforce_selection(func, 5, 1, 1, 5, best_metric, bestParList);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    bestParList.print();
    logerr("best metric %f took %.2f seconds to find",best_metric, time/1000);


    //create preapred tree
    {
        tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
        GrovePacker packer;
        Tree single_tree;
        tree_ggd.types[0].params->read_parameter_list(bestParList);
        gen->plant_tree(glm::vec3(100,0,0),&(tree_ggd.types[0]));
        while (gen->iterate(voxels))
        {
            
        }
        gen->finalize_generation(&single_tree,voxels);
        packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
        save_impostor_as_reference(scene.grove.impostors[1], 256, 256, 0, "imp_res");
    }
}