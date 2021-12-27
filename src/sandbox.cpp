#include "sandbox.h"
#include "generation/scene_generator.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"
#include "graphics_utils/texture_manager.h"
void sandbox_main(int argc, char **argv, Scene &scene)
{
    metainfoManager.reload_all();

    scene.heightmap = new Heightmap(glm::vec3(0,0,0),glm::vec2(100,100),10);
    scene.heightmap->fill_const(0);
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
    parList.continuousParameters.at("branch_angle_3").min_val = 0;
    parList.continuousParameters.at("branch_angle_3").max_val = PI/2;
    parList.print();
    bestParList = parList;
    float best_metric = 0;
    int max_iters = 100;
    
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
                        sum_dot += glm::dot(dir, ch_dir);
                        dot_cnt++;
                    }
                }
            }
        }
        float dt = sum_dot/dot_cnt;
        float metric = 1 - abs(dt - dst_dot)/(dt + dst_dot);
        logerr("%d dot metric %f %f",i, dt, metric);
        if (metric > best_metric)
        {
            best_metric = metric;
            bestParList = parList;
        }
    }

    bestParList.print();
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