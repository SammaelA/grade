#include "grove_generator.h"
#include "generation/grove_generation_utils.h"
#include "planter.h"
#include "trees_preprocessor.h"
#include "tree_generators/all_generators.h"
#include <map>


void GroveGenerator::prepare_patch(GrovePrototype &prototype, 
                                   const std::vector<TreeTypeData> &treeTypesCatalogue,
                                   Heightmap &hmap,
                                   GroveMask &mask,
                                   LightVoxelsCube &voxels,
                                   Tree *trees)
{
    std::map<std::string,AbstractTreeGenerator *> generators;
    std::map<std::string,int> t_counts;
    float w = 0;
    for (auto &p : prototype.possible_types)
    {
        std::string g_name = treeTypesCatalogue[p.first].generator_name;
        generators.emplace(g_name,get_generator(g_name));
        t_counts.emplace(g_name,0);
        w += p.second;
    }

    int trees_planted = 0;
    bool generating = true;
    Planter planter = Planter(&voxels, &hmap, &mask, prototype.biome_mask,
                              glm::vec3(prototype.pos.x,0,prototype.pos.y), prototype.size,
                              1,prototype.trees_count - prototype.preplanted_trees.size(),5);
    
    for (auto &p : prototype.preplanted_trees)
    {
        AbstractTreeGenerator *gen = generators.at(treeTypesCatalogue[p.first].generator_name);
        gen->plant_tree(p.second, &(treeTypesCatalogue[p.first]));
        t_counts.at(treeTypesCatalogue[p.first].generator_name)++;
        trees_planted++;
    }
    for (auto &p : generators)
    {
      p.second->set_seed(urandi(0, RAND_MAX));
    }    
    while (trees_planted < prototype.trees_count || generating)
    {
        if (trees_planted < prototype.trees_count)
        {
            std::vector<glm::vec3> seeds = planter.get_saplings();
            
            if (!seeds.empty())
            {
                for (auto &seed : seeds)
                {
                    float f = urand(0,w);
                    int type = 0;
                    for (auto &p : prototype.possible_types)
                    {
                        if (f < p.second)
                        {
                            type = p.first;
                            break;
                        }
                        else
                            f -= p.second;
                    }

                    AbstractTreeGenerator *gen = generators.at(treeTypesCatalogue[type].generator_name);
                    gen->plant_tree(seed, &(treeTypesCatalogue[type]));
                    t_counts.at(treeTypesCatalogue[type].generator_name)++;
                    trees_planted++;
                }
            }
            else
            {
                prototype.trees_count = trees_planted;
            }
        }

        generating = false;
        for (auto &p : generators)
        {
            generating = generating || p.second->iterate(voxels);
        }
    }

    int t_count = 0;
    for (auto &p : generators)
    {
        int t_cnt = t_counts.at(p.first);
        p.second->finalize_generation(trees + t_count,voxels);
        t_count += t_cnt;
    }

    TreePreprocessor t_prep;
    Block prep_settings;
    for (int i=0;i<prototype.trees_count;i++)
    {
        t_prep.preprocess_tree(trees[i],prep_settings);
    }
}