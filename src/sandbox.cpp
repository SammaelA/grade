#include "sandbox.h"
#include "generation/scene_generator.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"

void sandbox_main(int argc, char **argv, Scene &scene)
{
    metainfoManager.reload_all();

    scene.heightmap = new Heightmap(glm::vec3(0,0,0),glm::vec2(100,100),10);
    scene.heightmap->fill_const(0);
    GrovePacker packer;
    GroveGenerationData tree_ggd;

    tree_ggd.trees_count = 1;
    TreeTypeData type = metainfoManager.get_tree_type("simpliest_tree_default");
    AbstractTreeGenerator *gen = GroveGenerator::get_generator(type.generator_name);
    tree_ggd.types = {type};
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS | GenerationTask::MODELS;
    Tree single_tree;
    LightVoxelsCube voxels = LightVoxelsCube(glm::vec3(0,0,0),type.params->get_tree_max_size(), type.params->get_scale_factor());
    gen->plant_tree(glm::vec3(0,0,0),&type);
    while (gen->iterate(voxels))
    {
        
    }
    gen->finalize_generation(&single_tree,voxels);
    packer.add_trees_to_grove(tree_ggd, scene.grove, &single_tree, scene.heightmap, false);
}