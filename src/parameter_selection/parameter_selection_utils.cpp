#include "parameter_selection_utils.h"
#include "tree_generators/all_generators.h"
namespace ps_utils
{
    void gen_tree(LightVoxelsCube &voxels, const TreeTypeData *type, Tree *tree)
    {
        AbstractTreeGenerator *gen = get_generator(type->generator_name);
        voxels.fill(0);
        gen->plant_tree(glm::vec3(0, 0, 0), type);
        while (gen->iterate(voxels))
        {
        }
        gen->finalize_generation(tree, voxels);
        delete gen;
    }
    void gen_tree_task(int start_n, int stop_n, LightVoxelsCube *vox, const std::vector<TreeTypeData> *types, Tree *trees)
    {
        for (int i = start_n; i < stop_n; i++)
        {
            gen_tree(*vox, &((*types)[i]), trees + i);
        }
    }
};