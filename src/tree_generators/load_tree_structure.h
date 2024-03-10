#pragma once

#include <string>
#include <fstream>
#include "tree_generators/abstract_generator.h"
#include "common_utils/blk.h"


class TreeLoaderBlk : public AbstractTreeGenerator
{
public:
    virtual void plant_tree(float3 pos, const TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
    static void set_trees_directory(std::string &_trees_directory)
    {
        if (_trees_directory != trees_directory)
            trees_taken = 0;
        trees_directory = _trees_directory;
    }
    static std::string blks_base_name;
protected:
    void load_tree(const TreeTypeData *type, ::Tree &tree_external, Block &b);
    static int trees_taken;
    static std::string trees_directory;
    std::vector<std::pair<float3, const TreeTypeData *>> tree_saplings;
};
