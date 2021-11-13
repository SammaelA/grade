#pragma once

#include <string>
#include <fstream>
#include "tree_generators/abstract_generator.h"
#include "save_utils/blk.h"


class TreeLoaderBlk : public AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override;
    static void set_trees_directory(std::string &_trees_directory)
    {
        if (_trees_directory != trees_directory)
            trees_taken = 0;
        trees_directory = _trees_directory;
    }
    static std::string blks_base_name;
protected:
    void load_tree(GroveGenerationData ggd, ::Tree &tree_external, Heightmap &h, Block &b);
    static int trees_taken;
    static std::string trees_directory;
};
