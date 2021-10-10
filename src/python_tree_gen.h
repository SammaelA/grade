#pragma once

#include "load_tree_structure.h"

class PythonTreeGen : public TreeLoaderBlk
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) override;
};