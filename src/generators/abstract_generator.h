#pragma once
#include "../parameter.h"
#include "../tree.h"
#include "../terrain.h"

class AbstractTreeGenerator
{
public:
    virtual void create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h) = 0;
};