#pragma once

#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "glm/gtc/matrix_transform.hpp"
#include <unordered_set>
#include "volumetric_occlusion.h"
#include "tinyEngine/utility/model.h"
#include "billboard_cloud.h"
#include "tree.h"
#include "grove.h"
#include "grove_generation_utils.h"
class DebugVisualizer;
class Heightmap;

class GrovePacker
{
public:
    void pack_grove(GroveGenerationData ggd, GrovePacked &grove, DebugVisualizer &debug, 
                    ::Tree *trees_external, Heightmap *h, bool visualize_voxels);
};