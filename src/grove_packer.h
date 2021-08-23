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
#include "branch_clusterization.h"

class DebugVisualizer;
class Heightmap;
struct ClusterAdditionalData
{
    std::list<InstancedBranch>::iterator instanced_branch;
    std::vector<std::list<BillboardData>::iterator> small_billboards;
    std::vector<std::list<BillboardData>::iterator> large_billboards;
    bool is_presented = false;
};
struct ClusterPackingLayer
{
    std::vector<ClusterData> clusters;
    std::vector<ClusterAdditionalData> additional_data;
};
class GrovePacker
{
public:
    void pack_grove(GroveGenerationData ggd, GrovePacked &grove, DebugVisualizer &debug, 
                    ::Tree *trees_external, Heightmap *h, bool visualize_voxels);
    void add_trees_to_grove(GroveGenerationData ggd, GrovePacked &grove, ::Tree *trees_external, Heightmap *h);
private:
    void transform_all_according_to_root(GrovePacked &grove);

    std::vector<ClusterPackingLayer> packingLayersBranches = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrunks = {ClusterPackingLayer()};
    std::vector<ClusterPackingLayer> packingLayersTrees = {ClusterPackingLayer()};
};