#pragma once
#include "tree_generators/abstract_generator.h"
#include "generation/grove_generation_utils.h"

struct GrovePrototype
{
    friend class boost::serialization::access;

    int trees_count;
    float2 pos;//pos of center
    float2 size;//size from center, BBox is [pos - (size.x,0,size.z)] - [pos + size]
    std::vector<std::pair<int,float>> possible_types;//type id and chance to create tree of that type
    std::vector<std::pair<int, float3>> preplanted_trees;
    GroveMask *biome_mask = nullptr;//should be deleted manually

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & trees_count;
      ar & pos;
      ar & size;
      ar & possible_types;
      ar & preplanted_trees;
      ar & biome_mask;
    }
};
class GroveGenerator
{
public:
    void prepare_patch(GrovePrototype &prototype, 
                       const std::vector<TreeTypeData> &treeTypesCatalogue,
                       Heightmap &hmap,
                       GroveMask &mask,
                       LightVoxelsCube &voxels,
                       Tree *trees);
};