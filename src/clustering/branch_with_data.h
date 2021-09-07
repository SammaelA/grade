#pragma once
#include "clustering.h"
#include "helpers.h"
struct ClassicStructureSimilarityParams
{
    int bwd_rotations = 18;
    float delta = 0.2;
    float light_importance = 0.4;
    float voxels_size_mult = 1/2.5;
    bool voxelized_structure = false;
    bool different_types_tolerance = true;
    int leaf_size_mult = 0;
    float structure_voxels_size_mult = 1/2.5;
    int ignore_structure_level = 1000;
    int min_clusters = 1;
    float max_individual_dist = 0.95;
    std::vector<float> weights = std::vector<float>{5000,800,40,1,0.01};
    std::vector<float> light_weights = std::vector<float>{5000,800,40,1,0.01};
    std::vector<float> r_weights = std::vector<float>{0.5,0.2,0,0,0};
    void load_from_block(Block *b);
    void load(Block *b);
};
    struct BranchWithData : public BranchClusteringData
    {
        Branch *original;
        Branch *b;
        BranchHeap branchHeap;
        LeafHeap leafHeap;
        int base_cluster_id;
        float rot_angle = 0.0;
        std::vector<int> joint_counts;
        std::vector<LightVoxelsCube *> leavesDensity;
        std::vector<LightVoxelsCube *> voxelizedStructures;

        BranchWithData(ClassicStructureSimilarityParams &clusterizationParams, Branch *_original, int levels, int _id, 
                       glm::mat4 _transform, float r_transform);
        ~BranchWithData();
        virtual void clear() override;
    };