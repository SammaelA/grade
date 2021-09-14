#include "branch_with_data.h"

void calc_joints_count(Branch *b, std::vector<int> &counts)
{
    for (Joint &j : b->joints)
    {
        counts[b->level]++;
        for (Branch *br : j.childBranches)
        {
            calc_joints_count(br, counts);
        }
    }
}

BranchWithData::~BranchWithData()
{
}

void voxelize_branch(Branch *b, LightVoxelsCube *light, int level_to);
BranchWithData::BranchWithData(ClassicStructureSimilarityParams &clusterizationParams,
                               Branch *_original, int levels, int _id, glm::mat4 _transform, float _r_transform)
{
    original = _original;
    id = _id;
    transform = _transform;
    r_transform = _r_transform;
    b = branchHeap.new_branch();
    b->deep_copy(original, branchHeap, &leafHeap);
    auto tr = glm::inverse(transform);
    b->transform(tr, _r_transform);
    for (int i = 0; i < levels; i++)
        joint_counts.push_back(0);
    if (b)
        calc_joints_count(b, joint_counts);

    glm::vec3 axis = b->joints.back().pos - b->joints.front().pos;
    glm::mat4 rot = glm::rotate(glm::mat4(1.0f), 2 * PI / clusterizationParams.bwd_rotations, axis);

    for (int i = 0; i < clusterizationParams.bwd_rotations; i++)
    {
        b->transform(rot);

        leavesDensity.push_back(new LightVoxelsCube(
            glm::vec3(0.5f * canonical_bbox().x, 0, 0),
            glm::vec3(0.5f * canonical_bbox().x, canonical_bbox().y, canonical_bbox().z),
            1 / clusterizationParams.voxels_size_mult, 1, 1, 1));

        set_occlusion(b, leavesDensity.back());
        if (clusterizationParams.voxelized_structure)
        {
            voxelizedStructures.push_back(new LightVoxelsCube(
                glm::vec3(0.5f * canonical_bbox().x, 0, 0),
                glm::vec3(0.5f * canonical_bbox().x, canonical_bbox().y, canonical_bbox().z),
                1 / clusterizationParams.structure_voxels_size_mult, 1.0f));

            voxelize_original_branch(b, voxelizedStructures.back(), 1, 1);
        }
    }
}

void BranchWithData::clear()
{
    for (int i = 0; i < leavesDensity.size(); i++)
    {
        if (leavesDensity[i])
        {
            delete leavesDensity[i];
            leavesDensity[i] = nullptr;
        }
    }
    leavesDensity.clear();
}
void ClassicStructureSimilarityParams::load_from_block(Block *b)
{
    bwd_rotations = b->get_int("bwd_rotations",bwd_rotations);
    delta = b->get_double("delta",delta);
    light_importance = b->get_double("light_importance",light_importance);
    voxels_size_mult = b->get_double("voxels_size_mult",voxels_size_mult);
    ignore_structure_level = b->get_int("ignore_structure_level",ignore_structure_level);
    min_clusters = b->get_int("min_clusters",min_clusters);
    max_individual_dist = b->get_double("max_individual_dist",max_individual_dist);
    voxelized_structure = b->get_bool("voxelized_structure",voxelized_structure);
    different_types_tolerance = b->get_bool("different_types_tolerance",different_types_tolerance);
    structure_voxels_size_mult = b->get_double("structure_voxels_size_mult",structure_voxels_size_mult);
    leaf_size_mult = b->get_int("leaf_size_mult",leaf_size_mult);
    b->get_arr("weights",weights, true);
    b->get_arr("light_weights",light_weights, true);
    b->get_arr("r_weights",r_weights, true);
}
void ClassicStructureSimilarityParams::load(Block *b)
{
    if (!b)
        return;
    
    Block &def = get_default_block();
    load_from_block(&def);
    load_from_block(b);
}