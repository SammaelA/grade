#pragma once
#include "generated_tree.h"
#include "tinyEngine/utility.h"
#include "clustering/clustering.h"
#include "texture_manager.h"
#include "visualizer.h"
#include "distribution.h"
#include <math.h>
#include <algorithm>
#include "body.h"
#include <chrono>
#include "tinyEngine/save_utils/saver.h"
#include "impostor.h"
#include "terrain.h"
#include "field_2d.h"
#include "grove_generation_utils.h"

struct SyntheticTree
{
    ClusterData *trunk;
    InstanceDataArrays trunk_instance_data;//will always have vectors with length 1, used only for unification
    std::vector<ClusterData *> branches;//same size for three vectors
    std::vector<InstanceDataArrays> branches_instance_data;
    std::vector<int> joints_ns;//to which root joint this branch is attached
};
struct TransformStat
{
    Normal *phi;
    Normal *psi;
    Normal *r;
    Normal *rot_angle;
};
struct BranchStat
{
    DiscreteGeneral *typeStat = nullptr;
    TransformStat transformStat;
    bool valid = false;
};
struct RootStat
{
    std::vector<DiscreteGeneral *> branchExistanceStat;
    BranchStat selfBranchStat;
    DiscreteGeneral *childBranchesClusterStat;
    int child_branches_count;
};
struct FullStat
{
    DiscreteGeneral *rootClusterStat;
    std::vector<RootStat> rootStats;
    std::vector<BranchStat> branchStats;
    FullStat() {};
    ~FullStat();
};
class SyntheticTreeGenerator
{
public:
    SyntheticTreeGenerator(Seeder &seeder, std::vector<ClusterData> &trunks_clusters, 
                           std::vector<ClusterData> &branches_clusters, GroveGenerationData &ggd);
    void generate(Tree *_trees, int count, LightVoxelsCube *voxels);
private:
    void synt_tree_to_real(SyntheticTree &synt, Tree &t);
    glm::mat4 get_transform(Branch *base, glm::vec3 pos, BranchStat &stat);
    void make_synt_tree(SyntheticTree &synt);
    BranchStat get_branch_stat(ClusterData &cd);
    DiscreteGeneral *get_child_branches_cluster_stat(ClusterData &cd);
    void get_existance_stat(ClusterData &cd, std::vector<DiscreteGeneral *> &branchExistanceStat);
    void get_mean_and_stddev(std::vector<float> &values, double &mean, double &stddev);
    void collect_statistics();
    void add_shadow(Branch *b, glm::mat4 &transform);
    void add_planar_shadow(Tree &t);
    Seeder &seeder;
    GroveGenerationData &ggd;
    std::vector<ClusterData> &trunks_clusters;
    std::vector<ClusterData> &branches_clusters;
    Tree *trees;
    LightVoxelsCube *voxels = nullptr;
    FullStat stat;
};