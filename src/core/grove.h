#pragma once
#include "tree.h"
#include "common_utils/bit_vector.h"
#include "billboard_cloud_data.h"
#include "limits.h"

class GroveGenerationData;
struct BBox;
struct PackedLeaf
{
    std::vector<glm::vec3> edges;
};
struct PackedJoint
{
    glm::vec3 pos;
    float r;
    PackedJoint() { r = 0; }
    PackedJoint(glm::vec3 &_pos, float _r)
    {
        pos = _pos;
        r = _r;
    }
};
struct PackedBranch
{
    std::vector<PackedJoint> joints;
    std::vector<PackedLeaf> leaves;
    std::vector<std::vector<float>> r_mults;
    int level;
    uint type_id;
    glm::vec4 plane_coef;
};
struct BranchCatalogue
{
    std::vector<std::vector<PackedBranch>> branches;
    BranchCatalogue(int levels);

    int levels() const
    {
        return branches.size();
    }
    PackedBranch &get(unsigned pos);
    int add(PackedBranch &b, int level);
    int get_level_size(int level);
    void remove(unsigned id);
private:
    PackedBranch emptyBranch;
    std::vector<BitVector> occupied;
    std::vector<int> first_free;

};
struct BranchStructure
{
    unsigned pos;
    std::vector<BranchStructure> childBranches;
    std::vector<std::pair<glm::mat4, BranchStructure>> childBranchesInstanced;
    explicit BranchStructure(unsigned _pos = 0) {pos = _pos;}
};
struct InstancedBranch
{
    std::vector<unsigned> branches;
    InstanceDataArrays IDA;
    BBox bbox;
};
struct GrovePacked
{
    glm::vec3 center;
    std::string ggd_name;
    BranchCatalogue instancedCatalogue;

    std::list<InstancedBranch> instancedBranches;
    std::vector<BillboardCloudData> clouds;
    std::vector<ImpostorsData> impostors;
    GrovePacked() : instancedCatalogue(MAX_BRANCH_LEVELS){};
};