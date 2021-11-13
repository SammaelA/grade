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
    static const unsigned LEVEL_BITS = 4;
    std::vector<std::vector<PackedBranch>> branches;
    BranchCatalogue(int levels)
    {
        if (levels >= (1 << LEVEL_BITS))
        {
            logerr("Branch catalogue created with too many branch levels: %d. Max value is %d", levels, (1 << LEVEL_BITS) - 1);
            levels = (1 << LEVEL_BITS) - 1;
        }
        for (int i = 0; i < levels; i++)
        {
            branches.push_back(std::vector<PackedBranch>());
            occupied.push_back(BitVector());
            first_free.push_back(0);
        }
    }
    int levels() const
    {
        return branches.size();
    }
    PackedBranch &get(unsigned pos)
    {
        if (occupied[pos & ((1 << LEVEL_BITS) - 1)].get(pos >> LEVEL_BITS))
            return branches[pos & ((1 << LEVEL_BITS) - 1)][pos >> LEVEL_BITS];
        else
            return emptyBranch;
    }
    int add(PackedBranch &b, int level)
    {
        if (level < 0 || level >= branches.size())
            return -1;
        int pos = first_free[level];
        if (pos == branches[level].size())
        {
            //increase vector size
            first_free[level]++;
            occupied[level].push_back(true);
            branches[level].push_back(b);
        }
        else
        {
            branches[level][pos] = b;
            occupied[level].set_unsafe(pos,true);
            first_free[level] = branches[level].size(); 
            for (int i=pos;i<branches[level].size();i++)
            {
                if (!occupied[level].get_unsafe(i))
                {
                    first_free[level] = i;
                    break;
                }
            }
        }
        return ((pos) << LEVEL_BITS) + level;
    }
    int get_level_size(int level)
    {
        if (level < 0)
            level = 0;
        if (level >= branches.size())
            level = branches.size() - 1;
        int sz = 0;
        for (int i=0;i<occupied[level].size();i++)
        {
            sz += occupied[level].get_unsafe(i);
        }
        return sz;
    }
    void remove(unsigned id)
    {
        int level = id & ((1 << LEVEL_BITS) - 1);
        int pos = id >> LEVEL_BITS;
        occupied[level].set(pos, false);
        first_free[level] = MIN(first_free[level],pos);
    }
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