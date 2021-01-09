#pragma once
#include "tinyEngine/utility.h"
#include "tinyEngine/utility/model.h"
#include "tinyEngine/utility/instance.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility.h"
#include "billboard_cloud_data.h"
class BillboardCloudRenderer;
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
};
struct BranchCatalogue
{
    static const unsigned LEVEL_BITS = 3;
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
        }
    }
    int levels() const
    {
        return branches.size();
    }
    PackedBranch &get(unsigned pos)
    {
        return branches[pos & ((1 << LEVEL_BITS) - 1)][pos >> LEVEL_BITS];
    }
    int add(PackedBranch &b, int level)
    {
        if (level < 0 || level >= branches.size())
            return -1;
        branches[level].push_back(b);
        return ((branches[level].size() - 1) << LEVEL_BITS) + level;
    }
    std::vector<PackedBranch> &get_level(int level)
    {
        if (level < 0)
            level = 0;
        if (level >= branches.size())
            level = branches.size() - 1;
        return branches[level];
    }
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
    std::vector<glm::mat4> transforms;
};
struct GrovePacked
{
    std::vector<BranchStructure> roots;
    BranchCatalogue instancedCatalogue;
    BranchCatalogue uniqueCatalogue;

    std::vector<InstancedBranch> instancedBranches;
    std::vector<BillboardCloudData> clouds;
    GrovePacked() : uniqueCatalogue(7), instancedCatalogue(7){};
};
class GroveRenderer
{
public:
    struct LOD
    {
        Model *m;
        BillboardCloudRenderer *cloud;
        std::vector<Instance *> instances;
    };
    void render(int lod, glm::mat4 prc);
    GroveRenderer(GrovePacked *_source, int LODs_count);

private:
    void add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch);
    std::vector<LOD> LODs;
    Shader renderer;
    Shader rendererInstancing;
    GrovePacked *source;
    Texture wood;
    Texture leaf;
};