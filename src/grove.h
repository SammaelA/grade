#include "tree.h"

struct PackedJoint
{
    glm::vec3 pos;
    float r;
};
struct PackedBranch
{
    std::vector<PackedJoint> joints;
};
struct BranchCatalogue
{
    std::vector<std::vector<PackedBranch>> branches;
    PackedBranch &get_branch(unsigned char level, unsigned pos)
    { return branches[level][pos]; }
};
struct BranchStructure
{
    unsigned char level;
    unsigned pos;
    std::vector<BranchStructure *> childBranches;
};
class Grove
{

};