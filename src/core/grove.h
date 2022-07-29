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
    int level = 0;
    uint type_id = 0;
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

struct CompressedTree
{
  enum NodeType
  {
    MODEL,
    BILLBOARD,
    IMPOSTOR
  };
  struct Node
  {
    Node() = default;
    Node(NodeType t, int mn, int in)
    {
      type = t;
      model_num = mn;
      instance_num = in;
    }
    NodeType type = MODEL;
    int model_num = -1;
    int instance_num = -1;
    std::vector<Node> children;
  };

  glm::vec3 pos;
  std::vector<Node> LOD_roots;  //root node for every LOD that this tree has
  AABB bbox;
  int global_id;
};
struct GroveTexturesAtlas
{
    bool atlases_valid = false;
    bool maps_valid = false;
    TextureAtlas *woodAtlas = nullptr;
    std::map<int, int> wood_tex_map;//<type_id, tex n in atlas> 
    TextureAtlas *leavesAtlas = nullptr;
    std::map<int, int> leaves_tex_map;
};
struct GrovePacked
{
    glm::vec3 center;
    std::string ggd_name;
    BranchCatalogue instancedCatalogue;

    std::list<InstancedBranch> instancedBranches;
    std::vector<BillboardCloudData> clouds;
    std::vector<ImpostorsData> impostors;
    std::vector<CompressedTree> compressedTrees;
    std::map<int, int> trees_by_global_id;
    std::vector<std::list<InstancedBranch>::const_iterator> instancedBranchesDirect;
    GroveTexturesAtlas groveTexturesAtlas;
    GrovePacked() : instancedCatalogue(MAX_BRANCH_LEVELS){};
};