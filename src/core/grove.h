#pragma once
#include "tree.h"
#include "common_utils/bit_vector.h"
#include "billboard_cloud_data.h"
#include "limits.h"
#include "save_utils/serialization.h"

struct BBox;
struct PackedLeaf
{
    friend class boost::serialization::access;

    std::vector<float3> edges;

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & edges;
    }
};
struct PackedJoint
{
    friend class boost::serialization::access;

    float3 pos;
    float r;
    PackedJoint() { r = 0; }
    PackedJoint(float3 &_pos, float _r)
    {
        pos = _pos;
        r = _r;
    }

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & pos;
      ar & r;
    }
};
struct PackedBranch
{
    friend class boost::serialization::access;

    std::vector<PackedJoint> joints;
    std::vector<PackedLeaf> leaves;
    std::vector<std::vector<float>> r_mults;
    int level = 0;
    uint type_id = 0;
    float4 plane_coef;

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & joints;
      ar & leaves;
      ar & r_mults;
      ar & level;
      ar & type_id;
      ar & plane_coef;
    }
};
struct BranchCatalogue
{
    friend class boost::serialization::access;

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

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & branches;
      ar & occupied;
      ar & first_free;
    }
};

struct InstancedBranch
{
  friend class boost::serialization::access;

    std::vector<unsigned> branches;
    InstanceDataArrays IDA;
    BBox bbox;
    int id = -1;
private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & branches;
      ar & IDA;
      ar & bbox;
      ar & id;
    }
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
    friend class boost::serialization::access;

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

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & type;
      ar & model_num;
      ar & instance_num;
      ar & children;
    }
  };

  friend class boost::serialization::access;

  float3 pos;
  std::vector<Node> LOD_roots;  //root node for every LOD that this tree has
  AABB bbox;
  int global_id;

private:
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & pos;
    ar & LOD_roots;
    ar & bbox;
    ar & global_id;
  }
};
struct GroveTexturesAtlas
{
    friend class boost::serialization::access;

    bool atlases_valid = false;
    bool maps_valid = false;
    TextureAtlas *woodAtlas = nullptr;
    std::map<int, int> wood_tex_map;//<type_id, tex n in atlas> 
    TextureAtlas *leavesAtlas = nullptr;
    std::map<int, int> leaves_tex_map;
    void clear()
    {
      if (woodAtlas)
        delete woodAtlas;
      woodAtlas = nullptr;
      if (leavesAtlas)
        delete leavesAtlas;
      leavesAtlas = nullptr;
    }
    GroveTexturesAtlas& operator=(GroveTexturesAtlas &&s) = delete;
    ~GroveTexturesAtlas(){clear();}

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & atlases_valid;
      ar & maps_valid;
      ar & wood_tex_map;
      ar & leaves_tex_map;
      ar & woodAtlas;
      ar & leavesAtlas;
    }
};
struct GrovePacked
{
    friend class boost::serialization::access;

    BranchCatalogue instancedCatalogue;
    std::list<InstancedBranch> instancedBranches;
    std::vector<BillboardCloudData> clouds;
    std::vector<ImpostorsData> impostors;
    std::vector<CompressedTree> compressedTrees;
    std::map<int, int> trees_by_global_id;
    std::vector<std::list<InstancedBranch>::const_iterator> instancedBranchesDirect;
    GroveTexturesAtlas groveTexturesAtlas;
    GrovePacked() : instancedCatalogue(MAX_BRANCH_LEVELS){};
    ~GrovePacked() {clear();}
    void clear()
    {
      instancedCatalogue = BranchCatalogue(MAX_BRANCH_LEVELS);
      instancedBranches.clear();
      clouds.clear();
      impostors.clear();
      compressedTrees.clear();
      trees_by_global_id.clear();
      instancedBranchesDirect.clear();
      groveTexturesAtlas.clear();
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & instancedCatalogue;
      ar & instancedBranches;
      //ar & clouds;
      //ar & impostors;
      ar & compressedTrees;
      ar & trees_by_global_id;
      ar & groveTexturesAtlas;

      instancedBranchesDirect.clear();
      for (auto it = instancedBranches.begin(); it != instancedBranches.end(); it++)
        instancedBranchesDirect.push_back(it);
    }
};