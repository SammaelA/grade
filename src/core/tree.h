#pragma once
#include <list>
#include <vector>
#include <string>
#include <atomic>
#include "common_utils/utility.h"
#include "../tinyEngine/texture.h"
#include "common_utils/parameter.h"
#include "save_utils/serialization.h"

class Texture;
struct Joint;
struct Segment;
struct Branch;
struct Leaf;
struct BranchHeap;
struct LeafHeap;
struct PackedBranch;

extern std::atomic<int> br_h_cnt;

struct Segment
{
    friend class boost::serialization::access;

    glm::vec3 begin;
    glm::vec3 end;
    float rel_r_begin;
    float rel_r_end;
    std::vector<float> mults;
  
    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & begin;
      ar & end;
      ar & rel_r_begin;
      ar & rel_r_end;
      ar & mults;
    }
};
struct Joint
{
    friend class boost::serialization::access;

    Leaf *leaf = nullptr;
    glm::vec3 pos;
    std::list<Branch *> childBranches;
    short mark_A;

    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & leaf;
      ar & pos;
      ar & mark_A;

      if (Archive::is_loading::value)
      {
        std::vector<uint64_t> branch_ids;
        ar & branch_ids;
        for (uint64_t id : branch_ids)
        {
          childBranches.push_back(reinterpret_cast<Branch *>(id));//dirty hack to resolve pointer conflict 
          //it will be replaced in heap serialization
        }
      }
      else
      {
        std::vector<uint64_t> branch_ids;
        save_ids(branch_ids);
        ar & branch_ids;
      }
    }
    void save_ids(std::vector<uint64_t> &branch_ids);
};
struct Branch
{
    friend class boost::serialization::access;

    bool dead = false;
    int self_id = 0;
    int id = 0;
    ushort type_id = 0;
    short level;
    std::list<Segment> segments;
    std::list<Joint> joints;
    glm::vec4 plane_coef;//plane ax+by+cz+d = 0 len(a,b,c) = 1
    glm::vec3 center_par;
    glm::vec3 center_self;
    int mark_A;
    int mark_B;
    void deep_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void norecursive_copy(const Branch *b, BranchHeap &heap, LeafHeap *leaf_heap = nullptr);
    void transform(glm::mat4 &trans_matrix, float r_transform = 1.0);
    void pack(PackedBranch &branch);
    void mark_dead();
    static float get_r_mult(float phi, std::vector<float> &mults);

    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & dead;
      ar & self_id;
      ar & id;
      ar & type_id;
      ar & level;
      ar & segments;
      ar & joints;
      ar & plane_coef;
      ar & center_par;
      ar & center_self;
      ar & mark_A;
      ar & mark_B;
    }
};
struct Leaf
{
    friend class boost::serialization::access;

    glm::vec3 pos;
    std::vector<glm::vec3> edges;
    ushort type;

    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & pos;
      ar & edges;
      ar & type;
    }
};
struct LeafHeap
{
    friend class boost::serialization::access;

    std::list<Leaf> leaves;
    Leaf *new_leaf()
    {
        leaves.push_back(Leaf());
        return &leaves.back();
    }
    void clear_removed();
    ~LeafHeap()
    {
        leaves.clear();
    }
    LeafHeap(){};
    LeafHeap& operator=(LeafHeap&&h)
    {
        leaves = std::move(h.leaves);
    }
    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & leaves;
    }
};
struct BranchHeap
{   
    friend class boost::serialization::access;

    std::list<Branch> branches;
    Branch *new_branch()
    {
        branches.push_back(Branch());
        return &branches.back();
    }
    BranchHeap(){br_h_cnt.fetch_add(1);};
   ~BranchHeap(){br_h_cnt.fetch_add(-1);};
    BranchHeap& operator=(BranchHeap&&h)
    {
        branches = std::move(h.branches);
    }
    void clear_removed();
    private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & branches;
      if (Archive::is_loading::value)
      {
        //we need to restore pointers
        std::map<uint64_t, Branch *> branch_by_id;
        for (auto &b : branches)
          branch_by_id.emplace(b.self_id, &b);
        for (auto &b : branches)
        {
          for (auto &j : b.joints)
          {
            for (Branch *& chb: j.childBranches)
            {
              uint64_t br_id = reinterpret_cast<uint64_t>(chb);
              auto it = branch_by_id.find(br_id);
              if (it == branch_by_id.end())
              {
                logerr("cannot find branch with id %d", br_id);
              }
              else
              {
                chb = it->second;
              }
            }
          }
        }
      }
    }
};

struct TreeTypeData
{
    TreeTypeData();
    TreeTypeData(const TreeTypeData &t);
    TreeTypeData &operator=(const TreeTypeData &t);
    TreeTypeData(TreeTypeData &&t) = default;
    TreeTypeData &operator=(TreeTypeData &&t) = default;
    TreeTypeData(int id, ParameterSet *params, std::string wood_tex_name, std::string leaf_tex_name);
    ParameterSet *get_params() const;
    void set_params(ParameterSet *_params);
    ~TreeTypeData();
    int type_id;
    Texture wood;
    Texture leaf;
    int wood_id;
    int leaf_id;
    std::vector<Texture> additional_textures;
    std::string generator_name = "default";
    std::string wood_tex_name = "wood";
    std::string leaf_tex_name = "leaf";
private:
    ParameterSet *params = nullptr;
};
struct Tree
{
    std::vector<BranchHeap *> branchHeaps;
    LeafHeap *leaves = nullptr;
    glm::vec3 pos = glm::vec3(0,0,0);
    Branch *root= nullptr;
    uint id = 0;
    const TreeTypeData *type = nullptr;
    bool valid = false;
    Tree() {};
    ~Tree() {clear();};
    void clear()
    {
        if (leaves)
        {
            delete leaves;
            leaves = nullptr;
        }
        for (int i=0;i<branchHeaps.size();i++)
        {
            delete branchHeaps[i];
            branchHeaps[i] = nullptr;
        }
        branchHeaps.clear();
        root = nullptr;
        type = nullptr;
    }
};

enum Quality
{
    LOW_AS_F = 128,
    ULTRALOW = 256,
    LOW = 400,
    MEDIUM = 600,
    HIGH = 800,
    ULTRA = 1024
};
