#pragma once
#include "tree.h"
#include "billboard_cloud_data.h"

#include <list>
#include <vector>

struct BaseBranchClusteringData
{
    glm::mat4 transform = glm::mat4(1.0f);
    float r_transform = 1;
    bool can_be_center = true;
};

struct BranchClusteringData
{
    friend class boost::serialization::access;

    int base_cluster_id;
    int id;  
    unsigned short tree_type;
    glm::mat4 transform;
    float r_transform;
    bool can_be_center;
    glm::vec3 sizes;
    virtual void clear() {};
    virtual ~BranchClusteringData() {clear();}
    void set_base(BaseBranchClusteringData &base)
    {
        transform = base.transform;
        r_transform = base.r_transform;
        can_be_center = base.can_be_center;
    }

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & base_cluster_id;
      ar & id;
      ar & tree_type;
      ar & transform;
      ar & r_transform;
      ar & can_be_center;
      ar & sizes;
    }
};

struct AdditionalClusterDataArrays
{
  friend class boost::serialization::access;

  std::vector<BranchClusteringData *> clustering_data;

private:
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & clustering_data;
  }
};
struct ClusterData
{   
    friend class boost::serialization::access;

    long id = -1;//cluster id
    int base_pos = 0;
    Branch *base = nullptr;
    InstanceDataArrays IDA;
    AdditionalClusterDataArrays ACDA;
    bool is_valid();
    ClusterData();

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & id;
      ar & base_pos;
      ar & IDA;
      ar & ACDA;
    }
};