#pragma once
#include "core/tree.h"
#include "graphics_utils/modeling.h"
#include "graphics_utils/volumetric_occlusion.h"
#include <vector>
#include <map>
#include "common_utils/blk.h"
#include "common_utils/hash.h"
#include "default_clustering_params.h"

struct ClusteringSerializationHelper
{
  template<typename T>
  struct CSHUnit
  {
    std::list<T> *base_container = nullptr;
    std::vector<typename std::list<T>::iterator> base_container_direct_access;
    void set_container(std::list<T> *_base_container)
    {
      base_container = _base_container;
      base_container_direct_access.clear();
      if (base_container)
      {
        typename std::list<T>::iterator it = base_container->begin();
        while (it != base_container->end())
        {
          base_container_direct_access.push_back(it);
          it++;
        }
      }
    }
    int pos_by_it(typename std::list<T>::iterator it)
    {
      if (!base_container)
      {
        logerr("container is not set");
        return -1;
      }
      return std::distance(base_container->begin(), it);
    }
    typename std::list<T>::iterator it_by_pos(int pos)
    {
      if (!base_container)
      {
        logerr("container is not set");
        return typename std::list<T>::iterator();
      }
      if (pos < 0 || pos >= base_container->size())
      {
        logerr("wrong pos %d, base container size %d", pos, base_container->size());
        return typename std::list<T>::iterator();
      }
      return base_container_direct_access[pos];
    }
  };

  CSHUnit<InstancedBranch> instanced_branches;
  CSHUnit<Impostor> clust_metric_impostors;
};