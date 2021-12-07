#include "scene_generator_helper.h"
using namespace glm;

namespace SceneGenHelper
{
  uint64_t pack_id(unsigned _empty, unsigned category, unsigned type, unsigned id)
  {
    //[_empty][category][type][id]
    //[8][8][16][32]
    return ((_empty & ((1UL << 8) - 1UL)) << 56) |
           ((category & ((1UL << 8) - 1UL)) << 48) |
           ((type & ((1UL << 16) - 1UL)) << 32) |
           ((id & ((1UL << 32) - 1UL)) << 0);
  }

  void unpack_id(uint64_t packed_id, unsigned &_empty, unsigned &category, unsigned &type, unsigned &id)
  {
    //[_empty][category][type][id]
    //[8][8][16][32]
    _empty = (packed_id >> 56) & ((1UL << 8) - 1UL);
    category = (packed_id >> 48) & ((1UL << 8) - 1UL);
    type = (packed_id >> 32) & ((1UL << 16) - 1UL);
    id = (packed_id >> 0) & ((1UL << 32) - 1UL);
  }

  void get_AABB_list_from_instance(Model *m, glm::mat4 transform, std::vector<AABB> &boxes, int max_count, float inflation_q)
  {
    if (!m || m->positions.size() % 3 || max_count < 1)
    {
      logerr("cannot get AABB list - malformed model");
      return;
    }
    boxes.clear();
    vec3 min_pos = vec3(1e9, 1e9, 1e9);
    vec3 max_pos = vec3(-1e9, -1e9, -1e9);

    //find AABB for whole model
    std::vector<vec3> positions = std::vector<vec3>(m->positions.size() / 3, vec3(0, 0, 0));
    for (int i = 0; i < m->positions.size(); i += 3)
    {
      vec4 p = vec4(m->positions[i], m->positions[i + 1], m->positions[i + 2], 1);
      vec3 p3 = vec3(transform * p);
      positions[i / 3] = p3;
      min_pos = min(min_pos, p3);
      max_pos = max(max_pos, p3);
    }

    //inflate bbox
    vec3 full_sz = inflation_q * (max_pos - min_pos);
    vec3 center = 0.5f * (max_pos + min_pos);
    max_pos = center + 0.5f * full_sz;
    min_pos = center - 0.5f * full_sz;

    if (full_sz.x*full_sz.y*full_sz.z < 1e-9)
    {
      logerr("cannot get AABB list - malformed or too small model");
      return;
    }
    if (max_count == 1)
    {
      boxes.push_back(AABB(min_pos, max_pos));
      return;
    }
    //find main == longest axis
    int main_axis = 1;
    int sa[2] = {0,2};
    if (full_sz.x > full_sz.y && full_sz.x > full_sz.z)
    {
      main_axis = 0;
      sa[0] = 1;
      sa[1] = 2;
    }
    else if (full_sz.y > full_sz.x && full_sz.y > full_sz.z)
    {
      main_axis = 1;
      sa[0] = 0;
      sa[1] = 2;
    }
    else 
    {
      main_axis = 2;
      sa[0] = 0;
      sa[1] = 1;
    }
    //logerr("box center %f %f %f full_size %f %f %f %d %d %d",center.x, center.y, center.z, full_sz.x, full_sz.y, full_sz.z,
    //main_axis,sa[0],sa[1]);
    //slice AABB of object into max_count slices along the main axis and then find AABB for each slice
    float step = full_sz[main_axis]/max_count;
    float eps = 0.01*step;

    std::vector<vec3> maxes = std::vector<vec3>(max_count, vec3(-1e9, -1e9, -1e9));
    std::vector<vec3> mins = std::vector<vec3>(max_count, vec3(1e9, 1e9, 1e9)); 

    //set starting maxes and mins
    for (int i=0;i<max_count;i++)
    {
      maxes[i][sa[0]] = center[sa[0]];
      maxes[i][sa[1]] = center[sa[1]];
      maxes[i][main_axis] = min_pos[main_axis] + (i+1)*step;

      mins[i][sa[0]] = center[sa[0]];
      mins[i][sa[1]] = center[sa[1]];
      mins[i][main_axis] = min_pos[main_axis] + (i)*step;
    }

    //find maxes and mins for each slice
    
    for (auto &pos : positions)
    {
      int slice = CLAMP((pos[main_axis] - min_pos[main_axis])/step, 0, max_count-1);

      maxes[slice][sa[0]] = MAX(maxes[slice][sa[0]], pos[sa[0]]);
      maxes[slice][sa[1]] = MAX(maxes[slice][sa[1]], pos[sa[1]]);
      mins[slice][sa[0]] = MIN(mins[slice][sa[0]], pos[sa[0]]);
      mins[slice][sa[1]] = MIN(mins[slice][sa[1]], pos[sa[1]]);
    }

    //inflate them
    for (int i=0;i<max_count;i++)
    {
      vec3 center = 0.5f*(maxes[i] + mins[i]);
      vec3 full_size = maxes[i] - mins[i];
      full_size[sa[0]] *= inflation_q;
      full_size[sa[1]] *= inflation_q;

      maxes[i] = center + 0.5f*full_size;
      mins[i] = center - 0.5f*full_size;
    }

    AABB prev_bbox, cur_bbox;
    prev_bbox.min_pos = vec3(0,0,0);
    prev_bbox.max_pos = vec3(0,0,0);
    //transform them to bboxes
    for (int i=0;i<max_count;i++)
    {
      cur_bbox = AABB(mins[i], maxes[i]);

      if (prev_bbox.min_pos[sa[0]] == prev_bbox.max_pos[sa[0]] || 
          prev_bbox.min_pos[sa[1]] == prev_bbox.max_pos[sa[1]])
      {
        prev_bbox = cur_bbox;
      }
      else
      {
        if (abs(prev_bbox.min_pos[sa[0]] - cur_bbox.min_pos[sa[0]]) < eps &&
            abs(prev_bbox.min_pos[sa[1]] - cur_bbox.min_pos[sa[1]]) < eps &&
            abs(prev_bbox.max_pos[sa[0]] - cur_bbox.max_pos[sa[0]]) < eps &&
            abs(prev_bbox.max_pos[sa[1]] - cur_bbox.max_pos[sa[1]]) < eps)
        {
          prev_bbox.min_pos[sa[0]] = MIN(prev_bbox.min_pos[sa[0]], cur_bbox.min_pos[sa[0]]);
          prev_bbox.min_pos[sa[1]] = MIN(prev_bbox.min_pos[sa[1]], cur_bbox.min_pos[sa[1]]);
          prev_bbox.max_pos[sa[0]] = MAX(prev_bbox.max_pos[sa[0]], cur_bbox.max_pos[sa[0]]);
          prev_bbox.max_pos[sa[1]] = MAX(prev_bbox.max_pos[sa[1]], cur_bbox.max_pos[sa[1]]);
          
          prev_bbox.max_pos[main_axis] = cur_bbox.max_pos[main_axis];
        }
        else if (cur_bbox.min_pos[sa[0]] == cur_bbox.max_pos[sa[0]] || 
                 cur_bbox.min_pos[sa[1]] == cur_bbox.max_pos[sa[1]])
        {
          prev_bbox.max_pos[main_axis] = cur_bbox.max_pos[main_axis];   
        }
        else
        {
                    logerr("created new AABB for object %.1f %.1f %.1f -- %.1f %.1f %.1f", 
            prev_bbox.min_pos.x, prev_bbox.min_pos.y, prev_bbox.min_pos.z,
            prev_bbox.max_pos.x, prev_bbox.max_pos.y, prev_bbox.max_pos.z);
          boxes.push_back(prev_bbox); 
          prev_bbox = cur_bbox;

        }
      }
    }
    boxes.push_back(prev_bbox);
  }
};