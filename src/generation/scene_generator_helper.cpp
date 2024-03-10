#include "scene_generator_helper.h"


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

  void get_AABB_list_from_instance(ComplexModel &cm, float4x4 transform, std::vector<AABB> &boxes, float _column_size, float inflation_q)
  {
    for (Model *m : cm.models)
    {
      if (!m || m->positions.size() % 3 || _column_size < 0.1)
      {
        logerr("cannot get AABB list - malformed model");
        continue;
      }
      float3 min_pos = float3(1e9, 1e9, 1e9);
      float3 max_pos = float3(-1e9, -1e9, -1e9);

      //find AABB for whole model
      std::vector<float3> positions = std::vector<float3>(m->positions.size() / 3, float3(0, 0, 0));
      for (int i = 0; i < m->positions.size(); i += 3)
      {
        float4 p = float4(m->positions[i], m->positions[i + 1], m->positions[i + 2], 1);
        float3 p3 = to_float3(transform * p);
        positions[i / 3] = p3;
        min_pos = min(min_pos, p3);
        max_pos = max(max_pos, p3);
      }
      std::vector<int3> indices = std::vector<int3>(m->indices.size()/3, int3(0,0,0));
      for (int i=0;i<m->indices.size();i+=3)
      {
        indices[i / 3] = int3(m->indices[i], m->indices[i+1], m->indices[i+2]);
      }
      //inflate bbox
      float3 full_sz = inflation_q * (max_pos - min_pos);
      float3 center = 0.5f * (max_pos + min_pos);
      max_pos = center + 0.5f * full_sz;
      min_pos = center - 0.5f * full_sz;

      if (full_sz.x*full_sz.y*full_sz.z < 1e-9)
      {
        logerr("cannot get AABB list - malformed or too small model");
        continue;
      }
      AABB bounding = AABB(min_pos, max_pos);
      
      //we should be more precise in x,z coordinates than in y
      float2 column_size = float2(_column_size,_column_size);
      int2 columns_cnt = int2(ceil(full_sz.x/column_size.x), ceil(full_sz.z/column_size.y));

      float csz = std::max(column_size.x, column_size.y);

      int i = 0;
      while (i < indices.size())
      {
        float3 p1 = positions[indices[i].x];
        float3 p2 = positions[indices[i].y];
        float3 p3 = positions[indices[i].z];
        float d1 = length(p1 - p2);
        float d2 = length(p1 - p3);
        float d3 = length(p2 - p3);
        if (d1 > 2 * csz || d2 > 2 * csz || d3 > 2 * csz)
        {
          positions.push_back(0.5f * p1 + 0.5f * p3);
          positions.push_back(0.5f * p2 + 0.5f * p3);
          positions.push_back(0.5f * p1 + 0.5f * p2);
          int ni[6] = {indices[i].x, indices[i].y, indices[i].z, (int)(positions.size() - 3), (int)(positions.size() - 2), (int)(positions.size() - 1)};
          indices[i] = int3(ni[0], ni[5], ni[3]);
          indices.push_back(int3(ni[1], ni[4], ni[5]));
          indices.push_back(int3(ni[2], ni[3], ni[4]));
          indices.push_back(int3(ni[3], ni[5], ni[4]));
        }
        else
        {
          i++;
        }
      }
      column_size = float2(full_sz.x/columns_cnt.x, full_sz.z/columns_cnt.y);
      
      std::vector<int2> column_positions = std::vector<int2>(positions.size(), int2(0,0));
      for (int i=0;i<positions.size();i++)
      {
        column_positions[i] = to_int2(float2(positions[i].x - bounding.min_pos.x, positions[i].z - bounding.min_pos.z)/column_size); 
      }

      std::vector<float2> columns = std::vector<float2>(columns_cnt.x*columns_cnt.y, float2(1,-1));
      for (int i = 0;i<indices.size();i++)//assume that everything is triangles
      {
        int2 max_column = max(max(column_positions[indices[i].x], column_positions[indices[i].y]),column_positions[indices[i].z]);
        int2 min_column = min(min(column_positions[indices[i].x], column_positions[indices[i].y]),column_positions[indices[i].z]);
      
        float max_h = std::max(std::max(positions[indices[i].x].y, positions[indices[i].y].y), positions[indices[i].z].y);
        float min_h = std::max(std::max(positions[indices[i].x].y, positions[indices[i].y].y), positions[indices[i].z].y);

        for (int x = min_column.x; x<= max_column.x;x++)
        {
          for (int y = min_column.y; y<= max_column.y;y++)
          {
            //logerr("ind %d %d %d",i,x,y);
            columns[y*columns_cnt.x + x].x = std::min(columns[y*columns_cnt.x + x].x, min_h);
            columns[y*columns_cnt.x + x].y = std::max(columns[y*columns_cnt.x + x].y, max_h);
          }
        }
      }
      for (int x = 0; x< columns_cnt.x;x++)
      {
        for (int y=0;y<columns_cnt.y;y++)
        {
          float2 hh = columns[y*columns_cnt.x + x];
          if (hh.x < hh.y)
          {
            float h_range = inflation_q*(hh.y - hh.x);
            float h0 = 0.5*(hh.y + hh.x);
            hh.x = h0 - 0.5*h_range;
            hh.y = h0 + 0.5*h_range;
            AABB column = AABB(float3(bounding.min_pos.x + x*column_size.x, hh.x, bounding.min_pos.z + y*column_size.y),
                              float3(bounding.min_pos.x + (x+1)*column_size.x, hh.y, bounding.min_pos.z + (y+1)*column_size.y));
            boxes.push_back(column);
          }
        }
      }
    }
  }

  bool is_terrain(float4 world_pos_type)
  {
    return (int)(world_pos_type.w) == PIXEL_TYPE_TERRAIN;
  }
};