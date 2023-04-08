#pragma once
#include <glm/glm.hpp>
#include "voxelization.h"
#include <array>
#include <fstream>
#include <vector>

template <typename VoxelType>
class VoxelArray
{
public:
  struct TrilinearSampleData
  {
    bool out_of_array = true;
    std::array<float, 8> Qs;
    std::array<glm::ivec3, 8> tcs;
    std::array<VoxelType, 8> values;
  };
  VoxelArray() = default;
  VoxelArray(glm::vec3 _p0, glm::vec3 _p1, glm::ivec3 _vox_count, VoxelType _default_value,
             VoxelType *external_data = nullptr)
  {
    p0 = _p0;
    p1 = _p1;
    vox_count = _vox_count;
    default_value = _default_value;

    assert(p0.x < p1.x && p0.y < p1.y && p0.z < p1.z);
    assert(vox_count.x > 0 && vox_count.y > 0 && vox_count.z > 0);

    sz = p1 - p0;
    voxel_size = glm::vec3(sz.x / vox_count.x, sz.y / vox_count.y, sz.z / vox_count.z);
    total_vox_count = vox_count.x * vox_count.y * vox_count.z;

    if (external_data)
    {
      data = external_data;
      has_external_data = true;
    }
    else
    {
      data = new VoxelType[total_vox_count];
      has_external_data = false;
    }
  }
  ~VoxelArray()
  {
    if (data && !has_external_data)
      delete data;
  }

  VoxelType get(glm::vec3 pos)
  {
    glm::ivec3 voxel = pos_to_v(pos);
    if (in_array(voxel))
      return data[v_to_i(voxel)];
    else
      return default_value;
  };
  VoxelType get_trilinear(glm::vec3 pos)
  {
    glm::vec3 fv = (pos - p0) / voxel_size - glm::vec3(0.5);
    glm::ivec3 voxel = glm::ivec3(fv);
    glm::vec3 dp(fv.x - voxel.x, fv.y - voxel.y, fv.z - voxel.z);
    int x = voxel.x;
    int y = voxel.y;
    int z = voxel.z;

    bool in_1 = in_array(voxel);
    bool in_2 = in_array(voxel + glm::ivec3(1,1,1));

    if (in_1 && in_2)
    {
      #define C(i, j, k) data[v_to_i(x + i, y + j, z + k)]
      VoxelType c00 = C(0, 0, 0) * (1 - dp.x) + C(1, 0, 0) * dp.x;
      VoxelType c01 = C(0, 0, 1) * (1 - dp.x) + C(1, 0, 1) * dp.x;
      VoxelType c10 = C(0, 1, 0) * (1 - dp.x) + C(1, 1, 0) * dp.x;
      VoxelType c11 = C(0, 1, 1) * (1 - dp.x) + C(1, 1, 1) * dp.x;

      VoxelType c0 = c00 * (1 - dp.y) + c10 * dp.y;
      VoxelType c1 = c01 * (1 - dp.y) + c11 * dp.y;

      VoxelType c = c0 * (1 - dp.z) + c1 * dp.z;

      return c;
    }
    else if (in_1 || in_2)
    {
      #define C2(i, j, k) data[v_to_i_safe(x + i, y + j, z + k)]
      VoxelType c00 = C2(0, 0, 0) * (1 - dp.x) + C2(1, 0, 0) * dp.x;
      VoxelType c01 = C2(0, 0, 1) * (1 - dp.x) + C2(1, 0, 1) * dp.x;
      VoxelType c10 = C2(0, 1, 0) * (1 - dp.x) + C2(1, 1, 0) * dp.x;
      VoxelType c11 = C2(0, 1, 1) * (1 - dp.x) + C2(1, 1, 1) * dp.x;

      VoxelType c0 = c00 * (1 - dp.y) + c10 * dp.y;
      VoxelType c1 = c01 * (1 - dp.y) + c11 * dp.y;

      VoxelType c = c0 * (1 - dp.z) + c1 * dp.z;

      return c;
    }
    else
    {
      return default_value;
    }
  }
  void set(glm::vec3 pos, const VoxelType &value)
  {
    glm::ivec3 voxel = pos_to_v(pos);
    if (in_array(voxel))
      data[v_to_i(voxel)] = value;
  };
  void set_circle(glm::vec3 pos, float r, const VoxelType &value)
  {
    glm::ivec3 bbox_p0 = pos_to_v(pos + glm::vec3(-r));
    bbox_p0.x = voxelization::clamp(bbox_p0.x, 0, vox_count.x - 1);
    bbox_p0.y = voxelization::clamp(bbox_p0.y, 0, vox_count.y - 1);
    bbox_p0.z = voxelization::clamp(bbox_p0.z, 0, vox_count.z - 1);
    glm::ivec3 bbox_p1 = pos_to_v(pos + glm::vec3(r));
    bbox_p1.x = voxelization::clamp(bbox_p1.x, 0, vox_count.x - 1);
    bbox_p1.y = voxelization::clamp(bbox_p1.y, 0, vox_count.y - 1);
    bbox_p1.z = voxelization::clamp(bbox_p1.z, 0, vox_count.z - 1);

    for (unsigned k = bbox_p0.z; k <= bbox_p1.z; k++)
    {
      for (unsigned j = bbox_p0.y; j <= bbox_p1.y; j++)
      {
        for (unsigned i = bbox_p0.x; i <= bbox_p1.x; i++)
        {
          glm::vec3 p = v_to_pos(i, j, k) - pos;
          if (glm::dot(p, p) <= r * r)
            data[v_to_i(i, j, k)] = value;
        }
      }
    }
  }

  unsigned get_total_vox_count() 
  {
    return total_vox_count;
  }
  void set_direct(unsigned id, const VoxelType &value)
  {
    data[id] = value;
  }

  bool write_to_binary_file(std::string file_name)
  {
    std::ofstream fs(file_name, std::ios::out | std::ios::binary); 
    fs.write((const char *)data, total_vox_count*sizeof(VoxelType));
    fs.close();
    return !fs.fail();
  }
  bool read_from_binary_file(std::string file_name)
  {
    std::ifstream input(file_name, std::ios::binary);
    std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(input), {});
    if (total_vox_count*sizeof(VoxelType) == buffer.size())
    {
      memcpy(data, buffer.data(), buffer.size());
      return true;
    }
    else
      return false;
  }

private:
  inline glm::ivec3 pos_to_v(glm::vec3 pos)
  {
    return glm::ivec3((pos - p0) / voxel_size);
  }
  inline glm::vec3 v_to_pos(glm::ivec3 voxel)
  {
    return v_to_pos(voxel.x, voxel.y, voxel.z);
  }
  inline glm::vec3 v_to_pos(unsigned vx, unsigned vy, unsigned vz)
  {
    return p0 + voxel_size * glm::vec3(vx + 0.5, vy + 0.5, vz + 0.5);
  }
  inline bool in_array(glm::ivec3 voxel)
  {
    return voxel.x >= 0 && voxel.y >= 0 && voxel.z >= 0 &&
           voxel.x < vox_count.x && voxel.y < vox_count.y && voxel.z < vox_count.z;
  }
  inline unsigned v_to_i(glm::ivec3 voxel)
  {
    return v_to_i(voxel.x, voxel.y, voxel.z);
  }
  inline unsigned v_to_i(unsigned vx, unsigned vy, unsigned vz)
  {
    return vz * vox_count.x * vox_count.y + vy * vox_count.x + vx;
  }
  inline unsigned v_to_i_safe(unsigned vx, unsigned vy, unsigned vz)
  {
    return voxelization::clamp(vz, 0, vox_count.z - 1) * vox_count.x * vox_count.y + 
           voxelization::clamp(vy, 0, vox_count.y - 1) * vox_count.x + 
           voxelization::clamp(vx, 0, vox_count.x - 1);
  }

  VoxelType *data = nullptr;
  VoxelType default_value;
  glm::vec3 p0, p1, sz;     // corners of bounding box of voxel array
  glm::ivec3 vox_count;     // size of array in voxels i.e. (64,64,64)
  glm::vec3 voxel_size;     // size of each voxel
  unsigned total_vox_count; // number of all voxels in array
  bool has_external_data = false;
};