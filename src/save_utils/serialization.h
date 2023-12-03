#pragma once

/* This <boost/serialization/library_version_type.hpp> include guards against
 * an issue in boost::serialization from boost 1.74.0 that leads to compiler
 * error "'library_version_type' is not a member of 'boost::serialization'"
 * when including <boost/serialization/unordered_map.hpp>. More details
 * in ticket https://github.com/boostorg/serialization/issues/219
 */
#include <boost/serialization/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 == 74
#include <boost/serialization/library_version_type.hpp>
#endif
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/export.hpp>
#include <glm/glm.hpp>
#include "common_utils/bit_vector.h"
#include "common_utils/bbox.h"
#include "tinyEngine/model.h"
#include "common_utils/field_2d.h"
#include "graphics_utils/terrain.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "common_utils/blk.h"
#include "common_utils/bvh.h"
#include "texture_save_manager.h"
#include "graphics_utils/texture_atlas.h"
#include "tinyEngine/engine.h"

namespace boost
{
  namespace serialization
  {
    template <class Archive>
    void serialize(Archive &ar, glm::vec2 &v, const unsigned int version)
    {
      ar &v.x;
      ar &v.y;
    }

    template <class Archive>
    void serialize(Archive &ar, glm::vec3 &v, const unsigned int version)
    {
      ar &v.x;
      ar &v.y;
      ar &v.z;
    }

    template <class Archive>
    void serialize(Archive &ar, glm::vec4 &v, const unsigned int version)
    {
      ar &v.x;
      ar &v.y;
      ar &v.z;
      ar &v.w;
    }

    template <class Archive>
    void serialize(Archive &ar, glm::ivec2 &v, const unsigned int version)
    {
      ar &v.x;
      ar &v.y;
    }

    template <class Archive>
    void serialize(Archive &ar, glm::ivec3 &v, const unsigned int version)
    {
      ar &v.x;
      ar &v.y;
      ar &v.z;
    }

    template <class Archive>
    void serialize(Archive &ar, glm::ivec4 &v, const unsigned int version)
    {
      ar &v.x;
      ar &v.y;
      ar &v.z;
      ar &v.w;
    }

    template <class Archive>
    void serialize(Archive &ar, glm::mat4 &m, const unsigned int version)
    {
      ar &m[0][0];
      ar &m[0][1];
      ar &m[0][2];
      ar &m[0][3];

      ar &m[1][0];
      ar &m[1][1];
      ar &m[1][2];
      ar &m[1][3];

      ar &m[2][0];
      ar &m[2][1];
      ar &m[2][2];
      ar &m[2][3];

      ar &m[3][0];
      ar &m[3][1];
      ar &m[3][2];
      ar &m[3][3];
    }

    template <class Archive>
    void serialize(Archive &ar, BitVector &v, const unsigned int version)
    {
      ar &v.values;
      ar &v.default_value;
      ar &v.cur_size;
    }

    template <class Archive>
    void serialize(Archive &ar, BBox &b, const unsigned int version)
    {
      ar &b.position;
      ar &b.sizes;
      ar &b.a;
      ar &b.b;
      ar &b.c;
      ar &b.special;
    }

    template <class Archive>
    void serialize(Archive &ar, AABB &box, const unsigned int version)
    {
      ar &box.min_pos;
      ar &box.max_pos;
    }

    template <class Archive>
    void serialize(Archive &ar, AABB2D &box, const unsigned int version)
    {
      ar &box.min_pos;
      ar &box.max_pos;
    }

    template <class Archive>
    void serialize(Archive &ar, Sphere &s, const unsigned int version)
    {
      ar &s.pos;
      ar &s.r;
    }

    template <class Archive>
    void serialize(Archive &ar, Sphere2D &s, const unsigned int version)
    {
      ar &s.pos;
      ar &s.r;
    }

    template <class Archive>
    void serialize(Archive &ar, Mesh &m, const unsigned int version)
    {
      ar &m.positions;
      ar &m.normals;
      ar &m.colors;
      ar &m.indices;
      ar &m.mat_indicies;
      ar &m.tangents;
      ar &m.indexed;
    }

    template <class Archive>
    void serialize(Archive &ar, Model &m, const unsigned int version)
    {
      ar &m.positions;
      ar &m.normals;
      ar &m.colors;
      ar &m.indices;
      ar &m.mat_indicies;
      ar &m.tangents;
      ar &m.indexed;
      ar &m.model;
      m.update();
    }

    template <class Archive>
    void serialize(Archive &ar, Material &m, const unsigned int version)
    {
      ar &m.d;
      ar &m.illum;
      ar &m.Ka;
      ar &m.Kd;
      ar &m.Ks;
      ar & m.map_bump;
      ar & m.map_d;
      ar & m.map_Ka;
      ar & m.map_Kd;
      ar & m.map_Ks;
      ar & m.map_Ns;
      ar &m.Ni;
      ar &m.Ns;
    }

    template <class Archive>
    void serialize(Archive &ar, ComplexModel &m, const unsigned int version)
    {
      ar &m.materials;
      ar &m.models;
    }

    template <class Archive>
    void serialize(Archive &ar, Field_2d &f, const unsigned int version)
    {
      ar &f.base_val;
      ar &f.cell_size;
      ar &f.h;
      ar &f.max_val;
      ar &f.min_val;
      ar &f.pos;
      ar &f.size;
      ar &f.w;
      if (Archive::is_loading::value)
      {
        assert(f.data == nullptr);
        f.data = new float[(2 * f.w + 1) * (2 * f.h + 1)];
      }
      ar &make_array(f.data, (2 * f.w + 1) * (2 * f.h + 1));
    }

    template <class Archive>
    void serialize(Archive &ar, Heightmap &h, const unsigned int version)
    {
      ar &boost::serialization::base_object<Field_2d>(h);
    }

    template <class Archive>
    void serialize(Archive &ar, LightVoxelsCube &v, const unsigned int version)
    {
      ar &v.center;
      ar &v.voxel_size;
      ar &v.vox_x;
      ar &v.vox_y;
      ar &v.vox_z;
      ar &v.count;
      ar &v.block_size;
      ar &v.block_x;
      ar &v.block_y;
      ar &v.block_z;
      ar &v.block_cnt;
      ar &v.mip_levels;
      ar &v.mip_decrease;
      ar &v.mip_offsets;
      ar &v.mip_size_decrease;
      ar &v.mip_vox_xyz;
      if (Archive::is_loading::value)
      {
        assert(v.voxels == nullptr);
        v.voxels = new float[v.count];
      }
      ar &make_array(v.voxels, v.count);
    }

    template <class Archive>
    void serialize(Archive &ar, ::Block::DataArray &da, const unsigned int version);

    template <class Archive>
    void serialize(Archive &ar, ::Block &block, const unsigned int version);

    template <class Archive>
    void serialize(Archive &ar, ::Block::Value &v, const unsigned int version)
    {
      bool load = false;
      if (v.type == ::Block::ValueType::EMPTY)
        load = true;

      ar &v.type;
      auto &type = v.type;
      if (type == ::Block::ValueType::INT)
        ar &v.i;
      if (type == ::Block::ValueType::UINT64)
        ar &v.u;
      else if (type == ::Block::ValueType::BOOL)
        ar &v.b;
      else if (type == ::Block::ValueType::DOUBLE)
        ar &v.d;
      else if (type == ::Block::ValueType::VEC2)
        ar &v.v2;
      else if (type == ::Block::ValueType::VEC3)
        ar &v.v3;
      else if (type == ::Block::ValueType::VEC4)
        ar &v.v4;
      else if (type == ::Block::ValueType::IVEC2)
        ar &v.iv2;
      else if (type == ::Block::ValueType::IVEC3)
        ar &v.iv3;
      else if (type == ::Block::ValueType::IVEC4)
        ar &v.iv4;
      else if (type == ::Block::ValueType::MAT4)
        ar &v.m4;
      else if (type == ::Block::ValueType::STRING)
        {
          if (!v.s || load)
            v.s = new std::string("_");
          std::string s = *(v.s);
          ar &s;
          *(v.s) = s;
        }
      else if (type == ::Block::ValueType::BLOCK)
        ar &v.bl;
      else if (type == ::Block::ValueType::ARRAY)
        ar &v.a;
    }

    template <class Archive>
    void serialize(Archive &ar, ::Block::DataArray &da, const unsigned int version)
    {
      ar &da.type;
      ar &da.values;
    }

    template <class Archive>
    void serialize(Archive &ar, ::Block &block, const unsigned int version)
    {
      ar & block.names;
      ar & block.values;
    }

    template <class Archive>
    void serialize(Archive &ar, BVH::BVH_node &node, const unsigned int version)
    {
      ar & node.bbox_idx;
      ar & node.left_idx;
      ar & node.right_idx;
    }

    template <class Archive>
    void serialize(Archive &ar, BVH &bvh, const unsigned int version)
    {
      ar & bvh.simple_list;
      ar & bvh.obj_bboxes;
      ar & bvh.nodes;
      ar & bvh.root_node_idx;
      ar & bvh.added_boxes_cnt;
      ar & bvh.removed_boxes_cnt;
    }

    template <class Archive>
    void serialize(Archive &ar, Texture &t, const unsigned int version)
    {
      if (::texSaveManager)
      {
        bool loading = Archive::is_loading::value;

        ar & t.type;
        ar & t.format;
        ar & t.origin;
        ar & t.W;
        ar & t.H;
        ar & t.layers;
        ar & t.tag;
        ar & t.mip_levels;
        
        std::string token = "";
        if (!loading)
          token = ::texSaveManager->register_texture_to_save(t);
        
        ar & token;

        if (Archive::is_loading::value)
          ::texSaveManager->load_texture(token, t);
      }
    }

    template <class Archive>
    void serialize(Archive &ar, TextureAtlas &atlas, const unsigned int version)
    {
      ar & atlas.curNum;
      ar & atlas.width;
      ar & atlas.height;
      ar & atlas.layers;
      ar & atlas.gridWN;
      ar & atlas.gridHN;
      ar & atlas.isGrid;
      ar & atlas.resizable;
      ar & atlas.valid;
      ar & atlas.clearColor;
      ar & atlas.colorTex;
      ar & atlas.normalTex;
      ar & atlas.occupied;

      if (Archive::is_loading::value)
      {
        atlas.fbo = create_framebuffer();
        int prev_FBO = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
        glBindFramebuffer(GL_FRAMEBUFFER, atlas.fbo);
        atlas.bind(0,0);
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
        }
        glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
      }
    }

  } // namespace serialization
} // namespace boost