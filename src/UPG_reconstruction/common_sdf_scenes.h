#pragma once
#include "sdf_node.h"
#include "common_utils/distribution.h"

namespace upg
{
  using SceneDesc = std::pair<UPGStructure, UPGParametersRaw>;

  static SceneDesc scene_1_sphere()
  {
    SceneDesc desc;
    desc.first.s = {2,1};
    desc.second.p = {0,0,0,0.8};
    return desc;
  }

  static SceneDesc scene_8_spheres()
  {
    SceneDesc desc;
    desc.first.s = {3,3,3,2,1,2,1,3,2,1,2,1,3,3,2,1,2,1,3,2,1,2,1};
    desc.second.p = {0.6,0.6,0.6,0.5, -0.6,0.6,0.6,0.5, 0.6,-0.6,0.6,0.5, -0.6,-0.6,0.6,0.5, 
                     0.6,0  ,0  ,0.5, -0.6,0  ,0  ,0.5, 0  ,-0.6,0  ,0.5, 0   ,0.6 ,0  ,0.5};
    return desc;
  }

  static SceneDesc scene_8_rboxes()
  {
    SceneDesc desc;
    desc.first.s = {3,3,3,2,6,2,6,3,2,6,2,6,3,3,2,6,2,6,3,2,6,2,6};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3,0.1, -0.6,0.6,0.6,0.3,0.3,0.3,0.1, 0.6,-0.6,0.6,0.3,0.3,0.3,0.1, -0.6,-0.6,0.6,0.3,0.3,0.3,0.1,
                     0.6,0  ,0  ,0.3,0.3,0.3,0.1, -0.6,0  ,0  ,0.3,0.3,0.3,0.1, 0  ,-0.6,0  ,0.3,0.3,0.3,0.1, 0   ,0.6 ,0  ,0.3,0.3,0.3,0.1};
    return desc;
  }

  static SceneDesc scene_1_box()
  {
    SceneDesc desc;
    desc.first.s = {2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3};
    return desc;
  }

  static SceneDesc scene_2_boxes()
  {
    SceneDesc desc;
    desc.first.s = {3,2,4,2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3, -0.6,0.6,0.6,0.3,0.3,0.3};
    return desc;
  }

  static SceneDesc scene_4_boxes()
  {
    SceneDesc desc;
    desc.first.s = {3,3,2,4,2,4,3,2,4,2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3, -0.6,0.6,0.6,0.3,0.3,0.3, 0.6,-0.6,0.6,0.3,0.3,0.3, -0.6,-0.6,0.6,0.3,0.3,0.3};
    return desc;
  }

  static SceneDesc scene_8_boxes()
  {
    SceneDesc desc;
    desc.first.s = {3,3,3,2,4,2,4,3,2,4,2,4,3,3,2,4,2,4,3,2,4,2,4};
    desc.second.p = {0.6,0.6,0.6,0.3,0.3,0.3, -0.6,0.6,0.6,0.3,0.3,0.3, 0.6,-0.6,0.6,0.3,0.3,0.3, -0.6,-0.6,0.6,0.3,0.3,0.3,
                     0.6,0  ,0  ,0.3,0.3,0.3, -0.6,0  ,0  ,0.3,0.3,0.3, 0  ,-0.6,0  ,0.3,0.3,0.3, 0   ,0.6 ,0  ,0.3,0.3,0.3};
    return desc;
  }

  static SceneDesc scene_bubbles(int cnt_x, int cnt_z)
  {
    float3 p0(-2,1.2, 0.6);
    float3 p1(1, 1, -0.6);
    
    float base_r = std::min(abs(p1.x-p0.x)/(cnt_x+1), abs(p1.z-p0.z)/(cnt_z+1));

    std::vector<uint16_t> structure_inv;
    std::vector<float> params_inv;
    uint32_t num = 1;
    for (int i=0;i<cnt_z;i++)
    {
      for (int j=0;j<cnt_x;j++)
      {
        float3 p = float3(p0.x + (j+1+urand(-0.5,0.5))*(p1.x-p0.x)/(cnt_x+1), 
                                p1.y,
                                p0.z + (i+1+urand(-0.5,0.5))*(p1.z-p0.z)/(cnt_z+1));
        float rnd = urand(0.5,1);
        float r = rnd*base_r;
        structure_inv.push_back(1);
        structure_inv.push_back(2);
        params_inv.push_back(r);
        params_inv.push_back(p.x);
        params_inv.push_back(p.y);
        params_inv.push_back(p.z);

        uint32_t S = 1;
        while (num>= S && ((num & S) == 0))
        {
          structure_inv.push_back(3);
          S = S << 1;
        }
        num++;
      }
    }

    structure_inv.push_back(4);
    structure_inv.push_back(2);
    structure_inv.push_back(3);
    params_inv.push_back(0.5*abs(p0.x-p1.x));
    params_inv.push_back(0.5*abs(p0.y-p1.y));
    params_inv.push_back(0.5*abs(p0.z-p1.z));
    params_inv.push_back(0.5*(p0.x+p1.x));
    params_inv.push_back(0.5*(p0.y+p1.y));
    params_inv.push_back(0.5*(p0.z+p1.z));

    std::vector<uint16_t> structure = structure_inv;
    std::vector<float> params = params_inv;
    for (int i=0;i<structure.size();i++)
      structure[i] = structure_inv[structure.size()-i-1];
    for (int i=0;i<params.size();i++)
      params[i] = params_inv[params.size()-i-1];
    
    SceneDesc desc;
    desc.first.s = structure;
    desc.second.p = params;
    return desc;
  }

  static SceneDesc scene_1_cone()
  {
    SceneDesc desc;
    desc.first.s = {2,8};
    desc.second.p = {0.2,-0.1,0,0.2,sqrtf(1-0.2*0.2),1};
    return desc;
  }

  static SceneDesc scene_chair()
  {
    SceneDesc desc;
    desc.first.s = {3, 3,2,4,2,4, 3,3,2,5,2,5,3,2,5,2,5};
    desc.second.p = {0,0.0,0, 0.5,0.07,0.5, 0,0.45,-0.45, 0.5,0.45,0.07,
                     0.4,-0.3,0.4, 0.3,0.07,  -0.4,-0.3,0.4, 0.3,0.07,  0.4,-0.3,-0.4, 0.3,0.07,  -0.4,-0.3,-0.4, 0.3,0.07};
    return desc;    
  }

  static SceneDesc scene_complex_chair()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::CROTATE, SdfNodeType::CHAIR};
    desc.second.p = {1, 0, 1, 0.15, 0.15, 0.15, 0.05, 0.3, 0.5, 0.05, 0.5, 0.5};
    return desc;    
  }

  static SceneDesc scene_extrusion()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::ROTATE, SdfNodeType::EXTRUSION, SdfNodeType::CIRCLE};
    desc.second.p = {0.0, 0.0, 1.0, 0.7, 0.5};
    return desc;    
  }

  static SceneDesc scene_subtraction()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::SUBTRACT, SdfNodeType::MOVE, SdfNodeType::SPHERE, SdfNodeType::MOVE, SdfNodeType::SPHERE};
    desc.second.p = {0,0,0,0.6, 0.5,0.0,0.5,0.4};
    return desc;    
  }

  static SceneDesc scene_CSG_1()
  {
    SceneDesc desc;
    desc.first.s = {SdfNodeType::SCALE, SdfNodeType::ROTATE,
                    SdfNodeType::SUBTRACT, 
                    SdfNodeType::AND, SdfNodeType::OR, SdfNodeType::MOVE, SdfNodeType::BOX, SdfNodeType::MOVE, SdfNodeType::BOX,
                                      SdfNodeType::MOVE, SdfNodeType::SPHERE,
                    SdfNodeType::AND, SdfNodeType::MOVE, SdfNodeType::CYLINDER,
                                      SdfNodeType::MOVE, SdfNodeType::BOX};
    desc.second.p = {0.5, PI/2,0,-PI/6,
                    -1,1.5  ,0, 1,1.5,0.5,
                     1.5,1  ,0, 1.5,1,0.5,
                     0  ,0  ,0, 3,
                     0,0.5,0.5, 2,0.5,
                    -1  ,1  ,0, 0.75,1,3};
    return desc;       
  }
}