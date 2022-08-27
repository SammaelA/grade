
#include "graphics_utils/terrain.h"
#include "hydra_scene_exporter.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include <thread>
#include <iomanip>

#include <FreeImage.h>
#include <GLFW/glfw3.h>
#include "HydraAPI/hydra_api/HydraAPI.h"
#include "HydraAPI/hydra_api/HydraXMLVerify.h"
#include "HydraAPI/utils/mesh_utils.h"
#include "HydraAPI/hydra_api/LiteMath.h"

#include "core/scene.h"
#include "graphics_utils/modeling.h"
#include "tinyEngine/engine.h"

namespace hlm = LiteMath;
using pugi::xml_node;
extern GLFWwindow* g_window;
void initGLIfNeeded(int a_width, int a_height, const char* name);
namespace hydra
{
  static inline void WriteMatrix4x4(pugi::xml_node a_node, const wchar_t* a_attrib_name, float a_value[16])
  {
    std::wstringstream outStream;
    outStream << a_value[0]  << L" " << a_value[1]  << L" " << a_value[2]  << L" " << a_value[3]  << L" "
              << a_value[4]  << L" " << a_value[5]  << L" " << a_value[6]  << L" " << a_value[7]  << L" "
              << a_value[8]  << L" " << a_value[9]  << L" " << a_value[10] << L" " << a_value[11] << L" "
              << a_value[12] << L" " << a_value[13] << L" " << a_value[14] << L" " << a_value[15];

    a_node.attribute(a_attrib_name).set_value(outStream.str().c_str());
  }

void packed_branch_to_mesh(Mesh &model, GrovePacked *source, InstancedBranch &branch, 
                           int up_to_level, bool need_leaves, int wood_mat_id, int leaves_mat_id)
{
  if (branch.branches.empty())
        return;
    //clusterization process guarantees that type of all branches in instance
    //will be the same
    int type_id = branch.IDA.type_ids[0];

    uint ind_offset = model.indices.size();
    uint verts = model.colors.size();
    glm::vec4 w_tc_trans = glm::vec4(1,1,0,0);
    if (source->groveTexturesAtlas.maps_valid && source->groveTexturesAtlas.atlases_valid)
    {
      int tex_id = source->groveTexturesAtlas.wood_tex_map.at(type_id);
      glm::vec4 tctr = source->groveTexturesAtlas.woodAtlas->tc_transform(tex_id);
      w_tc_trans = tctr;
    }
    for (int id : branch.branches)
    {
        if (id < 0)
        {
            logerr("invalid id = %d", id);
            continue;//invalid id - TODO fix it
        }
        PackedBranch &b = source->instancedCatalogue.get(id);
        if (b.level <= up_to_level && !b.joints.empty())
            visualizer::packed_branch_to_model(b, &model, false, MAX(1, 3 - b.level), glm::vec2(1,0), false);
    }
      for (int i=verts;i<model.colors.size();i+=4)
      {
        model.colors[i] = (model.colors[i] + w_tc_trans.z)*w_tc_trans.x;
        model.colors[i + 1] = 1 - (model.colors[i + 1] + w_tc_trans.w)*w_tc_trans.y;
      }
    uint l_ind_offset = model.indices.size();
    uint l_verts = model.colors.size();

    verts = l_verts;
    w_tc_trans = glm::vec4(1,1,0,0);
    if (need_leaves)
    {
      if (source->groveTexturesAtlas.maps_valid && source->groveTexturesAtlas.atlases_valid)
      {
        int tex_id = source->groveTexturesAtlas.leaves_tex_map.at(type_id);
        glm::vec4 tctr = source->groveTexturesAtlas.leavesAtlas->tc_transform(tex_id);
        //logerr("type %d leaves tctr %f %f %f %f", type_id, tctr.x, tctr.y, tctr.z, tctr.w);
        w_tc_trans = tctr;
      }
        int type_slice = 0;
        for (int id : branch.branches)
        {
            if (id < 0)
            {
                logerr("invalid id = %d", id);
                continue;//invalid id - TODO fix it
            }
            PackedBranch &b = source->instancedCatalogue.get(id);
            if (!b.joints.empty())
                visualizer::packed_branch_to_model(b, &model, true, 4, glm::vec2(1,0), false);
        }
      for (int i=verts;i<model.colors.size();i+=4)
      {
        model.colors[i] = (model.colors[i] + w_tc_trans.z)*w_tc_trans.x;
        model.colors[i + 1] = 1 - (model.colors[i + 1] + w_tc_trans.w)*w_tc_trans.y;
      }
    }
    uint l_end = model.indices.size();
    l_verts = model.positions.size();
    int count = l_ind_offset - ind_offset;
    int l_count = l_end - l_ind_offset;

    if (model.positions.size() % 3 || model.normals.size() % 3 || model.indices.size() % 3
        || model.colors.size() % 4)
    {
      logerr("error creating model from packed branch");
    }
    else
    {
      model.mat_indicies = std::vector<int>(model.indices.size()/3,wood_mat_id);

      for (int i=count/3;i<model.indices.size()/3;i++)
      {
        model.mat_indicies[i] = leaves_mat_id;
      }
    }
}

void HrMesh_from_mesh(HRMeshRef &hr_mesh, Mesh &mesh, int mat_id = 0)
{
  std::vector<float> tc_2f = std::vector<float>(mesh.colors.size()/2,0);
  for (int i=0;i<mesh.colors.size()/4;i++)
  {
    tc_2f[2*i] = mesh.colors[4*i];
    tc_2f[2*i+1] = 1 - mesh.colors[4*i+1];
  } 
  auto ind_ui = std::vector<int>(mesh.indices.size(), 0);
  memcpy(ind_ui.data(), mesh.indices.data(), mesh.indices.size()*sizeof(int));
  //logerr("%d %d %d %d", mesh.positions.size(), mesh.normals.size(), tc_2f.size(), ind_ui.size());
  hrMeshOpen(hr_mesh, HR_TRIANGLE_IND3, HR_WRITE_DISCARD);
  {
    hrMeshVertexAttribPointer3f(hr_mesh, L"pos", mesh.positions.data());
    hrMeshVertexAttribPointer3f(hr_mesh, L"norm", mesh.normals.data());
    hrMeshVertexAttribPointer2f(hr_mesh, L"texcoord", tc_2f.data());
    if (mesh.mat_indicies.empty())
      hrMeshMaterialId(hr_mesh, mat_id);
    else
      hrMeshPrimitiveAttribPointer1i(hr_mesh, L"mind", mesh.mat_indicies.data());
    if (!mesh.tangents.empty())
      hrMeshVertexAttribPointer3f(hr_mesh, L"tangent", mesh.tangents.data());
    //logerr("%d %d", mesh.positions.size(), mesh.tangents.size());
    hrMeshAppendTriangles3(hr_mesh, ind_ui.size(), ind_ui.data());
  }
  hrMeshClose(hr_mesh);
}

bool export_internal2(std::string directory, Scene &scene, Block &export_settings)
{
  const int DEMO_WIDTH  = export_settings.get_int("image_width", 512);
  const int DEMO_HEIGHT = export_settings.get_int("image_height", 512);
  bool need_terrain = export_settings.get_bool("need_terrain",true);
  bool white_terrain = export_settings.get_bool("white_terrain", false);
  hrErrorCallerPlace(L"demo_02_load_obj");
  std::wstring dir = L"../../../../hydra_scenes/" + std::wstring(directory.begin(), directory.end());
  hrSceneLibraryOpen(dir.c_str(), HR_WRITE_DISCARD);
  std::wstring permanent_tex_dir = L"../../../../resources/textures/";
  std::string base_dir = "../../../../";
  std::wstring base_dir_w = L"../../../../";
  std::wstring temp_tex_dir = L"../../../../saves/";
  std::wstring g;
  g = temp_tex_dir + L"wood_atlas.png";
  HRTextureNodeRef texWood = hrTexture2DCreateFromFile(g.c_str());
  g = temp_tex_dir + L"leaves_atlas.png";
  HRTextureNodeRef texLeaf = hrTexture2DCreateFromFile(g.c_str());
  g = temp_tex_dir + L"grass_atlas.png";
  HRTextureNodeRef texGrass = hrTexture2DCreateFromFile(g.c_str());
  g = temp_tex_dir + L"leaves_atlas_alpha.png";
  HRTextureNodeRef texLeafOpacity = hrTexture2DCreateFromFile(g.c_str());
  g = temp_tex_dir + L"grass_atlas_alpha.png";
  HRTextureNodeRef texGrassOpacity = hrTexture2DCreateFromFile(g.c_str());
  HRTextureNodeRef texTerrain = hrTexture2DCreateFromFile(L"data/textures/terrain4.jpg");
  HRTextureNodeRef texTerrainWhite = hrTexture2DCreateFromFile(L"data/textures/white.bmp");

  HRTextureNodeRef cube_white[6] = {
      hrTexture2DCreateFromFile(L"data/textures/white.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/white.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/white.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/white.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/white.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/white.bmp"),
    };

  HRTextureNodeRef texEnvWhite = HRUtils::Cube2SphereLDR(cube_white);

  HRTextureNodeRef cube[6] = {
      hrTexture2DCreateFromFile(L"data/textures/clouds1/clouds1_east.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/clouds1/clouds1_west.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/clouds1/clouds1_up.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/clouds1/clouds1_down.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/clouds1/clouds1_north.bmp"),
      hrTexture2DCreateFromFile(L"data/textures/clouds1/clouds1_south.bmp"),
    };

  HRTextureNodeRef texEnv = HRUtils::Cube2SphereLDR(cube);

  HRMaterialRef mat_land = hrMaterialCreate(L"land_material");
  hrMaterialOpen(mat_land, HR_WRITE_DISCARD);
  {
    auto matNode = hrMaterialParamNode(mat_land);
    auto diff = matNode.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type").set_value(L"lambert");

    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"1 1 1");
    color.append_attribute(L"tex_apply_mode").set_value(L"multiply");
    auto texNode = hrTextureBind(white_terrain ? texTerrainWhite : texTerrain, color);

    VERIFY_XML(matNode);
  }
  hrMaterialClose(mat_land);

  HRMaterialRef mat_wood = hrMaterialCreate(L"mat_wood");
  hrMaterialOpen(mat_wood, HR_WRITE_DISCARD);
  {
    auto matNode = hrMaterialParamNode(mat_wood);
    auto diff = matNode.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type").set_value(L"lambert");

    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    color.append_attribute(L"tex_apply_mode").set_value(L"replace");
    auto texNode = hrTextureBind(texWood, color);

    VERIFY_XML(matNode);
  }
  hrMaterialClose(mat_wood);

  HRMaterialRef mat_leaf = hrMaterialCreate(L"mat_leaf");
  hrMaterialOpen(mat_leaf, HR_WRITE_DISCARD);
  {
    auto matNode = hrMaterialParamNode(mat_leaf);
  auto diff = matNode.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type").set_value(L"lambert");

    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    color.append_attribute(L"tex_apply_mode").set_value(L"multiply");
    auto texNode = hrTextureBind(texLeaf, color);

    auto opacity = matNode.append_child(L"opacity");
    opacity.append_child(L"skip_shadow").append_attribute(L"val").set_value(0);
    auto texNodeOp = hrTextureBind(texLeafOpacity, opacity);
   
    auto transl = matNode.append_child(L"translucency");
    auto colorTrans = transl.append_child(L"color");
    colorTrans.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    colorTrans.append_attribute(L"tex_apply_mode").set_value(L"multiply");

    auto texNodeTrans = hrTextureBind(texLeaf, transl);

    VERIFY_XML(matNode);
  }
  hrMaterialClose(mat_leaf);

  HRMaterialRef mat_grass = hrMaterialCreate(L"mat_grass");
  hrMaterialOpen(mat_grass, HR_WRITE_DISCARD);
  {
    auto matNode = hrMaterialParamNode(mat_grass);

    auto opacity = matNode.append_child(L"opacity");
    opacity.append_child(L"skip_shadow").append_attribute(L"val").set_value(0);
    auto texNodeOp = hrTextureBind(texGrassOpacity, opacity);
    
    auto diff = matNode.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type").set_value(L"lambert");
    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    color.append_attribute(L"tex_apply_mode").set_value(L"multiply");
    auto texNode = hrTextureBind(texGrass, diff);

    auto transl = matNode.append_child(L"translucency");
    auto colorTrans = transl.append_child(L"color");
    colorTrans.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    colorTrans.append_attribute(L"tex_apply_mode").set_value(L"multiply");
    auto texNodeTrans = hrTextureBind(texGrass, transl);

    VERIFY_XML(matNode);
  }
  hrMaterialClose(mat_grass);

  HRMeshRef terrainMeshRef = hrMeshCreate(L"terrain");
  {
    Model model;
    visualizer::heightmap_to_model(*(scene.heightmap), &model, glm::vec2(1000, 1000), glm::vec2(1000, 1000), 10, 0);
    HrMesh_from_mesh(terrainMeshRef, model, mat_land.id);
  }

  std::vector<HRMeshRef> instancedModels;
  std::vector<HRMaterialRef> instancedModelsMaterials;
  std::vector<HRTextureNodeRef> instancedModelsTextures;
  std::vector<std::vector<glm::mat4> *> instancedModelsTransforms;

  for (auto &im : scene.instanced_models)
  {
    for (int i=0;i<im.model.models.size();i++)
    {
      if (im.instances.empty())
        continue;
      Texture &tex = im.model.materials[i].map_Ka;
      std::wstring name = std::wstring(im.name.begin(), im.name.end()) + L"_" + std::to_wstring(i);
      std::wstring mat_name = L"im_mat_" + name;
      std::wstring mat_dir = base_dir_w + std::wstring(tex.origin.begin(), tex.origin.end());

      instancedModelsTextures.push_back(hrTexture2DCreateFromFile(mat_dir.c_str()));

      instancedModelsMaterials.push_back(hrMaterialCreate(mat_name.c_str()));
      HRMaterialRef &mat = instancedModelsMaterials.back();
      hrMaterialOpen(mat, HR_WRITE_DISCARD);
      {
        auto matNode = hrMaterialParamNode(mat);
        auto diff = matNode.append_child(L"diffuse");
        diff.append_attribute(L"brdf_type").set_value(L"lambert");

        auto color = diff.append_child(L"color");
        color.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
        color.append_attribute(L"tex_apply_mode").set_value(L"replace");
        auto texNode = hrTextureBind(instancedModelsTextures.back(), color);

        VERIFY_XML(matNode);
      }
      hrMaterialClose(mat);

      instancedModelsTransforms.push_back(&(im.instances));
      std::wstring model_name = L"model_" + name;
      instancedModels.push_back(hrMeshCreate(model_name.c_str()));
      HrMesh_from_mesh(instancedModels.back(), *(im.model.models[i]), mat.id);
    }
  }

  std::vector<HRMeshRef> branches;
  std::vector<InstanceDataArrays *> IDAs;

  for (auto &pb : scene.grove.instancedBranches)
  {
    Mesh m;
    packed_branch_to_mesh(m,&(scene.grove), pb, 1000, true, mat_wood.id, mat_leaf.id);
    std::wstring name = L"branch" + std::to_wstring(branches.size());
    IDAs.push_back(&(pb.IDA));
    branches.push_back(hrMeshCreate(name.c_str()));
    HrMesh_from_mesh(branches.back(), m);
  }

    std::vector<Model *> grass_models;
    std::vector<HRMeshRef> grasses;

    std::vector<int> inst_offsets;
    std::vector<int> inst_counts;

    glm::vec4 tex_transform = glm::vec4(1,1,0,0);
    Texture null = engine::textureManager->empty();
    std::string prev_base_dir = model_loader::base_path;
    model_loader::base_path = base_dir + prev_base_dir;
    int total_instances = 0;
    int type_n = 0;
    for (auto &p : scene.grass.grass_instances)
    {
        tex_transform = scene.grass.grass_textures.tc_transform(p.first);
        grass_models.push_back(model_loader::create_model_by_name(scene.grass.used_grass_types[type_n].model_name,null));
        for (int i=0;i<grass_models.back()->colors.size();i+=4)
        {
            grass_models.back()->colors[i] = tex_transform.x*(grass_models.back()->colors[i] + tex_transform.z);
            grass_models.back()->colors[i + 1] = 1 - tex_transform.y*(grass_models.back()->colors[i + 1] + tex_transform.w);
        }
        grass_models.back()->update();
        inst_offsets.push_back(total_instances);
        inst_counts.push_back(p.second.size());
        total_instances += p.second.size();
        type_n++;
    }
    model_loader::base_path = prev_base_dir;

    glm::mat4 *matrices = new glm::mat4[total_instances];
    int i=0;
    for (auto &p : scene.grass.grass_instances)
    {
        for (auto &in : p.second)
        {
            matrices[i] = glm::scale(
                          glm::rotate(glm::translate(glm::mat4(1.0f),glm::vec3(in.pos)),
                                      in.rot_y,glm::vec3(0,1,0)), glm::vec3(in.size));
            i++;
        }
    }
    
    for (auto *model : grass_models)
    {
      std::wstring name = L"grass" + std::to_wstring(grasses.size());
      grasses.push_back(hrMeshCreate(name.c_str()));
      HrMesh_from_mesh(grasses.back(), *model, mat_grass.id);
      delete model;
    }

  HRLightRef rectLight = hrLightCreate(L"my_area_light");
  hrLightOpen(rectLight, HR_WRITE_DISCARD);
  {
    pugi::xml_node lightNode = hrLightParamNode(rectLight);

    lightNode.attribute(L"type").set_value(L"area");
    lightNode.attribute(L"shape").set_value(L"rect");
    lightNode.attribute(L"distribution").set_value(L"diffuse");

    pugi::xml_node sizeNode = lightNode.append_child(L"size");

    sizeNode.append_attribute(L"half_length") = 1.0f;
    sizeNode.append_attribute(L"half_width")  = 1.0f;

    pugi::xml_node intensityNode = lightNode.append_child(L"intensity");

    intensityNode.append_child(L"color").append_attribute(L"val")      = L"1 1 1";
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 1.5f;

    VERIFY_XML(lightNode);
  }
  hrLightClose(rectLight);

  HRLightRef directLight = hrLightCreate(L"my_direct_light");
  hrLightOpen(directLight, HR_WRITE_DISCARD);
  {
    pugi::xml_node lightNode = hrLightParamNode(directLight);

    lightNode.attribute(L"type").set_value(L"directional");
    lightNode.attribute(L"shape").set_value(L"point");

    pugi::xml_node sizeNode = lightNode.append_child(L"size");

    sizeNode.append_child(L"inner_radius").append_attribute(L"val") = 1.0f;
    sizeNode.append_child(L"outer_radius").append_attribute(L"val") = 2.0f;

    pugi::xml_node intensityNode = lightNode.append_child(L"intensity");

    intensityNode.append_child(L"color").append_attribute(L"val") = L"0.4285 0.3095 0.262";
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 5.0f * PI;
  }
  hrLightClose(directLight);

  HRLightRef sky = hrLightCreate(L"sky");
  hrLightOpen(sky, HR_WRITE_DISCARD);
  {
    auto lightNode = hrLightParamNode(sky);
    lightNode.attribute(L"type").set_value(L"sky");
	  lightNode.attribute(L"distribution").set_value(L"map");
    auto intensityNode = lightNode.append_child(L"intensity");
    intensityNode.append_child(L"color").append_attribute(L"val").set_value(L"1 1 1");
    intensityNode.append_child(L"multiplier").append_attribute(L"val").set_value(L"1");

	  auto texNode = hrTextureBind(white_terrain ? texEnvWhite : texEnv, intensityNode.child(L"color"));
    //texNode.append_attribute(L"input_gamma").set_value(1.0f);
    //texNode.append_attribute(L"input_alpha").set_value(L"rgb");

	  VERIFY_XML(lightNode);
  }
  hrLightClose(sky);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // camera
  //
  glm::vec3 camera_pos = export_settings.get_vec3("camera_pos",glm::vec3(0,200,200));
  glm::vec3 camera_look_at = export_settings.get_vec3("camera_look_at",glm::vec3(0,0,0));
  glm::vec3 camera_up = export_settings.get_vec3("camera_up", glm::vec3(0,1,0));
  HRCameraRef camRef = hrCameraCreate(L"my camera");
  
  hrCameraOpen(camRef, HR_WRITE_DISCARD);
  {
    xml_node camNode = hrCameraParamNode(camRef);
    
    camNode.append_child(L"fov").text().set(L"90");
    camNode.append_child(L"nearClipPlane").text().set(L"0.01");
    camNode.append_child(L"farClipPlane").text().set(L"1000.0");
    
    std::wstring camera_up_s = std::to_wstring(camera_up.x)+L" "+std::to_wstring(camera_up.y)+L" "+std::to_wstring(camera_up.z);
    camNode.append_child(L"up").text().set(camera_up_s.c_str());
    std::wstring camera_pos_s = std::to_wstring(camera_pos.x)+L" "+std::to_wstring(camera_pos.y)+L" "+std::to_wstring(camera_pos.z);
    camNode.append_child(L"position").text().set(camera_pos_s.c_str());
    
    std::wstring camera_look_at_s = std::to_wstring(camera_look_at.x)+L" "+std::to_wstring(camera_look_at.y)+L" "+std::to_wstring(camera_look_at.z);
    camNode.append_child(L"look_at").text().set(camera_look_at_s.c_str());
  
    VERIFY_XML(camNode);
  }
  hrCameraClose(camRef);
  
  // set up render settings
  //
  HRRenderRef renderRef = hrRenderCreate(L"HydraModern"); // "HydraModern" is our main renderer name
  hrRenderEnableDevice(renderRef, 0, true);
  
  hrRenderOpen(renderRef, HR_WRITE_DISCARD);
  {
    pugi::xml_node node = hrRenderParamNode(renderRef);
    
    node.append_child(L"width").text()  = DEMO_WIDTH;
    node.append_child(L"height").text() = DEMO_HEIGHT;
    
    node.append_child(L"method_primary").text()   = L"pathtracing";
    node.append_child(L"method_secondary").text() = L"pathtracing";
    node.append_child(L"method_tertiary").text()  = L"pathtracing";
    node.append_child(L"method_caustic").text()   = L"pathtracing";
    
    node.append_child(L"trace_depth").text()      = 4;
    node.append_child(L"diff_trace_depth").text() = 3;
    node.append_child(L"maxRaysPerPixel").text()  = export_settings.get_int("max_rays_per_pixel", 64);
    node.append_child(L"qmc_variant").text()      = (HYDRA_QMC_DOF_FLAG | HYDRA_QMC_MTL_FLAG | HYDRA_QMC_LGT_FLAG); // enable all of them, results to '7'
  }
  hrRenderClose(renderRef);
  
  // create scene
  //
  HRSceneInstRef scnRef = hrSceneCreate(L"my scene");
  
  const float DEG_TO_RAD = float(3.14159265358979323846f) / 180.0f;
  
  hrSceneOpen(scnRef, HR_WRITE_DISCARD);
  {
    auto mind = hlm::scale4x4(hlm::float3(1,1,1));
    if (need_terrain)
    {
      hrMeshInstance(scnRef, terrainMeshRef, mind.L());
    }
    for (int i=0; i< instancedModels.size();i++)
    {
      for (auto &mat : *(instancedModelsTransforms[i]))
      {
        const float mat_vals2[16] = {mat[0][0],mat[1][0],mat[2][0],mat[3][0],
                                    mat[0][1],mat[1][1],mat[2][1],mat[3][1],
                                    mat[0][2],mat[1][2],mat[2][2],mat[3][2],
                                    mat[0][3],mat[1][3],mat[2][3],mat[3][3]};
        auto m = hlm::float4x4(mat_vals2);
        hrMeshInstance(scnRef, instancedModels[i], m.L());
      }
    }
    for (int i=0;i<branches.size();i++)
    {
      for (auto &mat : IDAs[i]->transforms)
      {
        const float mat_vals2[16] = {mat[0][0],mat[1][0],mat[2][0],mat[3][0],
                                    mat[0][1],mat[1][1],mat[2][1],mat[3][1],
                                    mat[0][2],mat[1][2],mat[2][2],mat[3][2],
                                    mat[0][3],mat[1][3],mat[2][3],mat[3][3]};
        auto m = hlm::float4x4(mat_vals2);
        hrMeshInstance(scnRef, branches[i], m.L());
      }
    }
    
    for (int i=0;i<inst_counts.size();i++)
    {
      for (int in = inst_offsets[i];in < inst_offsets[i] + inst_counts[i];in++)
      {
        auto &mat = matrices[in];
        const float mat_vals[16] = {mat[0][0],mat[1][0],mat[2][0],mat[3][0],
                                    mat[0][1],mat[1][1],mat[2][1],mat[3][1],
                                    mat[0][2],mat[1][2],mat[2][2],mat[3][2],
                                    mat[0][3],mat[1][3],mat[2][3],mat[3][3]};
        auto m = hlm::float4x4(mat_vals);
        hrMeshInstance(scnRef, grasses[i], m.L());
        //break;
      }
      logerr("added %d grass instances", inst_counts[i]);
    }
    
    delete matrices;
    hrLightInstance(scnRef, sky, mind.L());
    auto mres = hlm::mul(hlm::rotate4x4Z(30.0f * DEG_TO_RAD), hlm::translate4x4({0.0f, 3000.0f, 0.0f}));
    hrLightInstance(scnRef, directLight, mres.L());
  }
  hrSceneClose(scnRef);
  
  hrFlush(scnRef, renderRef, camRef);
  
  //////////////////////////////////////////////////////// opengl
  std::vector<int32_t> image(DEMO_WIDTH*DEMO_HEIGHT);
  initGLIfNeeded(DEMO_WIDTH,DEMO_HEIGHT, "load 'obj.' file demo");
  glViewport(0,0,DEMO_WIDTH,DEMO_HEIGHT);
  //////////////////////////////////////////////////////// opengl
  
  while (true)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    
    HRRenderUpdateInfo info = hrRenderHaveUpdate(renderRef);
    
    if (info.haveUpdateFB)
    {
      auto pres = std::cout.precision(2);
      std::cout << "rendering progress = " << info.progress << "% \r"; std::cout.flush();
      std::cout.precision(pres);
      
      hrRenderGetFrameBufferLDR1i(renderRef, DEMO_WIDTH, DEMO_HEIGHT, &image[0]);
  
      //////////////////////////////////////////////////////// opengl
      glDisable(GL_TEXTURE_2D);
      glDrawPixels(DEMO_WIDTH, DEMO_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, &image[0]);
      
      glfwSwapBuffers(g_window);
      glfwPollEvents();
      //////////////////////////////////////////////////////// opengl
    }
    
    if (info.finalUpdate)
      break;
  }
  
  std::wstring demo_dir = dir + std::wstring(L"/demo.png");
  std::string demo_copy_dir = export_settings.get_string("demo_copy_dir","");
  if (demo_copy_dir != "")
  {
    demo_dir = L"../../../../"+std::wstring(demo_copy_dir.begin(), demo_copy_dir.end())+ std::wstring(L".png"); 
  }
  hrRenderSaveFrameBufferLDR(renderRef, demo_dir.c_str());
  return true;
}

void get_default_settings(Block &b)
{
  b.set_vec3("camera_pos", glm::vec3(0,200,200));
  b.set_vec3("camera_look_at",glm::vec3(0,0,0));
  b.set_vec3("camera_up", glm::vec3(0,1,0));
  b.set_int("max_rays_per_pixel", 64);
  b.set_int("image_width", 512);
  b.set_int("image_height", 512);
  b.set_string("demo_copy_dir","");
  b.set_bool("need_terrain",true);
  b.set_bool("white_terrain", false);
}
}
