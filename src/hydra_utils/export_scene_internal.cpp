
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
namespace hlm = HydraLiteMath;
using pugi::xml_node;
extern GLFWwindow* g_window;

  static inline void WriteMatrix4x4(pugi::xml_node a_node, const wchar_t* a_attrib_name, float a_value[16])
  {
    std::wstringstream outStream;
    outStream << a_value[0]  << L" " << a_value[1]  << L" " << a_value[2]  << L" " << a_value[3]  << L" "
              << a_value[4]  << L" " << a_value[5]  << L" " << a_value[6]  << L" " << a_value[7]  << L" "
              << a_value[8]  << L" " << a_value[9]  << L" " << a_value[10] << L" " << a_value[11] << L" "
              << a_value[12] << L" " << a_value[13] << L" " << a_value[14] << L" " << a_value[15];

    a_node.attribute(a_attrib_name).set_value(outStream.str().c_str());
  }

void initGLIfNeeded(int a_width, int a_height, const char* name);

void heightmap_to_simple_mesh(Heightmap &h, SimpleMesh &mesh)
{
    glm::vec2 size = h.get_size();
    glm::vec3 pos = h.get_pos();
    glm::vec2 step = glm::vec2(10,10);
    SimpleMesh &flat_terrain = mesh;
            int x = (2*size.x/step.x) + 1;
            int y = (2*size.y/step.y) + 1;
            glm::vec2 range = h.get_height_range();
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    int ind = flat_terrain.vPos.size()/3;
                    glm::vec3 terr_pos = glm::vec3(pos.x - size.x + step.x*i,0,pos.z - size.y + step.y*j);
                glm::vec3 terr_pos1 = glm::vec3(pos.x - size.x + step.x*(i+1),0,pos.z - size.y + step.y*j);
                glm::vec3 terr_pos2 = glm::vec3(pos.x - size.x + step.x*(i),0,pos.z - size.y + step.y*(j+1));

                terr_pos.y = h.get_height(terr_pos);
                terr_pos1.y = h.get_height(terr_pos1);
                terr_pos2.y = h.get_height(terr_pos2);
                glm::vec3 n = glm::normalize(glm::cross(terr_pos1 - terr_pos,terr_pos2 - terr_pos));

                flat_terrain.vPos.push_back(terr_pos.x);
                flat_terrain.vPos.push_back(terr_pos.y);
                flat_terrain.vPos.push_back(terr_pos.z);
                //flat_terrain.vPos.push_back(1);

                flat_terrain.vNorm.push_back(n.x);
                flat_terrain.vNorm.push_back(n.y);
                flat_terrain.vNorm.push_back(n.z);

                flat_terrain.vTexCoord.push_back(0.5*(i % 2));
                flat_terrain.vTexCoord.push_back(0.5*(j % 2));

                if (i != x - 1 && j != y - 1)
                {
                    flat_terrain.triIndices.push_back(ind);
                    flat_terrain.triIndices.push_back(ind + y + 1);
                    flat_terrain.triIndices.push_back(ind + 1);
 
                    flat_terrain.triIndices.push_back(ind);
                    flat_terrain.triIndices.push_back(ind + y);
                    flat_terrain.triIndices.push_back(ind + y + 1);
                }
                }
            }
        }
void packed_branch_to_simple_mesh(SimpleMesh &mesh, GrovePacked *source, InstancedBranch &branch, 
                                  int up_to_level, bool need_leaves, int wood_mat_id, int leaves_mat_id)
{
  if (branch.branches.empty())
        return;
    //clusterization process guarantees that type of all branches in instance
    //will be the same
    Model model;
    //uint type = source->instancedCatalogue.get(branch.branches.front()).type_id;

    Visualizer v = Visualizer();
    uint ind_offset = model.indices.size();
    uint verts = model.positions.size();
    for (int id : branch.branches)
    {
        if (id < 0)
        {
            logerr("invalid id = %d", id);
            continue;//invalid id - TODO fix it
        }
        PackedBranch &b = source->instancedCatalogue.get(id);
        if (b.level <= up_to_level && !b.joints.empty())
            v.packed_branch_to_model(b, &model, false, up_to_level);
    }
    uint l_ind_offset = model.indices.size();
    uint l_verts = model.positions.size();

    verts = l_verts;
    if (need_leaves)
    {
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
                v.packed_branch_to_model(b, &model, true, up_to_level);
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
      mesh.vPos = std::move(model.positions);
      mesh.vNorm = std::move(model.normals);
      
      mesh.triIndices = std::vector<int>(model.indices.size(),0);
      memcpy(mesh.triIndices.data(),model.indices.data(),mesh.triIndices.size()*sizeof(int));
      
      mesh.vTexCoord = std::vector<GLfloat>(model.colors.size()/2, 0);
      for (int i = 0;i<model.colors.size()/4;i++)
      {
        mesh.vTexCoord[2*i] = model.colors[4*i];
        mesh.vTexCoord[2*i + 1] = model.colors[4*i + 1];
      }

      mesh.matIndices = std::vector<int>(model.indices.size()/3,wood_mat_id);

      for (int i=count/3;i<model.indices.size()/3;i++)
      {
        mesh.matIndices[i] = leaves_mat_id;
      }
    }
}
bool HydraSceneExporter::export_internal2(std::string directory, Scene &scene, Block &export_settings)
{
  const int DEMO_WIDTH  = 1024;
  const int DEMO_HEIGHT = 1024;
  
  hrErrorCallerPlace(L"demo_02_load_obj");
  hrSceneLibraryOpen(L"demos/demo_test", HR_WRITE_DISCARD);

  HRTextureNodeRef texWood = hrTexture2DCreateFromFile(L"data/textures/wood.jpg");
  HRTextureNodeRef texLeaf = hrTexture2DCreateFromFile(L"data/textures/leaf.png");
  HRTextureNodeRef texLeafOpacity = hrTexture2DCreateFromFile(L"data/textures/leaf_opacity.png");
  HRTextureNodeRef texTerrain = hrTexture2DCreateFromFile(L"data/textures/terrain2.jpg");

  HRTextureNodeRef cube[6] = {
      hrTexture2DCreateFromFile(L"data/textures/Meadow/posx.jpg"),
      hrTexture2DCreateFromFile(L"data/textures/Meadow/negx.jpg"),
      hrTexture2DCreateFromFile(L"data/textures/Meadow/posy.jpg"),
      hrTexture2DCreateFromFile(L"data/textures/Meadow/negy.jpg"),
      hrTexture2DCreateFromFile(L"data/textures/Meadow/posz.jpg"),
      hrTexture2DCreateFromFile(L"data/textures/Meadow/negz.jpg"),
    };

  HRTextureNodeRef texEnv = HRUtils::Cube2SphereLDR(cube);


  HRMaterialRef mat_land = hrMaterialCreate(L"land_material");
  hrMaterialOpen(mat_land, HR_WRITE_DISCARD);
  {
    auto matNode = hrMaterialParamNode(mat_land);
    auto diff = matNode.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type").set_value(L"lambert");

    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    color.append_attribute(L"tex_apply_mode").set_value(L"replace");
    auto texNode = hrTextureBind(texTerrain, color);

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
  /*auto diff = matNode.append_child(L"diffuse");
    diff.append_attribute(L"brdf_type").set_value(L"lambert");

    auto color = diff.append_child(L"color");
    color.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    color.append_attribute(L"tex_apply_mode").set_value(L"replace");
    auto texNode = hrTextureBind(texLeaf, color);*/

    auto opacity = matNode.append_child(L"opacity");
    auto texNodeOp = hrTextureBind(texLeafOpacity, opacity);
    
    auto transl = matNode.append_child(L"translucency");
    auto colorTrans = transl.append_child(L"color");
    colorTrans.append_attribute(L"val").set_value(L"1.0 1.0 1.0");
    colorTrans.append_attribute(L"tex_apply_mode").set_value(L"multiply");
    auto texNodeTrans = hrTextureBind(texLeaf, transl);

    VERIFY_XML(matNode);
  }
  hrMaterialClose(mat_leaf);

  SimpleMesh terrain;
  heightmap_to_simple_mesh(*(scene.heightmap), terrain);
  HRMeshRef terrainMeshRef = hrMeshCreate(L"terrain");
  hrMeshOpen(terrainMeshRef, HR_TRIANGLE_IND3, HR_WRITE_DISCARD);
  {
    hrMeshVertexAttribPointer3f(terrainMeshRef, L"pos",      &terrain.vPos[0]);
    hrMeshVertexAttribPointer3f(terrainMeshRef, L"norm",     &terrain.vNorm[0]);
    hrMeshVertexAttribPointer2f(terrainMeshRef, L"texcoord", &terrain.vTexCoord[0]);
    
    hrMeshMaterialId(terrainMeshRef, mat_land.id);
    hrMeshAppendTriangles3(terrainMeshRef, int(terrain.triIndices.size()), &terrain.triIndices[0]);
  }
  hrMeshClose(terrainMeshRef);

  std::vector<HRMeshRef> branches;
  std::vector<SimpleMesh> meshes;
  std::vector<InstanceDataArrays *> IDAs;

  for (auto &pb : scene.grove.instancedBranches)
  {
    meshes.emplace_back();
    SimpleMesh &br = meshes.back();
    packed_branch_to_simple_mesh(br,&(scene.grove), pb, 1000, true, mat_wood.id, mat_leaf.id);
    std::wstring name = L"branch" + std::to_wstring(branches.size());
    branches.push_back(hrMeshCreate(name.c_str()));
    IDAs.push_back(&(pb.IDA));
    hrMeshOpen(branches.back(), HR_TRIANGLE_IND3, HR_WRITE_DISCARD);
    logerr("open %d %d %d %d %d",meshes.size(),br.vPos.size(),br.vNorm.size(), 
           br.vTexCoord.size(), br.triIndices.size());
    {
      hrMeshVertexAttribPointer3f(branches.back(), L"pos",      &br.vPos[0]);
      hrMeshVertexAttribPointer3f(branches.back(), L"norm",     &br.vNorm[0]);
      hrMeshVertexAttribPointer2f(branches.back(), L"texcoord", &br.vTexCoord[0]);
      hrMeshPrimitiveAttribPointer1i(branches.back(), L"mind", br.matIndices.data());

      hrMeshAppendTriangles3(branches.back(), int(br.triIndices.size() - 3), br.triIndices.data());
    }
    hrMeshClose(branches.back());
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
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 8.0f;

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

    intensityNode.append_child(L"color").append_attribute(L"val") = L"1 1 1";
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 2.0f * PI;
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
    intensityNode.append_child(L"multiplier").append_attribute(L"val").set_value(L"0.75");

	  auto texNode = hrTextureBind(texEnv, intensityNode.child(L"color"));

	  VERIFY_XML(lightNode);
  }
  hrLightClose(sky);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // camera
  //
  HRCameraRef camRef = hrCameraCreate(L"my camera");
  
  hrCameraOpen(camRef, HR_WRITE_DISCARD);
  {
    xml_node camNode = hrCameraParamNode(camRef);
    
    camNode.append_child(L"fov").text().set(L"45");
    camNode.append_child(L"nearClipPlane").text().set(L"0.01");
    camNode.append_child(L"farClipPlane").text().set(L"1000.0");
    
    camNode.append_child(L"up").text().set(L"0 1 0");
    camNode.append_child(L"position").text().set(L"20 75 -200");
    camNode.append_child(L"look_at").text().set(L"0 0 0");
  
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
    node.append_child(L"maxRaysPerPixel").text()  = 512;
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
    hrMeshInstance(scnRef, terrainMeshRef, mind.L());
    
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

    hrLightInstance(scnRef, sky, mind.L());
    auto mres = hlm::mul(hlm::rotate_Z_4x4(1.0f * DEG_TO_RAD), hlm::translate4x4({0.0f, 100.0f, 0.0f}));
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
  
  hrRenderSaveFrameBufferLDR(renderRef, L"demos/demo_test/z_out.png");
  return true;
}

