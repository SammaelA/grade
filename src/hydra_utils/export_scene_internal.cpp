
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
namespace hlm = HydraLiteMath;
using pugi::xml_node;

extern GLFWwindow* g_window;
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
                
                flat_terrain.vNorm.push_back(n.x);
                flat_terrain.vNorm.push_back(n.y);
                flat_terrain.vNorm.push_back(n.z);

                flat_terrain.vTexCoord.push_back(0.5*(i % 2));
                flat_terrain.vTexCoord.push_back(0.5*(j % 2));

                if (i != x - 1 && j != y - 1)
                {
                    flat_terrain.triIndices.push_back(ind);
                    flat_terrain.triIndices.push_back(ind + 1);
                    flat_terrain.triIndices.push_back(ind + y + 1);
                    flat_terrain.triIndices.push_back(ind);
                    flat_terrain.triIndices.push_back(ind + y + 1);
                    flat_terrain.triIndices.push_back(ind + y);
                }
                }
            }
        }

bool HydraSceneExporter::export_internal2(std::string directory, Scene &scene, Block &export_settings)
{
  const int DEMO_WIDTH  = 512;
  const int DEMO_HEIGHT = 512;
  
  hrErrorCallerPlace(L"demo_02_load_obj");
  hrSceneLibraryOpen(L"demos/demo_test", HR_WRITE_DISCARD);
  
  SimpleMesh sphere;
  heightmap_to_simple_mesh(*(scene.heightmap), sphere);
  //   = CreateSphere(2.0f, 128);
  HRMaterialRef mat0 = hrMaterialCreate(L"mysimplemat");
  hrMaterialOpen(mat0, HR_WRITE_DISCARD);
  {
    xml_node matNode = hrMaterialParamNode(mat0);
    xml_node diff    = matNode.append_child(L"diffuse");
    
    diff.append_attribute(L"brdf_type").set_value(L"lambert");
    diff.append_child(L"color").append_attribute(L"val") = L"0.5 0.5 0.5";
    
    VERIFY_XML(matNode); // you can verify XML parameters for the main renderer "HydraModern" in this way
  }
  hrMaterialClose(mat0);

  HRMeshRef cubeOpenRef = hrMeshCreate(L"my_box");
  
  hrMeshOpen(cubeOpenRef, HR_TRIANGLE_IND3, HR_WRITE_DISCARD);
  {
    hrMeshVertexAttribPointer4f(cubeOpenRef, L"pos",      &sphere.vPos[0]);
    hrMeshVertexAttribPointer4f(cubeOpenRef, L"norm",     &sphere.vNorm[0]);
    hrMeshVertexAttribPointer2f(cubeOpenRef, L"texcoord", &sphere.vTexCoord[0]);
    
    hrMeshMaterialId(cubeOpenRef, 0);
    //hrMeshPrimitiveAttribPointer1i(cubeOpenRef, L"mind", cubeMatIndices);
    hrMeshAppendTriangles3(cubeOpenRef, int(sphere.triIndices.size()), &sphere.triIndices[0]);
  }
  hrMeshClose(cubeOpenRef);

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
    intensityNode.append_child(L"multiplier").append_attribute(L"val") = 25.0f;
  
    VERIFY_XML(lightNode);
  }
  hrLightClose(rectLight);
  
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
    camNode.append_child(L"farClipPlane").text().set(L"100.0");
    
    camNode.append_child(L"up").text().set(L"0 1 0");
    camNode.append_child(L"position").text().set(L"0 0 14");
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
    
    node.append_child(L"trace_depth").text()      = 8;
    node.append_child(L"diff_trace_depth").text() = 4;
    node.append_child(L"maxRaysPerPixel").text()  = 128;
    node.append_child(L"qmc_variant").text()      = (HYDRA_QMC_DOF_FLAG | HYDRA_QMC_MTL_FLAG | HYDRA_QMC_LGT_FLAG); // enable all of them, results to '7'
  }
  hrRenderClose(renderRef);
  
  // create scene
  //
  HRSceneInstRef scnRef = hrSceneCreate(L"my scene");
  
  const float DEG_TO_RAD = float(3.14159265358979323846f) / 180.0f;
  
  hrSceneOpen(scnRef, HR_WRITE_DISCARD);
  {
    // instance bynny and cornell box
    //
    auto mscale     = hlm::scale4x4(hlm::float3(3,3,3));
    auto mtranslate = hlm::translate4x4(hlm::float3(1, -4.2, 0));
    auto mres       = hlm::mul(mtranslate,mscale);
    
    
    auto mrot = hlm::rotate_Y_4x4(180.0f*DEG_TO_RAD);
    hrMeshInstance(scnRef, cubeOpenRef, mrot.L());
    
    //// instance light (!!!)
    //
    mtranslate = hlm::translate4x4(hlm::float3(0, 3.85f, 0));
    hrLightInstance(scnRef, rectLight, mtranslate.L());
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

