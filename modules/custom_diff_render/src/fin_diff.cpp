#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "myomp.h"

#include "LiteMath.h"
using namespace LiteMath;

#include <cassert>
#include <iomanip>

#include "dmesh.h"
#include "drender.h"

constexpr static int SAM_PER_PIXEL = 16;
namespace diff_render
{

void d_finDiff(const Scene &scene, const char* outFolder, const Img& origin, const Img& target, ::std::shared_ptr<IDiffRender> a_pDRImpl, const CamInfo& a_camData,
               DScene &d_scene, float dPos = 1.0f, float dCol = 0.01f) 
{
  int debug_mesh_id = 0;

  const TriangleMesh &mesh = scene.get_mesh(debug_mesh_id);
  Img img(origin.width(), origin.height());

  d_scene.reset(scene, a_pDRImpl->mode, {debug_mesh_id});
  auto &d_mesh = *d_scene.get_dmesh(debug_mesh_id);
  
  const float MSEOrigin = MSE(origin, target);
  const float scale = float(256*256*3);

  for(size_t i=0; i<mesh.vertex_count();i++)
  {
    Scene copy_scene;
    TriangleMesh copy;
    
    // dx
    //
    copy_scene = scene;
    copy = mesh;
    copy.vertices[i].x += dPos;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    auto diffToTarget = (MSE(img,target) - MSEOrigin)/dPos;
    d_mesh.pos(i)[0] += GradReal(diffToTarget*scale);
    
    // dy
    //
    copy_scene = scene;
    copy = mesh;
    copy.vertices[i].y += dPos;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);

    diffToTarget = (MSE(img,target) - MSEOrigin)/dPos;
    d_mesh.pos(i)[1] += GradReal(diffToTarget*scale);

    // dz 
    //
    copy_scene = scene;
    copy = mesh;
    copy.vertices[i].z += dPos;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    diffToTarget = (MSE(img,target) - MSEOrigin)/dPos;
    d_mesh.pos(i)[2] += GradReal(diffToTarget*scale);
  }
  
  size_t colrsNum = mesh.vertex_count();
  
  for(size_t i=0; i<colrsNum;i++)
  {
    Scene copy_scene;
    TriangleMesh copy;
    
    // d_red
    //
    copy_scene = scene;
    copy = mesh;
    copy.colors[i].x += dCol;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    auto diffToTarget = (MSE(img,target) - MSEOrigin)/dCol;
    d_mesh.color(i)[0] += GradReal(diffToTarget*scale);

    // d_green
    //
    copy_scene = scene;
    copy = mesh;
    copy.colors[i].y += dCol;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    diffToTarget = (MSE(img,target) - MSEOrigin)/dCol;
    d_mesh.color(i)[1] += GradReal(diffToTarget*scale);

    // d_blue
    //
    copy_scene = scene;
    copy = mesh;
    copy.colors[i].z += dCol;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    diffToTarget = (MSE(img,target) - MSEOrigin)/dCol;
    d_mesh.color(i)[2] += GradReal(diffToTarget*scale);
  }

}

void d_finDiff2(const Scene &scene, const char* outFolder, const Img& origin, const Img& target, ::std::shared_ptr<IDiffRender> a_pDRImpl, const CamInfo& a_camData,
                DScene &d_scene, float dPos = 1.0f, float dCol = 0.01f) 
{
  int debug_mesh_id = 0;

  const TriangleMesh &mesh = scene.get_mesh(debug_mesh_id);
  Img img(origin.width(), origin.height());

  d_scene.reset(scene, a_pDRImpl->mode, {debug_mesh_id});
  auto &d_mesh = *d_scene.get_dmesh(debug_mesh_id);
  
  const Img MSEOrigin = LiteImage::MSEImage(origin, target);

  for(size_t i=0; i<mesh.vertex_count();i++)
  {
    Scene copy_scene;
    TriangleMesh copy;
    
    // dx
    //
    copy_scene = scene;
    copy = mesh;
    copy.vertices[i].x += dPos;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    auto diffImageX = (LiteImage::MSEImage(img,target) - MSEOrigin)/dPos;   
    float3 summColor = SummOfPixels(diffImageX); 
    d_mesh.pos(i)[0] += GradReal(summColor.x + summColor.y + summColor.z);
    
    // dy
    //
    copy_scene = scene;
    copy = mesh;
    copy.vertices[i].y += dPos;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);

    auto diffImageY = (LiteImage::MSEImage(img,target) - MSEOrigin)/dPos;   
    summColor = SummOfPixels(diffImageY); 
    d_mesh.pos(i)[1] += GradReal(summColor.x + summColor.y + summColor.z);

    // dz 
    //
    copy_scene = scene;
    copy = mesh;
    copy.vertices[i].z += dPos;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);

    auto diffImageZ = (LiteImage::MSEImage(img,target) - MSEOrigin)/dPos;   
    Img diffImage(diffImageX.width(), diffImageX.height()); 
    for(int y=0;y<diffImageX.height();y++)
      for(int x=0;x<diffImageX.width();x++)
        diffImage[int2(x,y)] = float3(diffImageX[int2(x,y)].x, diffImageY[int2(x,y)].x, diffImageZ[int2(x,y)].x);

    if(outFolder != nullptr)
    {
      ::std::stringstream strOut;
      strOut << outFolder << "/" << "pos_xyz_" << i << ".bmp";
      auto path = strOut.str();
      LiteImage::SaveImage(path.c_str(), diffImage);
    }
    summColor = SummOfPixels(diffImageZ); 
    d_mesh.pos(i)[2] += GradReal(summColor.x + summColor.y + summColor.z);
  }
  
  size_t colrsNum = mesh.vertex_count();
  
  for(size_t i=0; i<colrsNum;i++)
  {
    Scene copy_scene;
    TriangleMesh copy;
    
    // d_red
    //
    copy_scene = scene;
    copy = mesh;
    copy.colors[i].x += dCol;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    auto diffToTargetX = (LiteImage::MSEImage(img,target) - MSEOrigin)/dCol;
    float3 summColor = SummOfPixels(diffToTargetX); 
    d_mesh.color(i)[0] += GradReal(summColor.x + summColor.y + summColor.z);

    // d_green
    //
    copy_scene = scene;
    copy = mesh;
    copy.colors[i].y += dCol;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    auto diffToTargetY = (LiteImage::MSEImage(img,target) - MSEOrigin)/dCol;
    summColor = SummOfPixels(diffToTargetY); 
    d_mesh.color(i)[1] += GradReal(summColor.x + summColor.y + summColor.z);

    // d_blue
    //
    copy_scene = scene;
    copy = mesh;
    copy.colors[i].z += dCol;
    img.clear(float3{0,0,0});

    copy_scene.set_mesh(copy, debug_mesh_id);
    a_pDRImpl->commit(copy_scene);
    a_pDRImpl->render(copy_scene, &a_camData, &img, 1);
    
    auto diffToTargetZ = (LiteImage::MSEImage(img,target) - MSEOrigin)/dCol;
    Img diffImage(diffToTargetX.width(), diffToTargetX.height()); 
    for(int y=0;y<diffToTargetX.height();y++)
      for(int x=0;x<diffToTargetX.width();x++)
        diffImage[int2(x,y)] = float3(diffToTargetX[int2(x,y)].x + diffToTargetX[int2(x,y)].y + diffToTargetX[int2(x,y)].z, 
                                      diffToTargetY[int2(x,y)].x + diffToTargetY[int2(x,y)].y + diffToTargetY[int2(x,y)].z, 
                                      diffToTargetZ[int2(x,y)].x + diffToTargetZ[int2(x,y)].y + diffToTargetZ[int2(x,y)].z)*0.3334f;

    if(outFolder != nullptr)
    {
      ::std::stringstream strOut;
      strOut << outFolder << "/" << "col_" << i << ".bmp";
      auto path = strOut.str();
      LiteImage::SaveImage(path.c_str(), diffImage); // 
    }
    summColor = SummOfPixels(diffToTargetZ); 
    d_mesh.color(i)[2] += GradReal(summColor.x + summColor.y + summColor.z);
  }

}

void PrintAndCompareGradients(const DScene& grad1, const DScene& grad2)
{
  double totalError = 0.0;
  double posError = 0.0;
  double colError = 0.0;
  double posLengthL1 = 0.0;
  double colLengthL1 = 0.0;

  const auto &dm1 = grad1.get_dmeshes();
  const auto &dm2 = grad2.get_dmeshes();
  assert(dm1.size() == dm2.size());

  for (int dmn = 0; dmn < dm1.size(); dmn++)
  {
    for(size_t i=0;i<3*dm1[dmn].vertices;i++) 
    {
      double diff = ::std::abs(double(dm1[dmn].pos_ptr[i] - dm2[dmn].pos_ptr[i]));
      posError    += diff;
      totalError  += diff;
      posLengthL1 += ::std::abs(dm2[dmn].pos_ptr[i]);
      ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << dm1[dmn].pos_ptr[i] << "\t";  
      ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << dm2[dmn].pos_ptr[i] << ::std::endl;
    }

    if (dm1[dmn].color_ptr)
    {
      ::std::cout << "--------------------------" << ::std::endl;
      for(size_t i=0;i<3*dm1[dmn].vertices;i++) 
      {
        double diff = ::std::abs(double(dm1[dmn].color_ptr[i] - dm2[dmn].color_ptr[i]));
        colError   += diff;
        totalError += diff;
        colLengthL1 += ::std::abs(dm2[dmn].color_ptr[i]);
        ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << dm1[dmn].color_ptr[i] << "\t";  
        ::std::cout << ::std::fixed << ::std::setw(8) << ::std::setprecision(4) << dm2[dmn].color_ptr[i] << ::std::endl;
      }
    }
    
    ::std::cout << "==========================" << ::std::endl;
    ::std::cout << "GradErr[L1](vpos ) = " << posError/double(dm1[dmn].vertices*3)    << "\t which is \t" << 100.0*(posError/posLengthL1) << "%" << ::std::endl;
    ::std::cout << "GradErr[L1](color) = " << colError/double(dm1[dmn].vertices*3)    << "\t which is \t" << 100.0*(colError/colLengthL1) << "%" << ::std::endl;
    ::std::cout << "GradErr[L1](total) = " << totalError/double(dm1[dmn].data.size()) << "\t which is \t" << 100.0*(totalError/(posLengthL1+colLengthL1)) << "%" << ::std::endl;
  }
}
}