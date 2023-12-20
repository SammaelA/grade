#pragma once
#include "scene.h"
#include "dmesh.h"
#include "camera.h"
#include "utils.h"
namespace diff_render
{
struct DiffRenderSettings
{
  SHADING_MODEL mode = SHADING_MODEL::SILHOUETTE;
  int spp = 1;
};

struct IDiffRender
{
  virtual ~IDiffRender() {};
  virtual void init(const DiffRenderSettings &settings) = 0;
  virtual void commit(const Scene &scene) = 0;
  virtual void render(const Scene &scene, const CamInfo* cams, Img *imgames, int a_viewsNum) = 0;
  virtual void d_render(const Scene &scene, const CamInfo* cams, const Img *adjoints, int a_viewsNum, const int edge_samples_in_total,
                        DScene &d_mesh,
                        Img* debugImages = nullptr, int debugImageNum = 0) = 0;
  virtual float d_render_and_compare(const Scene &scene, const CamInfo* cams, const Img *target_images, int a_viewsNum, 
                                     const int edge_samples_in_total, DScene &d_mesh, Img* outImages = nullptr) = 0;
  SHADING_MODEL mode;
};
}