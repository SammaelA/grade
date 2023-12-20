#pragma once
#include "utils.h"
namespace diff_render
{
struct DScene;
struct TriangleMesh;
struct Scene;

void PrintAndCompareGradients(const DScene& grad1, const DScene& grad2);
void d_finDiff(const Scene &mesh, const char* outFolder, const Img& origin, const Img& target, ::std::shared_ptr<IDiffRender> a_pDRImpl, const CamInfo& a_camData,
               DScene &d_mesh, float dPos = 1.0f, float dCol = 0.01f);


void d_finDiff2(const Scene &mesh, const char* outFolder, const Img& origin, const Img& target, ::std::shared_ptr<IDiffRender> a_pDRImpl, const CamInfo& a_camData,
               DScene &d_mesh, float dPos = 1.0f, float dCol = 0.01f);
}