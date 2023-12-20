#pragma once

#include "LiteMath.h"
#include "dmesh.h"
namespace diff_render
{
void scn01_TwoTrisFlat  (TriangleMesh& initial, TriangleMesh& target);
void scn02_TwoTrisSmooth(TriangleMesh& initial, TriangleMesh& target);
void scn03_Triangle3D_White(TriangleMesh& initial, TriangleMesh& target);
void scn04_Triangle3D_Colored(TriangleMesh& initial, TriangleMesh& target);

void scn05_Pyramid3D(TriangleMesh& initial, TriangleMesh& target);
void scn06_Cube3D_VColor(TriangleMesh& initial, TriangleMesh& target);
void scn08_Cube3D_Textured(TriangleMesh& initial, TriangleMesh& target);
void scn09_Sphere3D_Textured(TriangleMesh& initial, TriangleMesh& target);

void scn10_Teapot3D_Textured(TriangleMesh& initial, TriangleMesh& target);
void scn11_Teapot3D_Textured(TriangleMesh& initial, TriangleMesh& target);
}