#include "scenes.h"
#include <iostream>
#include "cmesh.h"

namespace diff_render
{
void scn01_TwoTrisFlat(TriangleMesh& initial, TriangleMesh& target)
{
  TriangleMesh mesh{
      // vertices
      {{50.0, 25.0, 0.0}, {200.0, 200.0, 0.0}, {15.0, 150.0, 0.0},
       {200.0, 15.0, 0.0}, {150.0, 250.0, 0.0}, {50.0, 100.0, 0.0}},
      // color
      {{0.3, 0.3, 0.3}, {0.3, 0.3, 0.3}},
      // indices
      {0, 1, 2, 
       3, 4, 5}
  };

  initial = mesh;

  TriangleMesh mesh2{
      // vertices
      {{50.0, 25.0+10.0, 0.0}, {200.0, 200.0+10.0, 0.0}, {15.0, 150.0+10.0, 0.0},
       {200.0-10.0 + 50.0, 15.0+5.0, 0.0}, {150.0+50.0+50.0, 250.0-25.0, 0.0}, {80.0, 100.0-25.0, 0.0}},
      // color
      {{0.3, 0.5, 0.3}, {0.3, 0.3, 0.5}},
      // indices
      {0, 1, 2, 
       3, 4, 5},
  };
  target = mesh2;
}

void scn02_TwoTrisSmooth(TriangleMesh& initial, TriangleMesh& target)
{
  TriangleMesh mesh{
      // vertices
      {{50.0, 25.0, 0.0}, {200.0, 200.0, 0.0}, {15.0, 150.0, 0.0},
       {200.0, 15.0, 0.0}, {150.0, 250.0, 0.0}, {50.0, 100.0, 0.0}},
       
      {{0.0, 0.0, 0.75}, {0.5, 0.1, 0.0}, {0.0, 0.75, 0.5},
       {0.3, 0.3, 0.5},  {0.3, 0.5, 0.3}, {0.3, 0.5, 0.3}}, 

      // indices
      {0, 1, 2, 
       3, 4, 5}
  };
  
  initial = mesh;
  ///////////////////////////////////////////////////////////////// 
  
  TriangleMesh mesh2{
      // vertices
      {{50.0, 25.0+10.0, 0.0}, {200.0, 200.0+10.0, 0.0}, {15.0, 150.0+10.0, 0.0},
       {200.0-10.0 + 50.0, 15.0+5.0, 0.0}, {150.0+50.0+50.0, 250.0-25.0, 0.0}, {80.0, 100.0-25.0, 0.0}},

       {{1,0,0}, {0,1,0}, {0,0,1},
        {1,1,0}, {1,1,0}, {1,1,0}}, 

      // indices
      {0, 1, 2, 
       3, 4, 5}
  };

  target = mesh2;
}

void scn03_Triangle3D_White(TriangleMesh& initial, TriangleMesh& target)
{
  TriangleMesh tridata;
  tridata.vertices = 
    {{0.0f, 0.0f, 0.0f},    
    {1.0f, 0.0f, 0.0f},
    {0.0f, 1.0f, 0.0f}};
  tridata.normals = 
   {{0.0f, 0.0f, -1.0f},    
    {0.0f, 0.0f, -1.0f},
    {0.0f, 0.0f, -1.0f}};
  tridata.indices = {0, 1, 2};

  initial = tridata;
  target  = tridata;
  
  // apply transforms
  //
  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.5f,0.0f));
  LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-10.0f);
  LiteMath::float4x4 mRotate2   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*+10.0f);
  
  auto mTransform1 = mTranslate*mRotate1;
  auto mTransform2 = mTranslate*mRotate2;

  transform(initial, mTransform1);
  transform(target, mTransform2);
  
  //::std::cout << "initial: [" << ::std::endl;
  //for(const auto& v : initial.vertices)
  //  ::std::cout << "[" << v[0] << ", " <<  v[1] << ", " << v[2] << "] "  << ::std::endl;
  //::std::cout << "]" << ::std::endl << ::std::endl;

  //::std::cout << "target: [" << ::std::endl;
  //for(const auto& v : target.vertices)
  //  ::std::cout << "[" << v[0] << ", " <<  v[1] << ", " << v[2] << "] "  << ::std::endl;
  //::std::cout << "]" << ::std::endl << ::std::endl;
}

void scn04_Triangle3D_Colored(TriangleMesh& initial, TriangleMesh& target)
{
   TriangleMesh tridata{
      // vertices
      {{0.0f, 1.0f, 0.0f},    
       {1.0f, -1.0f, 1.0f},
       {-1.0f, -1.0f, 1.0f},  
       },

      //// color
      {{1.0f, 0.0f, 0.0f}, 
       {1.0f, 1.0f, 0.0f}, 
       {0.0f, 0.0f, 1.0f},
       },

      // color
      //{{1.0f, 0.0f, 0.0f}, 
      // {1.0f, 0.0f, 0.0f}, 
      // {0.0f, 0.0f, 0.0f},
      // },

      // indices
      {0, 1, 2}
  };

  TriangleMesh tridata2{
      // vertices
      {{0.0f, 1.0f, 0.0f},    
       {1.0f, -1.0f, 1.0f},
       {-1.0f, -1.0f, 1.0f},  
       },

      // color
      {{0.1f, 1.0f, 0.1f}, 
       {1.0f, 0.1f, 1.0f}, 
       {0.25f, 0.25f, 0.25f},
       },
       

      // indices
      {0, 1, 2}
  };

  initial = tridata;
  target  = tridata2;
  
  // apply transforms
  //
  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.5f,0.0f));
  LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-40.0f);
  LiteMath::float4x4 mRotate2   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*+30.0f);
  
  auto mTransform1 = mTranslate*mRotate1;
  auto mTransform2 = mTranslate*mRotate2;

  transform(initial, mTransform1);
  transform(target, mTransform2);
  
  ::std::cout << "initial: [" << ::std::endl;
  for(const auto& v : initial.vertices)
    ::std::cout << "[" << v[0] << ", " <<  v[1] << ", " << v[2] << "] "  << ::std::endl;
  ::std::cout << "]" << ::std::endl << ::std::endl;

  ::std::cout << "target: [" << ::std::endl;
  for(const auto& v : target.vertices)
    ::std::cout << "[" << v[0] << ", " <<  v[1] << ", " << v[2] << "] "  << ::std::endl;
  ::std::cout << "]" << ::std::endl << ::std::endl;
}

void scn05_Pyramid3D(TriangleMesh& initial, TriangleMesh& target)
{
  TriangleMesh pyramid{
      // vertices
      {{0.0f, 1.0f, 0.0f},    
       {-1.0f, -1.0f, 1.0f},  
       {1.0f, -1.0f, 1.0f},
       {1.0f, -1.0f, -1.0f},
       {-1.0f, -1.0f, -1.0f},
       },

      // color
      {{1.0f, 0.0f, 0.0f}, 
       {1.0f, 1.0f, 0.0f}, 
       {0.0f, 0.0f, 1.0f},
       {0.0f, 1.0f, 0.0f}, 
       {0.0f, 1.0f, 1.0f}, 
       },

      // indices
      {0, 1, 2,
       0, 2, 3,
       0, 3, 4,
       0, 4, 1}
       
       //1, 2, 3,
       //1, 3, 4},
  };


  initial = pyramid;
  target  = pyramid;
  
  for(auto& c : initial.colors)
    c = float3(0.25f,0.25f,0.25f);

  // apply transforms
  //
  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.5f,0.0f));
  LiteMath::float4x4 mScale     = LiteMath::scale4x4(float3(0.5,0.5f,0.5f));
  LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-40.0f)*LiteMath::rotate4x4Z(LiteMath::DEG_TO_RAD*-30.0f);
  LiteMath::float4x4 mRotate2   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*+50.0f);
  
  auto mTransform1 = mScale*mTranslate*mRotate1;
  auto mTransform2 = mScale*mTranslate*mRotate2;

  transform(initial, mTransform1);
  transform(target, mTransform2);
}


void scn06_Cube3D_VColor(TriangleMesh& initial, TriangleMesh& target)
{
  TriangleMesh cube;
    // vertices
  cube.vertices = {{1.0f, 1.0f, -1.0f},    // Top
     {-1.0f, 1.0f, -1.0f},  
     {-1.0f, 1.0f, 1.0f},
     {1.0f, 1.0f, 1.0f},
      
     {1.0f, -1.0f, 1.0f},    // Bottom
     {-1.0f, -1.0f, 1.0f},
     {-1.0f, -1.0f, -1.0f},
     {1.0f, -1.0f, -1.0f},

     {1.0f, 1.0f, 1.0f},     // Front
     {-1.0f, 1.0f, 1.0f},
     {-1.0f, -1.0f, 1.0f},
     {1.0f, -1.0f, 1.0f},

     {1.0f, -1.0f, -1.0f},   // Back
     {-1.0f, -1.0f, -1.0f},
     {-1.0f, 1.0f, -1.0f},
     {1.0f, 1.0f, -1.0f},

     {-1.0f, 1.0f, 1.0f},    // Left
     {-1.0f, 1.0f, -1.0f},
     {-1.0f, -1.0f, -1.0f},
     {-1.0f, -1.0f, 1.0f},

     {1.0f, 1.0f, -1.0f},    // Right
     {1.0f, 1.0f, 1.0f}, 
     {1.0f, -1.0f, 1.0f},
     {1.0f, -1.0f, -1.0f},

    };
  cube.normals ={
     {0,1,0},    // Top
     {0,1,0},  
     {0,1,0},
     {0,1,0},
      
     {0,-1,0},    // Bottom
     {0,-1,0},    
     {0,-1,0},  
     {0,-1,0},

     {0, 0, 1.0f},     // Front
     {0, 0, 1.0f},
     {0, 0, 1.0f},
     {0, 0, 1.0f},

     {0, 0, -1.0f},   // Back
     {0, 0, -1.0f},
     {0, 0, -1.0f},
     {0, 0, -1.0f},

     {-1.0f, 0, 0},    // Left
     {-1.0f, 0, 0},
     {-1.0f, 0, 0},
     {-1.0f, 0, 0},

     {1.0f, 0, 0},    // Right
     {1.0f, 0, 0}, 
     {1.0f, 0, 0},
     {1.0f, 0, 0},

    };
    // color
    cube.colors = {{0.0f, 1.0f, 0.0f},    // Top
     {0.0f, 1.0f, 0.0f}, 
     {0.0f, 1.0f, 0.0f},
     {0.0f, 1.0f, 0.0f},

     {1.0f, 0.5f, 0.0f},    // Bottom
     {1.0f, 0.5f, 0.0f}, 
     {1.0f, 0.5f, 0.0f},
     {1.0f, 0.5f, 0.0f},

     {1.0f, 0.0f, 0.0f},    // Front
     {1.0f, 0.0f, 0.0f},
     {1.0f, 0.0f, 0.0f},
     {1.0f, 0.0f, 0.0f},

     {1.0f, 1.0f, 0.0f},    // Back
     {1.0f, 1.0f, 0.0f},
     {1.0f, 1.0f, 0.0f},
     {1.0f, 1.0f, 0.0f},

     {0.0f, 0.0f, 1.0f},    // Left
     {0.0f, 0.0f, 1.0f},
     {0.0f, 0.0f, 1.0f},
     {0.0f, 0.0f, 1.0f},

     {1.0f, 0.0f, 1.0f},    // Right
     {1.0f, 0.0f, 1.0f},
     {1.0f, 0.0f, 1.0f},
     {1.0f, 0.0f, 1.0f},

    };

  cube.indices.resize(6*2*3); // 6 faces, 2 triangles per face, 3 indices per triangle 

  for(int face=0; face < 6; face++) // GL_QUADS
  {  
    cube.indices[face*6+0] = face*4+0;
    cube.indices[face*6+1] = face*4+1;
    cube.indices[face*6+2] = face*4+2;

    cube.indices[face*6+3] = face*4+0;
    cube.indices[face*6+4] = face*4+2;
    cube.indices[face*6+5] = face*4+3;
  }

  initial = cube;
  target  = cube;
  
  for(auto& c : initial.colors)
    c = float3(0.25f, 0.25f, 0.25f);

  // apply transforms
  //
  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.0f,0.0f));
  LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-35.0f)*LiteMath::rotate4x4Z(LiteMath::DEG_TO_RAD*-40.0f);
  LiteMath::float4x4 mRotate2   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*60.0f)*LiteMath::rotate4x4Z(LiteMath::DEG_TO_RAD*-20.0f);;
  
  auto mTransform1 = mTranslate*mRotate1;
  auto mTransform2 = mTranslate*mRotate2;

  transform(initial, mTransform1);
  transform(target, mTransform2);
}

void scn08_Cube3D_Textured(TriangleMesh& initial, TriangleMesh& target)
{
   TriangleMesh cube{
    // vertices
    {{1.0f, 1.0f, -1.0f},    // 0 Top
     {-1.0f, 1.0f, -1.0f},   // 1
     {-1.0f, 1.0f, 1.0f},    // 2
     {1.0f, 1.0f, 1.0f},     // 3
      
     {1.0f, -1.0f, 1.0f},    // 4 Bottom
     {-1.0f, -1.0f, 1.0f},   // 5
     {-1.0f, -1.0f, -1.0f},  // 6
     {1.0f, -1.0f, -1.0f},   // 7
    },

    // color
    {{0.0f, 1.0f, 0.0f},    // Top
     {0.0f, 1.0f, 0.0f}, 
   
     {1.0f, 0.5f, 0.0f},    // Bottom
     {1.0f, 0.5f, 0.0f}, 

     {1.0f, 0.0f, 0.0f},    // Front
     {1.0f, 0.0f, 0.0f},

     {1.0f, 1.0f, 0.0f},    // Back
     {1.0f, 1.0f, 0.0f},

     {0.0f, 0.0f, 1.0f},    // Left
     {0.0f, 0.0f, 1.0f},

     {1.0f, 0.0f, 1.0f},    // Right
     {1.0f, 0.0f, 1.0f},
    },

  };

  cube.indices.resize(6*2*3); // 6 faces, 2 triangles per face, 3 indices per triangle 

  cube.indices[0] = 0; cube.indices[1] = 1; cube.indices[2] = 2;        cube.indices[3] = 0; cube.indices[4] = 2; cube.indices[5] = 3;
  cube.indices[6] = 4; cube.indices[7] = 5; cube.indices[8] = 6;        cube.indices[9] = 4; cube.indices[10] = 6; cube.indices[11] = 7;

  cube.indices[12] = 1; cube.indices[13] = 6; cube.indices[14] = 7;     cube.indices[15] = 1; cube.indices[16] = 7; cube.indices[17] = 0;
  cube.indices[18] = 3; cube.indices[19] = 4; cube.indices[20] = 5;     cube.indices[21] = 3; cube.indices[22] = 5; cube.indices[23] = 2; 

  cube.indices[24] = 2; cube.indices[25] = 1; cube.indices[26] = 6;     cube.indices[27] = 2; cube.indices[28] = 6; cube.indices[29] = 5;
  cube.indices[30] = 0; cube.indices[31] = 3; cube.indices[32] = 4;     cube.indices[33] = 0; cube.indices[34] = 4; cube.indices[35] = 7; 

  cube.tc = ::std::vector<float2>{
    {1.0f, 1.0f},    // 0 Top
    {0.0f, 1.0f},   // 1
    {0.0f, 1.0f},    // 2
    {1.0f, 1.0f},     // 3
      
    {1.0f, 0.0f},    // 4 Bottom
    {0.0f, 0.0f},   // 5
    {0.0f, 0.0f},  // 6
    {1.0f, 0.0f,},   // 7
  };

  {
    int w = 256;
    int h = 256;
    cube.textures.emplace_back();
    cube.textures.back().w = w;
    cube.textures.back().h = h;
    cube.textures.back().channels = 3;
    cube.textures.back().data.resize(w*h*3);
    for (int j=0;j<h;j++)
    {
      for (int i=0;i<w;i++)
      {
        float v = 0;
        if ((i/16)%2 != (j/16)%2)
          v = 1;
        cube.textures.back().data[3*(j*w+i)] = v;
        cube.textures.back().data[3*(j*w+i)+1] = 1-v;
        cube.textures.back().data[3*(j*w+i)+2] = 0;
      }
    }

  }

  initial = cube;
  target  = cube;
  
  initial.textures[0].data = ::std::vector<float>(target.textures[0].data.size(), 0.5);

  // testing texture reconstruction, so apply same transforms
  //
  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.0f,0.0f));
  LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-35.0f);
  
  auto mTransform1 = mTranslate*mRotate1;

  transform(initial, mTransform1);
  transform(target, mTransform1);
}

void CreateSphere(TriangleMesh &sphere, float radius, int numberSlices)
{
  int i, j;

  int numberParallels = numberSlices;
  int numberVertices  = (numberParallels + 1) * (numberSlices + 1);
  int numberIndices   = numberParallels * numberSlices * 3;

  float angleStep     = (2.0f * 3.14159265358979323846f) / ((float)numberSlices);
  //float helpVector[3] = { 0.0f, 1.0f, 0.0f };
 
  sphere.vertices.resize(numberVertices);
  sphere.normals.resize(numberVertices);
  sphere.tc.resize(numberVertices);
  sphere.indices.resize(numberIndices);


  for (i = 0; i < numberParallels + 1; i++)
  {
    for (j = 0; j < numberSlices + 1; j++)
    {
      int ind    = (i * (numberSlices + 1) + j);

      sphere.vertices[ind].x = radius * sinf(angleStep * (float)i) * sinf(angleStep * (float)j);
      sphere.vertices[ind].y = radius * cosf(angleStep * (float)i);
      sphere.vertices[ind].z = radius * sinf(angleStep * (float)i) * cosf(angleStep * (float)j);

      sphere.normals[ind] = sphere.vertices[ind] / radius;

      sphere.tc[ind].x = (float)j / (float)numberSlices;
      sphere.tc[ind].y = (1.0f - (float)i) / (float)(numberParallels - 1);
    }
  }

  unsigned* indexBuf = sphere.indices.data();

  for (i = 0; i < numberParallels; i++)
  {
    for (j = 0; j < numberSlices; j++)
    {
      *indexBuf++ = i * (numberSlices + 1) + j;
      *indexBuf++ = (i + 1) * (numberSlices + 1) + j;
      *indexBuf++ = (i + 1) * (numberSlices + 1) + (j + 1);

      *indexBuf++ = i * (numberSlices + 1) + j;
      *indexBuf++ = (i + 1) * (numberSlices + 1) + (j + 1);
      *indexBuf++ = i * (numberSlices + 1) + (j + 1);
      
      int diff = int(indexBuf - sphere.indices.data());
      if (diff >= numberIndices)
        break;
    }

    int diff = int(indexBuf - sphere.indices.data());
    if (diff >= numberIndices)
      break;
  }
}

CPUTexture RedGreenLinesTexture()
{
  CPUTexture tex;
  int w = 256;
  int h = 256;
  tex.w = w;
  tex.h = h;
  tex.channels = 3;
  tex.data.resize(w*h*3);
  
  for (int j=0;j<h;j++)
  {
    for (int i=0;i<w;i++)
    {
      float v = 0;
      if ((i/16)%2 != (j/16)%2)
        v = 1;
      tex.data[3*(j*w+i)] = v;
      tex.data[3*(j*w+i)+1] = 1-v;
      tex.data[3*(j*w+i)+2] = 0;
    }
  }
  return tex;
}


CPUTexture WhiteTexture()
{
  CPUTexture tex;
  int w = 256;
  int h = 256;
  tex.w = w;
  tex.h = h;
  tex.channels = 3;
  tex.data = ::std::vector<float>(w*h*3, 0.5);

  return tex;
}

void scn09_Sphere3D_Textured(TriangleMesh& initial, TriangleMesh& target)
{
  TriangleMesh sphere;
  CreateSphere(sphere, 0.5, 32);
  sphere.textures.emplace_back(RedGreenLinesTexture());

  initial = sphere;
  target  = sphere;
  
  initial.textures[0].data = ::std::vector<float>(target.textures[0].data.size(), 0.5);

  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.2f,0.0f));
  
  //transform(initial, mTranslate);
}

void scn10_Teapot3D_Textured(TriangleMesh& initial, TriangleMesh& target)
{
  //mesh_utils::SimpleMesh tmpMesh   = mesh_utils::LoadMeshFromVSGF("data/meshes/teapot_16K.vsgf");
  mesh_utils::SimpleMesh tmpMesh   = mesh_utils::LoadMeshFromVSGF("data/meshes/vase.vsgf");
  const size_t numberVertices = tmpMesh.VerticesNum();
  const size_t numberIndices  = tmpMesh.IndicesNum();

  initial.vertices.resize(numberVertices);
  initial.normals.resize(numberVertices);
  initial.tc.resize(numberVertices);
  initial.indices.resize(numberIndices);

  for(size_t i=0;i<numberVertices;i++)
  {
    initial.vertices[i] = LiteMath::to_float3(tmpMesh.vPos4f[i]);
    initial.normals [i] = LiteMath::to_float3(tmpMesh.vNorm4f[i]);
    initial.tc      [i] = tmpMesh.vTexCoord2f[i];
  }
  
  for(size_t i=0;i<numberIndices;i++)
    initial.indices[i] = tmpMesh.indices[i];

  
  initial.textures.emplace_back(RedGreenLinesTexture());

  target = initial;
  initial.textures[0].data = ::std::vector<float>(target.textures[0].data.size(), 0.5);

  // testing texture reconstruction, so apply same transforms
  //
  LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.0f,0.0f));
  LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-25.0f);
  
  auto mTransform1 = mTranslate*mRotate1;

  transform(initial, mTransform1);
  transform(target, mTransform1);
}

void scn11_Teapot3D_Textured(TriangleMesh& initial, TriangleMesh& target)
{
  mesh_utils::SimpleMesh tmpMesh   = mesh_utils::LoadMeshFromVSGF("data/meshes/teapot.vsgf");
  //mesh_utils::SimpleMesh tmpMesh   = mesh_utils::LoadMeshFromVSGF("data/meshes/vase.vsgf");
  const size_t numberVertices = tmpMesh.VerticesNum();
  const size_t numberIndices  = tmpMesh.IndicesNum();

  initial.vertices.resize(numberVertices);
  initial.normals.resize(numberVertices);
  initial.tangents.resize(numberVertices);
  initial.tc.resize(numberVertices);
  initial.indices.resize(numberIndices);

  ::std::vector<int> missed_normals;

  for(size_t i=0;i<numberVertices;i++)
  {
    initial.vertices[i] = LiteMath::to_float3(tmpMesh.vPos4f[i]);
    initial.tc      [i] = tmpMesh.vTexCoord2f[i];

    float3 n = LiteMath::to_float3(tmpMesh.vNorm4f[i]); 
    if (length(n) > 1e-4)
      initial.normals [i] = normalize(n);
    else
      missed_normals.push_back(i);
  }
  
  for(size_t i=0;i<numberIndices;i++)
    initial.indices[i] = tmpMesh.indices[i];

  ::std::vector<float3> accumulated_normals(missed_normals.size(),{0,0,0});
  int k = 0;
  for (auto &ni : missed_normals)
  {
    for (int i=0;i<initial.indices.size(); i+=3)
    {
      int A = initial.indices[i];
      int B = initial.indices[i+1];
      int C = initial.indices[i+2];
      if (ni == A || ni == B || ni == C)
        accumulated_normals[k] += normalize(cross(initial.vertices[B]-initial.vertices[A], initial.vertices[C]-initial.vertices[A]));
    }
    initial.normals[ni] = normalize(accumulated_normals[k]);
    k++;
  }
  for(size_t i=0;i<numberVertices;i++)
  {
    float3 tang = LiteMath::to_float3(tmpMesh.vTang4f[i]);
    if (length(tang) > 1e-4)
      initial.tangents[i] = normalize(tang);
    else if (fabsf(initial.normals[i].y) < 0.99999)
      initial.tangents[i] = normalize(cross(initial.normals[i], float3(0,1,0)));
    else
      initial.tangents[i] = normalize(cross(initial.normals[i], float3(1,0,0)));
  }
  
  initial.textures.emplace_back(RedGreenLinesTexture());

  target = initial;
  initial.textures[0].data = ::std::vector<float>(target.textures[0].data.size(), 0.5);

  // testing texture reconstruction, so apply same transforms
  //
  //LiteMath::float4x4 mTranslate = LiteMath::translate4x4(float3(0,+0.0f,0.0f));
  //LiteMath::float4x4 mRotate1   = LiteMath::rotate4x4Y(LiteMath::DEG_TO_RAD*-25.0f);
  
  //auto mTransform1 = mTranslate*mRotate1;

  //for(auto& v : initial.vertices)
  //  v = mTransform1*v;

  //for(auto& v : target.vertices)
  //  v = mTransform1*v;

}
}