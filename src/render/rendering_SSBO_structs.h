#pragma once
#include "common_utils/LiteMath_ext.h"

struct DrawElementsIndirectCommand
{
    uint count;         // Num elements (vertices)
    uint instanceCount; // Number of instances to draw (a.k.a primcount)
    uint firstIndex;    // Specifies a byte offset (cast to a pointer type) into the buffer bound to GL_ELEMENT_ARRAY_BUFFER to start reading indices from.
    uint baseVertex;    // Specifies a constant that should be added to each element of indices​ when chosing elements from the enabled vertex arrays.
    uint baseInstance;  // Specifies the base instance for use in fetching instanced vertex attributes.
    uint pad1,pad2,pad3;
};
struct LodData
{
    float2 min_max;
    float2 offset;
};
struct InstanceData
{
    float4 center_self;//xyz is center_self, w is cell_id
    float4 center_par;//xyz is center_par, w is texCoord.z in atlas
    float4x4 projection_camera;
};
struct ModelData
{
    uint LOD;
    uint type;
    uint vertexes;
    uint first_index;
    uint2 interval;
    uint culling;
    uint pad;

    float4 x_s;
    float4 y_s;
    float4 z_s;
};
struct TypeData
{
    uint offset;
    uint pad1,pad2,pad3;
};
struct currentInstancesData
{
    uint index;
    uint pad;
    float mn;
    float mx;
};
typedef uint4 currentModelsData;
typedef uint4 currentTypesData;
struct ImpostorData
{
    int slice_offset;
    int slice_verts;
    int slice_count;
    int pad1;
    float4 imp_center;
};
struct CellInfo
{
    uint lod_from;
    uint lod_to;
    uint pad1;
    uint pad2;
};