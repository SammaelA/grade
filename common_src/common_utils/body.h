#pragma once
#include "common_utils/bbox.h"
class Body
{
public:
    Body(float3 _pos, float3 _a, float3 _b, float3 _c);
    virtual ~Body() {}
    virtual bool in_body(float3 pos) = 0;
    BBox get_Bbox() {return bbox;}
    float3 pos,a,b,c;
    float4x4 transform;
    BBox bbox;
};
class Box: public Body
{
public:
    Box(float3 _pos, float3 _a, float3 _b, float3 _c);
    bool in_body(float3 pos);
};
class Ellipsoid: public Body
{
public:
    Ellipsoid(float3 _pos, float3 _a, float3 _b, float3 _c);
    bool in_body(float3 pos);
};
class Cylinder: public Body
{
public:
    Cylinder(float3 _pos, float3 _a, float3 _b, float3 _c);
    bool in_body(float3 pos);
};