#include "common_utils/body.h"
#include "common_utils/utility.h"
Body::Body(float3 _pos, float3 _a, float3 _b, float3 _c)
{
    pos = _pos;
    a = _a;
    b = _b;
    c = _c;
}
Box::Box(float3 _pos, float3 _a, float3 _b, float3 _c):
Body(_pos,_a,_b,_c)
{
    float3 vmin = pos;
    float3 vmax = pos;
    for (int i = 0; i <= 1; i++)
    {
        for (int j = 0; j <= 1; j++)
        {
            for (int k = 0; k <= 1; k++)
            {
                vmin = LiteMath::min(vmin,pos + (float)i*a + (float)j*b + (float)k*c);
                vmax = LiteMath::max(vmax,pos + (float)i*a + (float)j*b + (float)k*c);
            }
        }
    }
    bbox.a = float3(1,0,0);
    bbox.b = float3(0,1,0);
    bbox.c = float3(0,0,1);
    bbox.position = vmin;
    bbox.sizes = vmax - vmin;
    transform = LiteMath::inverse4x4(to_float4x4(to_float4(a,0),to_float4(b,0),to_float4(c,0),to_float4(pos,1)));
            
}
bool Box::in_body(float3 vec)
{
    float3 rp = vec - bbox.position;
    if (rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > bbox.sizes.x || rp.y > bbox.sizes.y || rp.z > bbox.sizes.z)
        return false;
    rp = to_float3(transform*to_float4(vec,1));
    return !(rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > 1 || rp.y > 1 || rp.z > 1);
}
Ellipsoid::Ellipsoid(float3 _pos, float3 _a, float3 _b, float3 _c):
Body(_pos,_a,_b,_c)
{
    float3 dir = a + b + c;
    bbox.a = float3(1,0,0);
    bbox.b = float3(0,1,0);
    bbox.c = float3(0,0,1);
    bbox.position = pos - dir;
    bbox.sizes = 2.0f*dir;
    transform = LiteMath::inverse4x4(to_float4x4(to_float4(a,0),to_float4(b,0),to_float4(c,0),to_float4(pos,1)));
    
}
bool Ellipsoid::in_body(float3 vec)
{
    float3 rp = vec - bbox.position;
    if (rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > bbox.sizes.x || rp.y > bbox.sizes.y || rp.z > bbox.sizes.z)
        return false;
    rp = to_float3(transform*to_float4(vec,1));
    return length(rp)<=1;
}
Cylinder::Cylinder(float3 _pos, float3 _a, float3 _b, float3 _c):
Body(_pos,_a,_b,_c)
{
    float3 dir = a + b + c;
    bbox.a = float3(1,0,0);
    bbox.b = float3(0,1,0);
    bbox.c = float3(0,0,1);
    bbox.position = pos - dir;
    bbox.sizes = 2.0f*dir;
    transform = LiteMath::inverse4x4(to_float4x4(to_float4(a,0),to_float4(b,0),to_float4(c,0),to_float4(pos,1)));
}
bool Cylinder::in_body(float3 vec)
{
    float3 rp = vec - bbox.position;
    if (rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > bbox.sizes.x || rp.y > bbox.sizes.y || rp.z > bbox.sizes.z)
        return false;
    rp = to_float3(transform*to_float4(vec,1));
    return ((rp.x*rp.x + rp.y*rp.y)<=1) && (rp.z>=-1) && (rp.z<=1);
}
