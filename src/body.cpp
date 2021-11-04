#include "body.h"
#include "tinyEngine/utility.h"
Body::Body(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c)
{
    pos = _pos;
    a = _a;
    b = _b;
    c = _c;
}
Box::Box(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c):
Body(_pos,_a,_b,_c)
{
    glm::vec3 vmin = pos;
    glm::vec3 vmax = pos;
    for (int i = 0; i <= 1; i++)
    {
        for (int j = 0; j <= 1; j++)
        {
            for (int k = 0; k <= 1; k++)
            {
                vmin = glm::min(vmin,pos + (float)i*a + (float)j*b + (float)k*c);
                vmax = glm::max(vmax,pos + (float)i*a + (float)j*b + (float)k*c);
            }
        }
    }
    bbox.a = glm::vec3(1,0,0);
    bbox.b = glm::vec3(0,1,0);
    bbox.c = glm::vec3(0,0,1);
    bbox.position = vmin;
    bbox.sizes = vmax - vmin;
    transform = glm::inverse(glm::mat4(glm::vec4(a,0),glm::vec4(b,0),glm::vec4(c,0),glm::vec4(pos,1)));
            
}
bool Box::in_body(glm::vec3 vec)
{
    glm::vec3 rp = vec - bbox.position;
    if (rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > bbox.sizes.x || rp.y > bbox.sizes.y || rp.z > bbox.sizes.z)
        return false;
    rp = glm::vec3(transform*glm::vec4(vec,1));
    return !(rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > 1 || rp.y > 1 || rp.z > 1);
}
Ellipsoid::Ellipsoid(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c):
Body(_pos,_a,_b,_c)
{
    glm::vec3 dir = a + b + c;
    bbox.a = glm::vec3(1,0,0);
    bbox.b = glm::vec3(0,1,0);
    bbox.c = glm::vec3(0,0,1);
    bbox.position = pos - dir;
    bbox.sizes = 2.0f*dir;
    transform = glm::inverse(glm::mat4(glm::vec4(a,0),glm::vec4(b,0),glm::vec4(c,0),glm::vec4(pos,1)));
    
}
bool Ellipsoid::in_body(glm::vec3 vec)
{
    glm::vec3 rp = vec - bbox.position;
    if (rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > bbox.sizes.x || rp.y > bbox.sizes.y || rp.z > bbox.sizes.z)
        return false;
    rp = glm::vec3(transform*glm::vec4(vec,1));
    return length(rp)<=1;
}
Cylinder::Cylinder(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c):
Body(_pos,_a,_b,_c)
{
    glm::vec3 dir = a + b + c;
    bbox.a = glm::vec3(1,0,0);
    bbox.b = glm::vec3(0,1,0);
    bbox.c = glm::vec3(0,0,1);
    bbox.position = pos - dir;
    bbox.sizes = 2.0f*dir;
    transform = glm::inverse(glm::mat4(glm::vec4(a,0),glm::vec4(b,0),glm::vec4(c,0),glm::vec4(pos,1)));
}
bool Cylinder::in_body(glm::vec3 vec)
{
    glm::vec3 rp = vec - bbox.position;
    if (rp.x < 0 || rp.y < 0 || rp.z < 0 || rp.x > bbox.sizes.x || rp.y > bbox.sizes.y || rp.z > bbox.sizes.z)
        return false;
    rp = glm::vec3(transform*glm::vec4(vec,1));
    return ((rp.x*rp.x + rp.y*rp.y)<=1) && (rp.z>=-1) && (rp.z<=1);
}
