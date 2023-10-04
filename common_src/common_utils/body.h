#pragma once
#include "common_utils/bbox.h"
class Body
{
public:
    Body(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c);
    virtual bool in_body(glm::vec3 pos) = 0;
    BBox get_Bbox() {return bbox;}
    glm::vec3 pos,a,b,c;
    glm::mat4 transform;
    BBox bbox;
};
class Box: public Body
{
public:
    Box(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c);
    bool in_body(glm::vec3 pos);
};
class Ellipsoid: public Body
{
public:
    Ellipsoid(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c);
    bool in_body(glm::vec3 pos);
};
class Cylinder: public Body
{
public:
    Cylinder(glm::vec3 _pos, glm::vec3 _a, glm::vec3 _b, glm::vec3 _c);
    bool in_body(glm::vec3 pos);
};