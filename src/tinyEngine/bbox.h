#pragma once
#include <glm/glm.hpp>
struct BBox
{
    glm::vec3 position;
    glm::vec3 sizes;
    glm::vec3 a;
    glm::vec3 b;
    glm::vec3 c;
    float V() { return sizes.x * sizes.y * sizes.z; }
    bool special = false;
};