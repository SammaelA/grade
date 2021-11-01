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

struct AABB
{
    glm::vec3 max_pos;
    glm::vec3 min_pos;

    AABB &operator=(AABB &aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
        return *this;
    }
    AABB(glm::vec3 _min_pos, glm::vec3 _max_pos)
    {
        min_pos = _min_pos;
        max_pos = _max_pos;
    }
    inline bool contains(glm::vec3 &p)
    {
        return (p.x>=min_pos.x) && (p.x<max_pos.x) && 
               (p.x>=min_pos.x) && (p.y<max_pos.y) && 
               (p.y>=min_pos.y) && (p.z<max_pos.z); 
    }
    inline bool intersects(AABB &aabb)
    {
        return (aabb.min_pos.x <= max_pos.x) && 
               (aabb.max_pos.x >= min_pos.x) && 
               (aabb.min_pos.y <= max_pos.y) && 
               (aabb.max_pos.y >= min_pos.y) && 
               (aabb.min_pos.z <= max_pos.z) && 
               (aabb.max_pos.z >= min_pos.z);
    }
};