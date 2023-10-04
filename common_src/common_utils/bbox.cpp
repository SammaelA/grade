#include "bbox.h"

bool AABB::intersects(const Sphere &s) const
{
    float dmin = 0;
    for (int i = 0; i < 3; i++)
    {
        if (s.pos[i] < min_pos[i])
            dmin += SQR(s.pos[i] - min_pos[i]);
        else if (s.pos[i] > max_pos[i])
            dmin += SQR(s.pos[i] - max_pos[i]);
    }
    return (dmin <= s.r * s.r);
}

bool AABB2D::intersects(const Sphere2D &s) const
{
    float dmin = 0;
    for (int i = 0; i < 2; i++)
    {
        if (s.pos[i] < min_pos[i])
            dmin += SQR(s.pos[i] - min_pos[i]);
        else if (s.pos[i] > max_pos[i])
            dmin += SQR(s.pos[i] - max_pos[i]);
    }
    return (dmin <= s.r * s.r);
}

glm::vec3 Barycentric(glm::vec3 p, glm::vec3 a, glm::vec3 b, glm::vec3 c)
{    
    glm::vec3 coords;
    glm::vec3 v0 = b - a, v1 = c - a, v2 = p - a;    
    float d00 = glm::dot(v0, v0);    
    float d01 = glm::dot(v0, v1);    
    float d11 = glm::dot(v1, v1);    
    float d20 = glm::dot(v2, v0);    
    float d21 = glm::dot(v2, v1);    
    float denom = d00 * d11 - d01 * d01;    
    coords.y = (d11 * d20 - d01 * d21) / denom;    
    coords.z = (d00 * d21 - d01 * d20) / denom;   
    coords.x = 1.0f - coords.y - coords.z;

    return coords;
}