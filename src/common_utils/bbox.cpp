#include "bbox.h"

bool AABB::intersects(Sphere &s)
{
    float dmin = 0;
    for (int i = 0; i < 3; i++)
    {
        if (s.pos[i] < min_pos[i])
            dmin += SQR(s.pos[i] - min_pos[i]);
        else if (s.pos[i] > max_pos[i])
            dmin += SQR(s.pos[i] - max_pos[i]);
    }
    if (dmin <= s.r * s.r)
        return true;
}

bool AABB2D::intersects(Sphere2D &s)
{
    float dmin = 0;
    for (int i = 0; i < 2; i++)
    {
        if (s.pos[i] < min_pos[i])
            dmin += SQR(s.pos[i] - min_pos[i]);
        else if (s.pos[i] > max_pos[i])
            dmin += SQR(s.pos[i] - max_pos[i]);
    }
    if (dmin <= s.r * s.r)
        return true;
}