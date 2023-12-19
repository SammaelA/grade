#pragma once
#include <glm/glm.hpp>
#include "utility.h"
struct BBox
{
    glm::vec3 position;
    glm::vec3 sizes;
    glm::vec3 a;
    glm::vec3 b;
    glm::vec3 c;
    float V() const { return sizes.x * sizes.y * sizes.z; }
    bool special = false;
};

struct Sphere;
struct AABB
{
    glm::vec3 max_pos = glm::vec3(0, 0, 0);
    glm::vec3 min_pos = glm::vec3(0, 0, 0);

    AABB()
    {
    }
    AABB(const AABB &aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
    }
    AABB(AABB &&aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
    }
    AABB &operator=(const AABB &aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
        return *this;
    }
    AABB &operator=(AABB &&aabb)
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
    inline bool contains(glm::vec3 &p) const
    {
        return (p.x >= min_pos.x) && (p.x < max_pos.x) &&
               (p.y >= min_pos.y) && (p.y < max_pos.y) &&
               (p.z >= min_pos.z) && (p.z < max_pos.z);
    }
    inline bool intersects(const AABB &aabb) const
    {
        return (aabb.min_pos.x <= max_pos.x) &&
               (aabb.max_pos.x >= min_pos.x) &&
               (aabb.min_pos.y <= max_pos.y) &&
               (aabb.max_pos.y >= min_pos.y) &&
               (aabb.min_pos.z <= max_pos.z) &&
               (aabb.max_pos.z >= min_pos.z);
    }
    inline bool intersects(const glm::vec3 &origin, const glm::vec3 &dir)
    {
        glm::vec3 tMin = (min_pos - origin)/glm::max(glm::vec3(1e-9f),dir);
        glm::vec3 tMax = (max_pos - origin)/glm::max(glm::vec3(1e-9f),dir);
        glm::vec3 t1 = glm::min(tMin, tMax);
        glm::vec3 t2 = glm::max(tMin, tMax);
        float tNear  = std::max(t1.x, std::max(t1.y, t1.z));
        float tFar   = std::min(t2.x, std::min(t2.y, t2.z));

        return tNear <= tFar;
    }
    inline bool empty() const
    {
        return min_pos.x >= max_pos.x ||
               min_pos.y >= max_pos.y ||
               min_pos.z >= max_pos.z;
    }
    bool intersects(const Sphere &s) const;
    inline AABB intersect_bbox(const AABB &bbox) const
    {
        return AABB(max(min_pos, bbox.min_pos), min(max_pos, bbox.max_pos));
    }
};
struct Sphere
{
    glm::vec3 pos;
    float r;
    Sphere()
    {
    }
    Sphere(const Sphere &s)
    {
        pos = s.pos;
        r = s.r;
    }
    Sphere(Sphere &&s)
    {
        pos = s.pos;
        r = s.r;
    }
    Sphere &operator=(const Sphere &s)
    {
        pos = s.pos;
        r = s.r;
        return *this;
    }
    Sphere &operator=(Sphere &&s)
    {
        pos = s.pos;
        r = s.r;
        return *this;
    }
    Sphere(glm::vec3 _pos, float _r)
    {
        pos = _pos;
        r = _r;
    }
    Sphere(glm::vec4 _pos_r)
    {
        pos = glm::vec3(_pos_r.x, _pos_r.y, _pos_r.z);
        r = _pos_r.w;
    }
    inline glm::vec4 compact() const { return glm::vec4(pos.x, pos.y, pos.z, r); }
    inline bool empty() const { return r <= 0; }
    inline bool contains(glm::vec3 _pos) const
    {
        return (SQR(_pos.x - pos.x) + SQR(_pos.y - pos.y) + SQR(_pos.z - pos.z) <= r * r);
    }
    inline bool intersects(const Sphere &s) const
    {
        return (SQR(s.pos.x - pos.x) + SQR(s.pos.y - pos.y) + SQR(s.pos.z - pos.z) <= SQR(r + s.r));
    }
    inline bool intersects(const AABB &bbox) const
    {
        return bbox.intersects(*this);
    }
    inline AABB get_bbox() const
    {
        return AABB(pos - glm::vec3(r, r, r), pos + glm::vec3(r, r, r));
    }
    inline AABB intersect_bbox(const AABB &bbox) const
    {
        if (intersects(bbox))
        {
            return get_bbox().intersect_bbox(bbox);
        }
        else
        {
            return AABB();
        }
    }
};

struct Sphere2D;
struct AABB2D
{
    glm::vec2 max_pos = glm::vec2(0, 0);
    glm::vec2 min_pos = glm::vec2(0, 0);

    AABB2D()
    {
    }
    AABB2D(const AABB2D &aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
    }
    AABB2D(AABB2D &&aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
    }
    AABB2D &operator=(const AABB2D &aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
        return *this;
    }
    AABB2D &operator=(AABB2D &&aabb)
    {
        min_pos = aabb.min_pos;
        max_pos = aabb.max_pos;
        return *this;
    }
    AABB2D(glm::vec2 _min_pos, glm::vec2 _max_pos)
    {
        min_pos = _min_pos;
        max_pos = _max_pos;
    }
    AABB2D(glm::vec4 compact)
    {
        min_pos = glm::vec2(compact.x, compact.y);
        max_pos = glm::vec2(compact.z, compact.w);
    }
    inline bool contains(const glm::vec2 &p) const
    {
        return (p.x >= min_pos.x) && (p.x < max_pos.x) &&
               (p.x >= min_pos.x) && (p.y < max_pos.y);
    }
    inline bool intersects(const AABB2D &aabb) const
    {
        return (aabb.min_pos.x <= max_pos.x) &&
               (aabb.max_pos.x >= min_pos.x) &&
               (aabb.min_pos.y <= max_pos.y) &&
               (aabb.max_pos.y >= min_pos.y);
    }
    inline bool intersectsXZ(const AABB &aabb) const
    {
        return (aabb.min_pos.x <= max_pos.x) &&
               (aabb.max_pos.x >= min_pos.x) &&
               (aabb.min_pos.z <= max_pos.y) &&
               (aabb.max_pos.z >= min_pos.y);
    }
    inline bool empty() const
    {
        return min_pos.x >= max_pos.x ||
               min_pos.y >= max_pos.y;
    }
    bool intersects(const Sphere2D &s) const;
    inline AABB2D intersect_bbox(const AABB2D &bbox) const
    {
        return AABB2D(max(min_pos, bbox.min_pos), min(max_pos, bbox.max_pos));
    }
};
struct Sphere2D
{
    glm::vec2 pos;
    float r;
    Sphere2D()
    {
    }
    Sphere2D(const Sphere2D &s)
    {
        pos = s.pos;
        r = s.r;
    }
    Sphere2D(Sphere2D &&s)
    {
        pos = s.pos;
        r = s.r;
    }
    Sphere2D &operator=(const Sphere2D &s)
    {
        pos = s.pos;
        r = s.r;
        return *this;
    }
    Sphere2D &operator=(Sphere2D &&s)
    {
        pos = s.pos;
        r = s.r;
        return *this;
    }
    Sphere2D(glm::vec2 _pos, float _r)
    {
        pos = _pos;
        r = _r;
    }
    Sphere2D(glm::vec3 _pos_r)
    {
        pos = glm::vec2(_pos_r.x, _pos_r.y);
        r = _pos_r.z;
    }
    inline glm::vec3 compact() { return glm::vec3(pos.x, pos.y, r); }
    inline bool empty() { return r <= 0; }
    inline bool contains(glm::vec2 _pos)
    {
        return (SQR(_pos.x - pos.x) + SQR(_pos.y - pos.y) <= r * r);
    }
    inline bool intersects(const Sphere2D &s)
    {
        return (SQR(s.pos.x - pos.x) + SQR(s.pos.y - pos.y) <= SQR(r + s.r));
    }
    inline bool intersects(const AABB2D &bbox)
    {
        return bbox.intersects(*this);
    }
    inline AABB2D get_bbox()
    {
        return AABB2D(pos - glm::vec2(r, r), pos + glm::vec2(r, r));
    }
    inline AABB2D intersect_bbox(const AABB2D &bbox)
    {
        if (intersects(bbox))
        {
            return get_bbox().intersect_bbox(bbox);
        }
        else
        {
            return AABB2D();
        }
    }
};

glm::vec3 Barycentric(glm::vec3 p, glm::vec3 a, glm::vec3 b, glm::vec3 c);