#pragma once
#include "common_utils/LiteMath_ext.h"
#include "utility.h"
struct BBox
{
    float3 position;
    float3 sizes;
    float3 a;
    float3 b;
    float3 c;
    float V() const { return sizes.x * sizes.y * sizes.z; }
    bool special = false;
};

struct Sphere;
struct AABB
{
    float3 max_pos = float3(0, 0, 0);
    float3 min_pos = float3(0, 0, 0);

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
    AABB(float3 _min_pos, float3 _max_pos)
    {
        min_pos = _min_pos;
        max_pos = _max_pos;
    }
    inline float3 size() const
    {
      return max_pos - min_pos;
    }
    inline float3 center() const 
    {
      return 0.5f*(min_pos + max_pos);
    }
    inline float volume() const
    {
      return (max_pos.x-min_pos.x)*(max_pos.y-min_pos.y)*(max_pos.z-min_pos.z);
    }
    inline AABB expand(float ratio) const
    {
      return AABB(center()-0.5f*ratio*size(), center()+0.5f*ratio*size());
    }

    inline bool contains(const float3 &p) const
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
    inline bool intersects(const float3 &origin, const float3 &dir, float *t_near = nullptr, float *t_far = nullptr) const
    {
        float3 safe_dir = LiteMath::sign(dir)*LiteMath::max(float3(1e-9f),LiteMath::abs(dir));
        float3 tMin = (min_pos - origin)/safe_dir;
        float3 tMax = (max_pos - origin)/safe_dir;
        float3 t1 = LiteMath::min(tMin, tMax);
        float3 t2 = LiteMath::max(tMin, tMax);
        float tNear  = std::max(t1.x, std::max(t1.y, t1.z));
        float tFar   = std::min(t2.x, std::min(t2.y, t2.z));

        if (t_near)
          *t_near = tNear;
        if (t_far)
          *t_far = tFar;

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
    float3 pos;
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
    Sphere(float3 _pos, float _r)
    {
        pos = _pos;
        r = _r;
    }
    Sphere(float4 _pos_r)
    {
        pos = float3(_pos_r.x, _pos_r.y, _pos_r.z);
        r = _pos_r.w;
    }
    inline float4 compact() const { return float4(pos.x, pos.y, pos.z, r); }
    inline bool empty() const { return r <= 0; }
    inline bool contains(float3 _pos) const
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
        return AABB(pos - float3(r, r, r), pos + float3(r, r, r));
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
    float2 max_pos = float2(0, 0);
    float2 min_pos = float2(0, 0);

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
    AABB2D(float2 _min_pos, float2 _max_pos)
    {
        min_pos = _min_pos;
        max_pos = _max_pos;
    }
    AABB2D(float4 compact)
    {
        min_pos = float2(compact.x, compact.y);
        max_pos = float2(compact.z, compact.w);
    }
    inline bool contains(const float2 &p) const
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
    float2 pos;
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
    Sphere2D(float2 _pos, float _r)
    {
        pos = _pos;
        r = _r;
    }
    Sphere2D(float3 _pos_r)
    {
        pos = float2(_pos_r.x, _pos_r.y);
        r = _pos_r.z;
    }
    inline float3 compact() { return float3(pos.x, pos.y, r); }
    inline bool empty() { return r <= 0; }
    inline bool contains(float2 _pos)
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
        return AABB2D(pos - float2(r, r), pos + float2(r, r));
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

float3 Barycentric(float3 p, float3 a, float3 b, float3 c);