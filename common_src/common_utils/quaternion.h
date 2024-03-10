#pragma once
#include "LiteMath_ext.h"

namespace LiteMath
{
  struct quat
  {
    inline quat() : x(0), y(0), z(0), w(0) {}
    inline quat(float a_x, float a_y, float a_z, float a_w) : x(a_x), y(a_y), z(a_z), w(a_w) {}
    
    inline float& operator[](int i)       { return M[i]; }
    inline float  operator[](int i) const { return M[i]; }

    union
    {
      struct { float x, y, z, w; };
      float M[4];
    };
  };

	static inline quat conjugate(quat const& q)
	{
		return quat(q.w, -q.x, -q.y, -q.z);
	}

	static inline float dot(quat const& q1, quat const& q2)
	{
		return q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w;
	}

	static inline quat operator/(quat const& q, float const& s)
	{
		return quat(q.w / s, q.x / s, q.y / s, q.z / s);
	}

	static inline quat inverse(quat const& q)
	{
		return conjugate(q) / dot(q, q);
	}

  static inline float3 operator*(quat const& q, float3 const& v)
	{
		float3 const QuatVector(q.x, q.y, q.z);
		float3 const uv(cross(QuatVector, v));
		float3 const uuv(cross(QuatVector, uv));

		return v + ((uv * q.w) + uuv) * 2.0f;
	}


	static inline float3 operator*(float3 const& v, quat const& q)
	{
		return inverse(q) * v;
	}

  static inline float3 rotate(quat const& q, float3 const& v)
	{
		return q * v;
	}

  static inline quat angleAxis(float angle, float3 const& v)
	{
		float a = angle;
		float s = sin(a * 0.5f);
    float3 xyz = v * s;

		return quat(xyz.x, xyz.y, xyz.z, cos(a * 0.5f));
	}

  static inline float3x3 to_mat3(quat const& q)
	{
		float3x3 Result;
		float qxx(q.x * q.x);
		float qyy(q.y * q.y);
		float qzz(q.z * q.z);
		float qxz(q.x * q.z);
		float qxy(q.x * q.y);
		float qyz(q.y * q.z);
		float qwx(q.w * q.x);
		float qwy(q.w * q.y);
		float qwz(q.w * q.z);

		Result.row[0] = float3(1.0f - 2.0f*(qyy +  qzz), 2.0f*(qxy + qwz)        , 2.0f*(qxz - qwy)        );
    Result.row[1] = float3(2.0f*(qxy - qwz)        , 1.0f - 2.0f*(qxx +  qzz), 2.0f*(qyz + qwx)        );
    Result.row[2] = float3(2.0f*(qxz + qwy)        , 2.0f*(qyz - qwx)        , 1.0f - 2.0f*(qxx +  qyy));
		return Result;
	}

  static inline quat operator*(quat const& p, quat const& q)
	{
		float w = p.w * q.w - p.x * q.x - p.y * q.y - p.z * q.z;
		float x = p.w * q.x + p.x * q.w + p.y * q.z - p.z * q.y;
		float y = p.w * q.y + p.y * q.w + p.z * q.x - p.x * q.z;
		float z = p.w * q.z + p.z * q.w + p.x * q.y - p.y * q.x;
		return quat(x,y,z,w);
	}
}