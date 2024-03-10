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

  static inline glm::vec3 operator*(quat const& q, glm::vec3 const& v)
	{
		glm::vec3 const QuatVector(q.x, q.y, q.z);
		glm::vec3 const uv(glm::cross(QuatVector, v));
		glm::vec3 const uuv(glm::cross(QuatVector, uv));

		return v + ((uv * q.w) + uuv) * 2.0f;
	}


	static inline glm::vec3 operator*(glm::vec3 const& v, quat const& q)
	{
		return inverse(q) * v;
	}

  static inline glm::vec3 rotate(quat const& q, glm::vec3 const& v)
	{
		return q * v;
	}

  static inline quat angleAxis(float angle, glm::vec3 const& v)
	{
		float a = angle;
		float s = sin(a * 0.5f);
    glm::vec3 xyz = v * s;

		return quat(xyz.x, xyz.y, xyz.z, glm::cos(a * 0.5f));
	}

  static inline glm::mat3x3 to_mat3(quat const& q)
	{
		glm::mat3x3 Result(1.0f);
		float qxx(q.x * q.x);
		float qyy(q.y * q.y);
		float qzz(q.z * q.z);
		float qxz(q.x * q.z);
		float qxy(q.x * q.y);
		float qyz(q.y * q.z);
		float qwx(q.w * q.x);
		float qwy(q.w * q.y);
		float qwz(q.w * q.z);

		Result[0][0] = float(1) - float(2) * (qyy +  qzz);
		Result[0][1] = float(2) * (qxy + qwz);
		Result[0][2] = float(2) * (qxz - qwy);

		Result[1][0] = float(2) * (qxy - qwz);
		Result[1][1] = float(1) - float(2) * (qxx +  qzz);
		Result[1][2] = float(2) * (qyz + qwx);

		Result[2][0] = float(2) * (qxz + qwy);
		Result[2][1] = float(2) * (qyz - qwx);
		Result[2][2] = float(1) - float(2) * (qxx +  qyy);
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