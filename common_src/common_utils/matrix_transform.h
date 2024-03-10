#pragma once
#define GLM_ENABLE_EXPERIMENTAL 1
#include "LiteMath_ext.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/euler_angles.hpp>

namespace LiteMath
{
  static inline glm::mat4 eulerAngleXYZ(float x_angle, float y_angle, float z_angle)
  {
    float c1 = cos(-x_angle);
    float c2 = cos(-y_angle);
    float c3 = cos(-z_angle);
    float s1 = sin(-x_angle);
    float s2 = sin(-y_angle);
    float s3 = sin(-z_angle);

    glm::mat4 v1;
    v1[0][0] = c2 * c3;
    v1[0][1] = -c1 * s3 + s1 * s2 * c3;
    v1[0][2] = s1 * s3 + c1 * s2 * c3;
    v1[0][3] = 0;
    v1[1][0] = c2 * s3;
    v1[1][1] = c1 * c3 + s1 * s2 * s3;
    v1[1][2] = -s1 * c3 + c1 * s2 * s3;
    v1[1][3] = 0;
    v1[2][0] = -s2;
    v1[2][1] = s1 * c2;
    v1[2][2] = c1 * c2;
    v1[2][3] = 0;
    v1[3][0] = 0;
    v1[3][1] = 0;
    v1[3][2] = 0;
    v1[3][3] = 1;

    return v1;
  }

  static inline glm::mat4 lookAtRH(const glm::vec3 &eye, const glm::vec3 &center, const glm::vec3 &up)
	{
		const glm::vec3  f(normalize(center - eye));
		const glm::vec3  s(normalize(cross(f, up)));
		const glm::vec3  u(cross(s, f));

		glm::mat4 Result(1);
		Result[0][0] = s.x;
		Result[1][0] = s.y;
		Result[2][0] = s.z;
		Result[0][1] = u.x;
		Result[1][1] = u.y;
		Result[2][1] = u.z;
		Result[0][2] =-f.x;
		Result[1][2] =-f.y;
		Result[2][2] =-f.z;
		Result[3][0] =-dot(s, eye);
		Result[3][1] =-dot(u, eye);
		Result[3][2] = dot(f, eye);
		return Result;
	}

	static inline glm::mat4 lookAtLH(const glm::vec3 & eye, const glm::vec3 & center, const glm::vec3 & up)
	{
		const glm::vec3  f(normalize(center - eye));
		const glm::vec3  s(normalize(cross(up, f)));
		const glm::vec3  u(cross(f, s));

		glm::mat4 Result(1);
		Result[0][0] = s.x;
		Result[1][0] = s.y;
		Result[2][0] = s.z;
		Result[0][1] = u.x;
		Result[1][1] = u.y;
		Result[2][1] = u.z;
		Result[0][2] = f.x;
		Result[1][2] = f.y;
		Result[2][2] = f.z;
		Result[3][0] = -dot(s, eye);
		Result[3][1] = -dot(u, eye);
		Result[3][2] = -dot(f, eye);
		return Result;
	}

  //LookAt matrix is right-handed by default
  static inline glm::mat4 lookAt(const glm::vec3 & eye, const glm::vec3 & center, const glm::vec3 & up)
  {
    return lookAtRH(eye, center, up);
  }

  //perspective projection with right-handed coordinate system to unit cube [-1,1]^3
  //fov in radians
  static inline glm::mat4 perspectiveRH_NO(float fovy, float aspect, float zNear, float zFar)
	{
    assert(aspect > 0.001f);
		assert(fovy > 0.001f);
    assert(zNear > 0);
    assert(zFar > zNear);

		float tanHalfFovy = tan(fovy / 2.0f);

		glm::mat4 Result(0.0f);
		Result[0][0] = 1.0f / (aspect * tanHalfFovy);
		Result[1][1] = 1.0f / (tanHalfFovy);
		Result[2][2] = - (zFar + zNear) / (zFar - zNear);
		Result[2][3] = - 1.0f;
		Result[3][2] = - (2.0f * zFar * zNear) / (zFar - zNear);
		return Result;
	}

  //perspective projection
  static inline glm::mat4 perspective(float fovy, float aspect, float zNear, float zFar)
  {
    return perspectiveRH_NO(fovy, aspect, zNear, zFar);
  }

  //orthographic projection with right-handed coordinate system to unit cube [-1,1]^3
  static inline glm::mat4 orthoRH_NO(float left, float right, float bottom, float top, float zNear, float zFar)
	{
		glm::mat4 Result(1.0f);
		Result[0][0] =  2.0f / (right - left);
		Result[1][1] =  2.0f / (top - bottom);
		Result[2][2] = -2.0f / (zFar - zNear);
		Result[3][0] = - (right + left) / (right - left);
		Result[3][1] = - (top + bottom) / (top - bottom);
		Result[3][2] = - (zFar + zNear) / (zFar - zNear);
		return Result;
	}

  //orthographic projection
  static inline glm::mat4 ortho(float left, float right, float bottom, float top, float zNear, float zFar)
	{
    return orthoRH_NO(left, right, bottom, top, zNear, zFar);
  }

  static inline glm::mat4 rotate(const glm::mat4 &m, float angle, const glm::vec3 &v)
	{
		float a = angle;
		float c = cos(a);
		float s = sin(a);

		glm::vec3 axis(normalize(v));
		glm::vec3 temp((1.0f - c) * axis);

		glm::mat4 Rotate;
		Rotate[0][0] = c + temp[0] * axis[0];
		Rotate[0][1] = temp[0] * axis[1] + s * axis[2];
		Rotate[0][2] = temp[0] * axis[2] - s * axis[1];

		Rotate[1][0] = temp[1] * axis[0] - s * axis[2];
		Rotate[1][1] = c + temp[1] * axis[1];
		Rotate[1][2] = temp[1] * axis[2] + s * axis[0];

		Rotate[2][0] = temp[2] * axis[0] + s * axis[1];
		Rotate[2][1] = temp[2] * axis[1] - s * axis[0];
		Rotate[2][2] = c + temp[2] * axis[2];

		glm::mat4 Result;
		Result[0] = m[0] * Rotate[0][0] + m[1] * Rotate[0][1] + m[2] * Rotate[0][2];
		Result[1] = m[0] * Rotate[1][0] + m[1] * Rotate[1][1] + m[2] * Rotate[1][2];
		Result[2] = m[0] * Rotate[2][0] + m[1] * Rotate[2][1] + m[2] * Rotate[2][2];
		Result[3] = m[3];
		
    return Result;
	}

  static inline glm::mat4 translate(glm::mat4 const& m, glm::vec3 const& v)
	{
		glm::mat4 Result(m);
		Result[3] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2] + m[3];
    return Result;
	}

  static inline glm::mat4 scale(glm::mat4 const& m, glm::vec3 const& v)
	{
		glm::mat4 Result;
		Result[0] = m[0] * v[0];
		Result[1] = m[1] * v[1];
		Result[2] = m[2] * v[2];
		Result[3] = m[3];
    return Result;
	}
}