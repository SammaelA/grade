#pragma once
#include <vector>
#include <array>
#include <cppad/cppad.hpp>

namespace CppAD
{
  template <class Base>
  class AD;
}
namespace dgen
{
  template <typename T>
  struct g_vec2
  {
    union
    {
      T data[2];
      struct
      {
        T x;
        T y;
      };
    };
    g_vec2<T>(T _x, T _y)
    {
      x = _x;
      y = _y;
    }
    g_vec2<T>(const g_vec2<T> &v) = default;
    g_vec2<T>(): g_vec2<T>(0,0)
    {

    }
    ~g_vec2<T>()
    {
      (&x)->T::~T();
      (&y)->T::~T();
    }
    g_vec2<T> &operator=(const g_vec2<T> &v) = default;
    g_vec2<T> &operator=(g_vec2<T> &&v) = default;

    T &operator[](int index)
    {
      return data[index];
    }
    const T &operator[](int index) const
    {
      return data[index];
    }
    g_vec2<T> operator+(const g_vec2<T> &rhs) const
    {
      g_vec2<T> res = *this;
      res.x += rhs.x;
      res.y += rhs.y;
      return res;
    }
    g_vec2<T> operator-(const g_vec2<T> &rhs) const
    {
      g_vec2<T> res = *this;
      res.x -= rhs.x;
      res.y -= rhs.y;
      return res;
    }
    g_vec2<T> operator*(const g_vec2<T> &rhs) const
    {
      g_vec2<T> res = *this;
      res.x *= rhs.x;
      res.y *= rhs.y;
      return res;
    }
    g_vec2<T> operator*(T a) const
    {
      g_vec2<T> res = *this;
      res.x *= a;
      res.y *= a;
      return res;
    }
    g_vec2<T> operator/(const g_vec2<T> &rhs) const
    {
      g_vec2<T> res = *this;
      res.x /= rhs.x;
      res.y /= rhs.y;
      return res;
    }
    g_vec2<T> operator/(T a) const
    {
      g_vec2<T> res = *this;
      res.x /= a;
      res.y /= a;
      return res;
    }
    g_vec2<T> &operator+=(const g_vec2<T> &rhs)
    {
      x += rhs.x;
      y += rhs.y;
      return *this;
    }
    g_vec2<T> &operator-=(const g_vec2<T> &rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      return *this;
    }
    g_vec2<T> &operator*=(const g_vec2<T> &rhs)
    {
      x *= rhs.x;
      y *= rhs.y;
      return *this;
    }
    g_vec2<T> &operator*=(const T a)
    {
      x *= a;
      y *= a;
      return *this;
    }
    g_vec2<T> &operator/=(const g_vec2<T> &rhs)
    {
      x /= rhs.x;
      y /= rhs.y;
      return *this;
    }
    g_vec2<T> &operator/=(const T a)
    {
      x /= a;
      y /= a;
      return *this;
    }
    inline bool operator==(const g_vec2<T> &rhs) const
    {
      return (x == rhs.x) && (y == rhs.y);
    }
    inline bool operator!=(const g_vec2<T> &rhs) const
    {
      return !(*this == rhs);
    }
    g_vec2<T> operator-()
    {
      return g_vec2<T>(-x, -y);
    }
  };

  template <typename T>
  g_vec2<T> operator*(T a, const g_vec2<T> &rhs)
  {
    g_vec2<T> res = rhs;
    res.x *= a;
    res.y *= a;
    return res;
  }

  template <typename T>
  T length(const g_vec2<T> &v)
  {
    return CppAD::sqrt(v.x * v.x + v.y * v.y);
  }

  template <typename T>
  T dot(const g_vec2<T> &v1, const g_vec2<T> &v2)
  {
    return v1.x * v2.x + v1.y * v2.y;
  }

  template <typename T>
  g_vec2<T> normalize(const g_vec2<T> &v)
  {
    T l = length(v);
    if (l > 1e-9)
      return v / l;
    else
      return g_vec2<T>(0,0);
  }

  template <typename T>
  struct g_vec3
  {
    union
    {
      T data[3];
      struct
      {
        T x;
        T y;
        T z;
      };
    };
    g_vec3<T>(T _x, T _y, T _z)
    {
      x = _x;
      y = _y;
      z = _z;
    }
    g_vec3<T>(const g_vec3<T> &v) = default;
    g_vec3<T>(): g_vec3<T>(0,0,0)
    {

    }
    ~g_vec3<T>()
    {
      (&x)->T::~T();
      (&y)->T::~T();
      (&z)->T::~T();
    }
    g_vec3<T> &operator=(const g_vec3<T> &v) = default;
    g_vec3<T> &operator=(g_vec3<T> &&v) = default;

    T &operator[](int index)
    {
      return data[index];
    }
    const T &operator[](int index) const
    {
      return data[index];
    }
    g_vec3<T> operator+(const g_vec3<T> &rhs) const
    {
      g_vec3<T> res = *this;
      res.x += rhs.x;
      res.y += rhs.y;
      res.z += rhs.z;
      return res;
    }
    g_vec3<T> operator-(const g_vec3<T> &rhs) const
    {
      g_vec3<T> res = *this;
      res.x -= rhs.x;
      res.y -= rhs.y;
      res.z -= rhs.z;
      return res;
    }
    g_vec3<T> operator*(const g_vec3<T> &rhs) const
    {
      g_vec3<T> res = *this;
      res.x *= rhs.x;
      res.y *= rhs.y;
      res.z *= rhs.z;
      return res;
    }
    g_vec3<T> operator*(T a) const
    {
      g_vec3<T> res = *this;
      res.x *= a;
      res.y *= a;
      res.z *= a;
      return res;
    }
    g_vec3<T> operator/(const g_vec3<T> &rhs) const
    {
      g_vec3<T> res = *this;
      res.x /= rhs.x;
      res.y /= rhs.y;
      res.z /= rhs.z;
      return res;
    }
    g_vec3<T> operator/(T a) const
    {
      g_vec3<T> res = *this;
      res.x /= a;
      res.y /= a;
      res.z /= a;
      return res;
    }
    g_vec3<T> &operator+=(const g_vec3<T> &rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;
      return *this;
    }
    g_vec3<T> &operator-=(const g_vec3<T> &rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      z -= rhs.z;
      return *this;
    }
    g_vec3<T> &operator*=(const g_vec3<T> &rhs)
    {
      x *= rhs.x;
      y *= rhs.y;
      z *= rhs.z;
      return *this;
    }
    g_vec3<T> &operator*=(const T a)
    {
      x *= a;
      y *= a;
      z *= a;
      return *this;
    }
    g_vec3<T> &operator/=(const g_vec3<T> &rhs)
    {
      x /= rhs.x;
      y /= rhs.y;
      z /= rhs.z;
      return *this;
    }
    g_vec3<T> &operator/=(const T a)
    {
      x /= a;
      y /= a;
      z /= a;
      return *this;
    }
    inline bool operator==(const g_vec3<T> &rhs) const
    {
      return (x == rhs.x) && (y == rhs.y) && (z == rhs.z);
    }
    inline bool operator!=(const g_vec3<T> &rhs) const
    {
      return !(*this == rhs);
    }
    g_vec3<T> operator-()
    {
      return g_vec3<T>(-x, -y, -z);
    }
  };

  template <typename T>
  g_vec3<T> operator*(T a, const g_vec3<T> &rhs)
  {
    g_vec3<T> res = rhs;
    res.x *= a;
    res.y *= a;
    res.z *= a;
    return res;
  }

  template <typename T>
  T length(const g_vec3<T> &v)
  {
    return CppAD::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  }

  template <typename T>
  T dot(const g_vec3<T> &v1, const g_vec3<T> &v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  template <typename T>
  g_vec3<T> normalize(const g_vec3<T> &v)
  {
    T l = length(v);
    if (l > 1e-9)
      return v / l;
    else
      return g_vec3<T>(0,0,0);
  }

  template <typename T>
  g_vec3<T> cross(const g_vec3<T> &a, const g_vec3<T> &b)
  {
    g_vec3<T> res;

    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];

    return res;
  }

template <typename T>
  struct g_vec4
  {
    union
    {
      T data[4];
      struct
      {
        T x;
        T y;
        T z;
        T w;
      };
    };
    g_vec4<T>(T _x, T _y, T _z, T _w)
    {
      x = _x;
      y = _y;
      z = _z;
      w = _w;
    }
    g_vec4<T>(const g_vec4<T> &v) = default;
    g_vec4<T>(): g_vec4<T>(0,0,0,0)
    {

    }
    ~g_vec4<T>()
    {
      (&x)->T::~T();
      (&y)->T::~T();
      (&z)->T::~T();
      (&w)->T::~T();
    }
    g_vec4<T> &operator=(const g_vec4<T> &v) = default;
    g_vec4<T> &operator=(g_vec4<T> &&v) = default;

    T &operator[](int index)
    {
      return data[index];
    }
    const T &operator[](int index) const
    {
      return data[index];
    }
    g_vec4<T> operator+(const g_vec4<T> &rhs) const
    {
      g_vec4<T> res = *this;
      res.x += rhs.x;
      res.y += rhs.y;
      res.z += rhs.z;
      res.w += rhs.w;
      return res;
    }
    g_vec4<T> operator-(const g_vec4<T> &rhs) const
    {
      g_vec4<T> res = *this;
      res.x -= rhs.x;
      res.y -= rhs.y;
      res.z -= rhs.z;
      res.w -= rhs.w;
      return res;
    }
    g_vec4<T> operator*(const g_vec4<T> &rhs) const
    {
      g_vec4<T> res = *this;
      res.x *= rhs.x;
      res.y *= rhs.y;
      res.z *= rhs.z;
      res.w *= rhs.w;
      return res;
    }
    g_vec4<T> operator*(T a) const
    {
      g_vec4<T> res = *this;
      res.x *= a;
      res.y *= a;
      res.z *= a;
      res.w *= a;
      return res;
    }
    g_vec4<T> operator/(const g_vec4<T> &rhs) const
    {
      g_vec4<T> res = *this;
      res.x /= rhs.x;
      res.y /= rhs.y;
      res.z /= rhs.z;
      res.w /= rhs.w;
      return res;
    }
    g_vec4<T> operator/(T a) const
    {
      g_vec4<T> res = *this;
      res.x /= a;
      res.y /= a;
      res.z /= a;
      res.w /= a;
      return res;
    }
    g_vec4<T> &operator+=(const g_vec4<T> &rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;
      w += rhs.z;
      return *this;
    }
    g_vec4<T> &operator-=(const g_vec4<T> &rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      z -= rhs.z;
      w -= rhs.z;
      return *this;
    }
    g_vec4<T> &operator*=(const g_vec4<T> &rhs)
    {
      x *= rhs.x;
      y *= rhs.y;
      z *= rhs.z;
      w *= rhs.w;
      return *this;
    }
    g_vec4<T> &operator*=(const T a)
    {
      x *= a;
      y *= a;
      z *= a;
      w *= a;
      return *this;
    }
    g_vec4<T> &operator/=(const g_vec4<T> &rhs)
    {
      x /= rhs.x;
      y /= rhs.y;
      z /= rhs.z;
      w /= rhs.z;
      return *this;
    }
    g_vec4<T> &operator/=(const T a)
    {
      x /= a;
      y /= a;
      z /= a;
      w /= a;
      return *this;
    }
    inline bool operator==(const g_vec4<T> &rhs) const
    {
      return (x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (w == rhs.w);
    }
    inline bool operator!=(const g_vec4<T> &rhs) const
    {
      return !(*this == rhs);
    }
    g_vec4<T> operator-()
    {
      return g_vec4<T>(-x, -y, -z, -w);
    }
  };

  template <typename T>
  g_vec4<T> operator*(T a, const g_vec4<T> &rhs)
  {
    g_vec4<T> res = rhs;
    res.x *= a;
    res.y *= a;
    res.z *= a;
    res.w *= a;
    return res;
  }
    
  template <typename T>
  T length(const g_vec4<T> &v)
  {
    return CppAD::sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
  }

  template <typename T>
  T dot(const g_vec4<T> &v1, const g_vec4<T> &v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
  }

  template <typename T>
  g_vec4<T> normalize(const g_vec4<T> &v)
  {
    T l = length(v);
    if (l > 1e-9)
      return v / l;
    else
      return g_vec4<T>(0,0,0,0);
  }

  template <typename T>
  g_vec4<T> get_dvec4(const g_vec3<T> &xyz, T w)
  {
    return g_vec4<T>(xyz.x, xyz.y, xyz.z, w);
  }

  typedef CppAD::AD<float> dfloat;
  typedef g_vec2<dfloat> dvec2;
  typedef g_vec3<dfloat> dvec3;
  typedef g_vec4<dfloat> dvec4;
  typedef std::array<dfloat, 12> dmat43;//4 vec3 

  //matrices
  dmat43 get_mat43(float *data);
  dmat43 get_mat43(const dvec3 &a, const dvec3 &b, const dvec3 &c, const dvec3 &tr);
  dmat43 ident();
  dmat43 translate(const dmat43 &input_mat, const dfloat x, const dfloat y, const dfloat z);
  dmat43 translate(const dmat43 &mat, const dvec3 &tr);
  dmat43 scale(const dmat43 &input_mat, const dfloat x, const dfloat y, const dfloat z);
  dmat43 scale(const dmat43 &mat, const dvec3 &sc);
  dmat43 rotate(const dmat43 &input_mat, const dvec3 &axis, dfloat angle);
  dmat43 mul(const dmat43 &a, const dmat43 &b);
  dvec3 mulp(const dmat43 &mat, const dvec3 &vec);
  void mulp(const dmat43 &mat, dfloat &x, dfloat &y, dfloat &z);
  dvec3 mulv(const dmat43 &mat, const dvec3 &vec);
  void mulv(const dmat43 &mat, dfloat &x, dfloat &y, dfloat &z);
  dvec3 mul4(const dmat43 &mat, const dvec4 &vec);
  dmat43 transpose3x3(const dmat43 &b);
  dmat43 transposedInverse3x3(const dmat43 &m);
  dmat43 inverse3x4(const dmat43 &m);
}