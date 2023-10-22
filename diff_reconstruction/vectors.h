#pragma once
#include <vector>
#include <array>
#include <cmath>

namespace dgen
{
  //a set of functions should be defined for T to use g_vec2/3/4<T> 
  //these functions are sqrt, sin, cos
  //here are required definitions for float and double

  inline float sqrt(float x) { return ::sqrtf(x); }
  inline float sin(float x)  { return ::sinf(x);  }
  inline float cos(float x)  { return ::cosf(x);  }

  inline double sqrt(double x) { return ::sqrt(x); }
  inline double sin(double x)  { return ::sin(x);  }
  inline double cos(double x)  { return ::cos(x);  }

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
    g_vec2<T> operator-() const
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
    return sqrt(v.x * v.x + v.y * v.y);
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
    g_vec3<T> operator-() const
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
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  }

  template <typename T>
  T dot(const g_vec3<T> &v1, const g_vec3<T> &v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  template <typename T>
  g_vec3<T> normalize_with_default(const g_vec3<T> &v, const g_vec3<T> &def_val)
  {
    T l = length(v);
    if (l > 1e-9)
      return v / l;
    else
      return def_val;
  }

  template <typename T>
  g_vec3<T> normalize(const g_vec3<T> &v)
  {
    return normalize_with_default(v, g_vec3<T>(0,0,0));
  }

  template <typename T>
  g_vec3<T> normalize_unsafe(const g_vec3<T> &v)
  {
    return v / length(v);
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
    g_vec4<T> operator-() const
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
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
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

  template <typename T>
  using g_mat43 = std::array<T, 12>;

  //matrices
  template <typename T>
  g_mat43<T> get_mat43(float *data)
  {
    g_mat43<T> mat;

    for (int i = 0; i < 12; i++)
      mat[i] = data[i];
    
    return mat;
  }

  template <typename T>
  g_mat43<T> get_mat43(const g_vec3<T> &a, const g_vec3<T> &b, const g_vec3<T> &c, const g_vec3<T> &tr)
  {
    g_mat43<T> mat;

    mat[0] = a[0];
    mat[1] = a[1];
    mat[2] = a[2];

    mat[3] = b[0];
    mat[4] = b[1];
    mat[5] = b[2];

    mat[6] = c[0];
    mat[7] = c[1];
    mat[8] = c[2];

    mat[9] = tr[0];
    mat[10] = tr[1];
    mat[11] = tr[2];

    return mat;
  }

  template <typename T>
  g_mat43<T> ident()
  {
    g_mat43<T> mat;

    mat[0] = 1;
    mat[1] = 0;
    mat[2] = 0;

    mat[3] = 0;
    mat[4] = 1;
    mat[5] = 0;

    mat[6] = 0;
    mat[7] = 0;
    mat[8] = 1;

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;

    return mat;
  }

  template <typename T>
  g_mat43<T> mul(const g_mat43<T> &a, const g_mat43<T> &b)
  {
    g_mat43<T> mat;
    mat[0] = a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
    mat[1] = a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
    mat[2] = a[2] * b[0] + a[5] * b[1] + a[8] * b[2];

    mat[3] = a[0] * b[3] + a[3] * b[4] + a[6] * b[5];
    mat[4] = a[1] * b[3] + a[4] * b[4] + a[7] * b[5];
    mat[5] = a[2] * b[3] + a[5] * b[4] + a[8] * b[5];

    mat[6] = a[0] * b[6] + a[3] * b[7] + a[6] * b[8];
    mat[7] = a[1] * b[6] + a[4] * b[7] + a[7] * b[8];
    mat[8] = a[2] * b[6] + a[5] * b[7] + a[8] * b[8];

    mat[9] = a[0] * b[9] + a[3] * b[10] + a[6] * b[11] + a[9];
    mat[10] = a[1] * b[9] + a[4] * b[10] + a[7] * b[11] + a[10];
    mat[11] = a[2] * b[9] + a[5] * b[10] + a[8] * b[11] + a[11];
    
    return mat;
  }

  template <typename T>
  g_mat43<T> translate(const g_mat43<T> &input_mat, const T x, const T y, const T z)
  {
    g_mat43<T> mat;
    
    mat[0] = 1;
    mat[1] = 0;
    mat[2] = 0;

    mat[3] = 0;
    mat[4] = 1;
    mat[5] = 0;

    mat[6] = 0;
    mat[7] = 0;
    mat[8] = 1;

    mat[9] = x;
    mat[10] = y;
    mat[11] = z;

    return mul(input_mat, mat);
  }

  template <typename T>
  g_mat43<T> translate(const g_mat43<T> &mat, const g_vec3<T> &tr)
  {
    return translate(mat, tr[0], tr[1], tr[2]);
  }

  template <typename T>
  g_mat43<T> scale(const g_mat43<T> &input_mat, const T x, const T y, const T z)
  {
    g_mat43<T> mat;

    mat[0] = x;
    mat[1] = 0;
    mat[2] = 0;

    mat[3] = 0;
    mat[4] = y;
    mat[5] = 0;

    mat[6] = 0;
    mat[7] = 0;
    mat[8] = z;

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;

    return mul(input_mat, mat);
  }

  template <typename T>
  g_mat43<T> scale(const g_mat43<T> &mat, const g_vec3<T> &sc)
  {
    return scale(mat, sc[0], sc[1], sc[2]);
  }

  template <typename T>
  g_mat43<T> rotate(const g_mat43<T> &input_mat, const g_vec3<T> &axis, T angle)
  {
    g_mat43<T> mat;
    g_vec3<T> u = normalize(axis);
    T sn = sin(angle);
    T cs = cos(angle);

    mat[0] = cs + u[0] * u[0] * (1 - cs);
    mat[1] = u[0] * u[1] * (1 - cs) + u[2] * sn;
    mat[2] = u[0] * u[2] * (1 - cs) - u[1] * sn;

    mat[3] = u[0] * u[1] * (1 - cs) - u[2] * sn;
    mat[4] = cs + u[1] * u[1] * (1 - cs);
    mat[5] = u[1] * u[2] * (1 - cs) + u[0] * sn;

    mat[6] = u[0] * u[2] * (1 - cs) + u[1] * sn;
    mat[7] = u[1] * u[2] * (1 - cs) - u[0] * sn;
    mat[8] = cs + u[2] * u[2] * (1 - cs);

    mat[9] = 0;
    mat[10] = 0;
    mat[11] = 0;

    return mul(input_mat, mat);
  }

  template <typename T>
  g_vec3<T> mulp(const g_mat43<T> &mat, const g_vec3<T> &vec) // mul (vec.x, vec.y, vec.z, 1)
  {
    g_vec3<T> res;
    res[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2] + mat[9];
    res[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2] + mat[10];
    res[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2] + mat[11];
    return res;
  }

  template <typename T>
  void mulp(const g_mat43<T> &mat, T &x, T &y, T &z) // mul (vec.x, vec.y, vec.z, 1)
  {
    T x1 = mat[0] * x + mat[3] * y + mat[6] * z + mat[9];
    T y1 = mat[1] * x + mat[4] * y + mat[7] * z + mat[10];
    T z1 = mat[2] * x + mat[5] * y + mat[8] * z + mat[11];

    x = x1;
    y = y1;
    z = z1;
  }

  template <typename T>
  g_vec3<T> mulv(const g_mat43<T> &mat, const g_vec3<T> &vec) // mul (vec.x, vec.y, vec.z, 0)
  {
    g_vec3<T> res;
    res[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2];
    res[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2];
    res[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2];
    return res;
  }

  template <typename T>
  void mulv(const g_mat43<T> &mat, T &x, T &y, T &z) // mul (vec.x, vec.y, vec.z, 0)
  {
    T x1 = mat[0] * x + mat[3] * y + mat[6] * z;
    T y1 = mat[1] * x + mat[4] * y + mat[7] * z;
    T z1 = mat[2] * x + mat[5] * y + mat[8] * z;

    x = x1;
    y = y1;
    z = z1;
  }

  template <typename T>
  g_vec3<T> mul4(const g_mat43<T> &mat, const g_vec4<T> &vec)
  {
    g_vec3<T> res;
    res[0] = mat[0] * vec[0] + mat[3] * vec[1] + mat[6] * vec[2] + mat[9] * vec[3];
    res[1] = mat[1] * vec[0] + mat[4] * vec[1] + mat[7] * vec[2] + mat[10] * vec[3];
    res[2] = mat[2] * vec[0] + mat[5] * vec[1] + mat[8] * vec[2] + mat[11] * vec[3];
    return res;
  }

  template <typename T>
  g_mat43<T> transpose3x3(const g_mat43<T> &b)
  {
    g_mat43<T> mat;
    mat[0] = b[0];
    mat[1] = b[3];
    mat[2] = b[6];

    mat[3] = b[1];
    mat[4] = b[4];
    mat[5] = b[7];

    mat[6] = b[2];
    mat[7] = b[5];
    mat[8] = b[8];

    mat[9] = b[9];
    mat[10] = b[10];
    mat[11] = b[11];
    return mat;
  }
  
  template <typename T>
  g_mat43<T> transposedInverse3x3(const g_mat43<T> &m)
  {
    #define MAT(i,j) m[3*i + j] 
    g_mat43<T> result;
    T determinant =+MAT(0,0)*(MAT(1,1)*MAT(2,2)-MAT(2,1)*MAT(1,2))
                        -MAT(0,1)*(MAT(1,0)*MAT(2,2)-MAT(1,2)*MAT(2,0))
                        +MAT(0,2)*(MAT(1,0)*MAT(2,1)-MAT(1,1)*MAT(2,0));
    T invdet = 1/(determinant+1e-19);
    result[0] =  (MAT(1,1)*MAT(2,2)-MAT(2,1)*MAT(1,2))*invdet;
    result[3] = -(MAT(0,1)*MAT(2,2)-MAT(0,2)*MAT(2,1))*invdet;
    result[6] =  (MAT(0,1)*MAT(1,2)-MAT(0,2)*MAT(1,1))*invdet;
    result[1] = -(MAT(1,0)*MAT(2,2)-MAT(1,2)*MAT(2,0))*invdet;
    result[4] =  (MAT(0,0)*MAT(2,2)-MAT(0,2)*MAT(2,0))*invdet;
    result[7] = -(MAT(0,0)*MAT(1,2)-MAT(1,0)*MAT(0,2))*invdet;
    result[2] =  (MAT(1,0)*MAT(2,1)-MAT(2,0)*MAT(1,1))*invdet;
    result[5] = -(MAT(0,0)*MAT(2,1)-MAT(2,0)*MAT(0,1))*invdet;
    result[8] =  (MAT(0,0)*MAT(1,1)-MAT(1,0)*MAT(0,1))*invdet;
    result[9] =  m[9];
    result[10] = m[10];
    result[11] = m[11];

    return result;
  }

  template <typename T>
  g_mat43<T> inverse3x4(const g_mat43<T> &m)
  {
    #define MAT(i,j) m[3*i + j] 
    g_mat43<T> result;
    T determinant =+MAT(0,0)*(MAT(1,1)*MAT(2,2)-MAT(2,1)*MAT(1,2))
                        -MAT(0,1)*(MAT(1,0)*MAT(2,2)-MAT(1,2)*MAT(2,0))
                        +MAT(0,2)*(MAT(1,0)*MAT(2,1)-MAT(1,1)*MAT(2,0));
    T invdet = 1/(determinant+1e-19);
    result[0] =  (MAT(1,1)*MAT(2,2)-MAT(2,1)*MAT(1,2))*invdet;
    result[1] = -(MAT(0,1)*MAT(2,2)-MAT(0,2)*MAT(2,1))*invdet;
    result[2] =  (MAT(0,1)*MAT(1,2)-MAT(0,2)*MAT(1,1))*invdet;
    result[3] = -(MAT(1,0)*MAT(2,2)-MAT(1,2)*MAT(2,0))*invdet;
    result[4] =  (MAT(0,0)*MAT(2,2)-MAT(0,2)*MAT(2,0))*invdet;
    result[5] = -(MAT(0,0)*MAT(1,2)-MAT(1,0)*MAT(0,2))*invdet;
    result[6] =  (MAT(1,0)*MAT(2,1)-MAT(2,0)*MAT(1,1))*invdet;
    result[7] = -(MAT(0,0)*MAT(2,1)-MAT(2,0)*MAT(0,1))*invdet;
    result[8] =  (MAT(0,0)*MAT(1,1)-MAT(1,0)*MAT(0,1))*invdet;
    g_vec3<T> inv_tr;
    g_vec3<T> tr{m[9],m[10],m[11]};
    inv_tr = mulv(result, tr);
    result[9] =  -inv_tr[0];
    result[10] = -inv_tr[1];
    result[11] = -inv_tr[2];

    return result;
  }
}