typedef float my_float;

namespace u_g
{
  struct vec2
  {
    my_float x;
    my_float y;
  };

  struct vec3
  {
    my_float x;
    my_float y;
    my_float z;
  };

  struct vec4
  {
    my_float x;
    my_float y;
    my_float z;
    my_float w;
  };

  struct mat3
  {
    my_float data[3][3];
  };

  vec2 operator-(const vec2 &a)
  {
    return vec2{-a.x, -a.y};
  }

  vec2 operator+(const vec2 &a, const vec2 &b)
  {
    return vec2{a.x + b.x, a.y + b.y};
  }

  vec3 operator-(const vec3 &a)
  {
    return vec3{-a.x, -a.y, -a.z};
  }

  vec3 operator+(const vec3 &a, const vec3 &b)
  {
    return vec3{a.x + b.x, a.y + b.y, a.z + b.z};
  }

  vec4 operator-(const vec4 &a)
  {
    return vec4{-a.x, -a.y, -a.z, -a.w};
  }

  vec4 operator+(const vec4 &a, const vec4 &b)
  {
    return vec4{a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
  }
}