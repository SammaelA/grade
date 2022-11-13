#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
uniform int radius;
uniform vec3 base_color;
out vec4 fragColor;

void main(void) 
{
  ivec2 tc = ivec2(tex_size*ex_Tex);
  vec4 val = texelFetch(tex, tc, 0);
  if (val.a == 0)
  {
    ivec2 e_dir = ivec2(0,0);
    for (int i=1;i<radius;i++)
    {
      if (texelFetch(tex, tc+ivec2(i,0), 0).a > 0)
      {
        e_dir = ivec2(i, 0);
        break;
      }
      if (texelFetch(tex, tc+ivec2(-i,0), 0).a > 0)
      {
        e_dir = ivec2(-i, 0);
        break;
      }
      if (texelFetch(tex, tc+ivec2(0,i), 0).a > 0)
      {
        e_dir = ivec2(0, i);
        break;
      }
      if (texelFetch(tex, tc+ivec2(0,-i), 0).a > 0)
      {
        e_dir = ivec2(0, -i);
        break;
      }
    }
    if (e_dir == vec2(0,0))
      fragColor = vec4(base_color, 1);
    else
      fragColor = vec4(texelFetch(tex, tc+2*e_dir, 0).xyz, 1);
  }
  else
    fragColor = vec4(val.xyz, 1);
}