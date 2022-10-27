#version 330

in vec2 ex_Tex;
uniform sampler2D tex_mask;
uniform int search_radius;
uniform vec2 tex_size;
out vec4 fragColor;

float is_hole(ivec2 pixel_coords, ivec2 search_step)
{
  ivec2 a_step = ivec2(-search_step.y, search_step.x);
  for (int i=1;i<=search_radius;i++)
  {
    float v1 = texelFetch(tex_mask, pixel_coords + i*search_step, 0).x +
               texelFetch(tex_mask, pixel_coords + i*search_step - a_step, 0).x +
               texelFetch(tex_mask, pixel_coords + i*search_step + a_step, 0).x;
    float v2 = texelFetch(tex_mask, pixel_coords - i*search_step, 0).x +
               texelFetch(tex_mask, pixel_coords - i*search_step - a_step, 0).x +
               texelFetch(tex_mask, pixel_coords - i*search_step + a_step, 0).x;
    if (v1 > 1.5 && v2 > 1.5)
      return 1;
  }
  return 0;
}

void main(void) 
{
  vec4 res = vec4(0, 0, 0, 1);
  ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
  vec3 color = texelFetch(tex_mask, pixel_coords, 0).xyz;
  res.xyz = color;
  if (color.x < 0.9)
    res.xyz = vec3(1,1,1)*(is_hole(pixel_coords, ivec2(1,0)) + is_hole(pixel_coords, ivec2(0,1)));
  fragColor = res;
}