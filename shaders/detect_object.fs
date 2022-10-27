#version 330

in vec2 ex_Tex;
uniform sampler2D tex_edges;
uniform sampler2D tex_color;
uniform float color_thr;
uniform vec2 tex_size;
out vec4 fragColor;

void main(void) 
{
  vec4 res = vec4(0, 0, 0, 1);
  ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
  //res.y = texelFetch(tex_edges, pixel_coords, 0).x;
  vec3 color = texelFetch(tex_color, pixel_coords, 0).xyz;
  vec3 bc1 = texelFetch(tex_color, ivec2(1, pixel_coords.y), 0).xyz;
  vec3 bc2 = texelFetch(tex_color, ivec2(pixel_coords.x, 1), 0).xyz;
  if (max(length(bc1 - color),length(bc2 - color)) > color_thr)
    res.xyz = vec3(1,1,1);
  fragColor = res;
}