#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform sampler2D tex_reference;
uniform vec2 tex_size;
out vec4 fragColor;

void main(void) 
{
  ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
  float diff = texelFetch(tex, pixel_coords, 0).x - texelFetch(tex_reference, pixel_coords, 0).x;
  fragColor = vec4(diff,diff,diff,1);
}