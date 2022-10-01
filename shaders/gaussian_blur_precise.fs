#version 430

layout(std140, binding=0) readonly buffer KernelBuf
{
  vec4 kernel[];
};

in vec2 ex_Tex;
uniform sampler2D tex;
uniform int pass;
uniform int steps;
uniform vec2 tex_size;
out vec4 fragColor;

void main(void) 
{
    vec4 res = vec4(0, 0, 0, 0);
    ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
    if (pass == 0)
    {
      for (int i = -steps; i <= steps; i++)
        res += texelFetch(tex, clamp(pixel_coords + ivec2(i, 0),ivec2(0,0),ivec2(tex_size-1)), 0) * kernel[abs(i)].x;
    }
    else
    {
      for (int i = -steps; i <= steps; i++)
        res += texelFetch(tex, clamp(pixel_coords + ivec2(0, i),ivec2(0,0),ivec2(tex_size-1)), 0) * kernel[abs(i)].x;
    }
    fragColor = res;
}
