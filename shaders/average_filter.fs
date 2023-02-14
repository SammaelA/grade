#version 330

in vec2 ex_Tex;

out vec4 fragColor;

uniform sampler2D tex;
uniform int radius;
uniform vec2 tex_size;

void main(void) 
{
  vec4 c = vec4(0,0,0,0);
  for (int i=-radius; i<=radius; i++)
  {
    for (int j=-radius; j<=radius; j++)
    {
      c += texelFetch(tex, ivec2(ex_Tex*tex_size) + ivec2(i,j), 0);
    }
  }
  fragColor = c/float((2*radius + 1)*(2*radius + 1));
}
