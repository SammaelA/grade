#version 330

in vec2 ex_Tex;

out vec4 fragColor;

uniform sampler2D tex;
uniform int radius;
uniform vec2 tex_size;

void main(void) 
{
  float c = 1;
  for (int i=0; i<radius; i++)
  {
    for (int j=0; j<radius; j++)
    {
      c = min(c, texelFetch(tex, ivec2(ex_Tex*tex_size) + ivec2(i,j), 0).a);
    }
  }
  fragColor = vec4(c,c,c,c);
}
