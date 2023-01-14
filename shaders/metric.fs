#version 330

in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D sil;
uniform float max_rad;
uniform float size_x;
uniform float size_y;

void main(void) 
{
  fragColor = vec4(1, 1, 1, 1);
  if (texture(sil,ex_Tex).x > 0)
  {
    fragColor = vec4(0, 2 * max_rad, texture(sil,ex_Tex).x, 1);
  }
  else
  {
    fragColor = vec4(1, 1, 0, 1);
    int brk = 0;
    for (float i = 1; i < max_rad && brk == 0; i=i+1)
    {
      for (float j = 0; j <= i && brk == 0; j=j+1)
      {
        float k = sqrt(i * i - j * j);
        if (texture(sil,ex_Tex + vec2(j / size_x, k / size_y)).x > 0)
        {
          fragColor = vec4(i / max_rad, max_rad / i, 0, 1);
          brk = 1;
        }
        if (texture(sil,ex_Tex + vec2(-j / size_x, k / size_y)).x > 0)
        {
          fragColor = vec4(i / max_rad, max_rad / i, 0, 1);
          brk = 1;
        }
        if (texture(sil,ex_Tex + vec2(j / size_x, -k / size_y)).x > 0)
        {
          fragColor = vec4(i / max_rad, max_rad / i, 0, 1);
          brk = 1;
        }
        if (texture(sil,ex_Tex + vec2(-j / size_x, -k / size_y)).x > 0)
        {
          fragColor = vec4(i / max_rad, max_rad / i, 0, 1);
          brk = 1;
        }
      }
    } 
  }
}
