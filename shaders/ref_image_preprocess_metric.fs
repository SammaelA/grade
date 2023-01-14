#version 330

in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D sil;
uniform float max_rad;
uniform float size_x;
uniform float size_y;

void main(void) 
{
  fragColor = texture(sil,ex_Tex);
  if (fragColor.x > 0.9)
  {
    fragColor = vec4(1, 1, 1, 1);
  }
  else
  {
    fragColor = vec4(0, 0, 0, 1);
    int brk = 0;
    for (float r = 2; r < max_rad; r=r+1)
    {
      float sectors = max(4, r);
      for (float i = 0; i < sectors; i=i+1)
      {
        #define PI 3.14159
        float s = 2*PI*i/sectors;
        float hx = r*sin(s);
        float hy = r*cos(s);
        float k = max(2, sqrt(hx*hx + hy*hy));
        if (texture(sil,ex_Tex + vec2(hx/size_x, hy/size_x)).x > 0)
        {
          fragColor = vec4(1/k, 1/k, 1/k, 1);
          return;
        }
      }
    }
  }
}
