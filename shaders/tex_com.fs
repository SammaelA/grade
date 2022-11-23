#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform int x_sym;
uniform int y_sym;

void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (fragColor.a == 0)
  {
    float cnt = 0;
    fragColor = vec4(0, 0, 0, 0);
    for(float i = 0; i <= abs(x_sym); i += 1)
    {
      if (x_sym == 0) break;
      if (y_sym == 0) break;
      for(float j = 0; j <= abs(y_sym); j += 1)
      {
        if (texture(tex,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))).a > 0)
        {
          fragColor += texture(tex,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym)));
          cnt += 1;
        }
      }
    }
    fragColor = fragColor / cnt;
  }
}
