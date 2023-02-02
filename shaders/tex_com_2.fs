#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex1;
uniform sampler2D tex2;
uniform sampler2D mask1;
uniform sampler2D mask2;
uniform int x_sym;
uniform int y_sym;

void main(void) 
{
  float cnt = 0;
  fragColor = texture(tex1,ex_Tex) + texture(tex2,ex_Tex);
  float norm1 = 0;
  float norm2 = 0;
  for(float i = 0; i <= abs(x_sym); i += 1)
  {
    if (x_sym == 0) return;
    if (y_sym == 0) return;
    for(float j = 0; j <= abs(y_sym); j += 1)
    {
      norm1 += texture(mask1,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))).r;
      norm2 += texture(mask2,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))).r;
    }
  }
  if (norm1 == 0) norm1 = 1;
  if (norm2 == 0) norm2 = 1;
  fragColor = vec4(0, 0, 0, 0);
  for(float i = 0; i <= abs(x_sym); i += 1)
  {
    for(float j = 0; j <= abs(y_sym); j += 1)
    {
      fragColor += texture(tex1,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))) * texture(mask1,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))).r / norm1;
      fragColor += texture(tex2,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))) * texture(mask2,ex_Tex + vec2(i / abs(x_sym), j / abs(y_sym))).r / norm2;
    }
  }
}
