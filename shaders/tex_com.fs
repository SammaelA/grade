#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform int sym;
uniform int is_x;
uniform ivec2 tex_size;

void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (is_x != 0)
  {
    if (fragColor.a == 0)
    {
      float cnt = 0;
      fragColor = vec4(0, 0, 0, 0);
      for(float i = 0; i < abs(sym); i += 1)
      {
        if (sym > 0)
        {
          if (texture(tex,ex_Tex + vec2(i / sym, 0)).a > 0)
          {
            fragColor += texture(tex,ex_Tex + vec2(i / sym, 0));
            cnt += 1;
          }
        }
        else if (sym < 0)
        {
          if (texture(tex,-ex_Tex + vec2(i / sym, 0)).a > 0)
          {
            fragColor += texture(tex,vec2(-ex_Tex.x, ex_Tex.y) + vec2(i / sym, 0));
            cnt += 1;
          }
        }
      }
      fragColor = fragColor / cnt;
    }
  }
  else
  {
    if (fragColor.a == 0)
    {
      float cnt = 0;
      fragColor = vec4(0, 0, 0, 0);
      for(float i = 0; i < abs(sym); i += 1)
      {
        if (sym > 0)
        {
          if (texture(tex,ex_Tex + vec2(0, i / sym)).a > 0)
          {
            fragColor += texture(tex,ex_Tex + vec2(0, i / sym));
            cnt += 1;
          }
        }
        else if (sym < 0)
        {
          if (texture(tex,-ex_Tex + vec2(0, i / sym)).a > 0)
          {
            fragColor += texture(tex,vec2(ex_Tex.x, -ex_Tex.y) + vec2(0, i / sym));
            cnt += 1;
          }
        }
      }
      fragColor = fragColor / cnt;
    }
  }
}
