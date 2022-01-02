#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform vec3 wood_color;
uniform vec3 leaves_color;
uniform vec3 background_color;
void main(void) 
{
  vec4 color = texture(tex,ex_Tex);

  fragColor = vec4(0,0,0,0);
  if (color.a > 0.5)
  {
      float a = dot(wood_color - color.xyz, wood_color - color.xyz)/dot(wood_color + color.xyz, wood_color + color.xyz);
      float b = dot(leaves_color - color.xyz, leaves_color - color.xyz)/dot(leaves_color + color.xyz, leaves_color + color.xyz);
      float c = dot(background_color - color.xyz, background_color - color.xyz);
      float g_part = color.g/(color.r + color.b + 1e-4);
      if (c < 1e-4)
        fragColor = vec4(0,0,0,0);
      else if (g_part < 0.5)
        fragColor = vec4(1,0,0,1);
      else
        fragColor = vec4(0,1,0,1);    
  }
}