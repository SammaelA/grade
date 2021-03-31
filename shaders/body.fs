#version 330

in vec4 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;

bool need_coord;
uniform vec3 main_color;
void main(void) 
{
  if (need_coord)
    fragColor = vec4(main_color,1);
  else 
    fragColor = vec4(ex_Tex.xyz,1);
}
