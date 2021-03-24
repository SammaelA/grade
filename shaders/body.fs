#version 330

in vec4 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;


uniform vec3 main_color;
void main(void) 
{
  fragColor = vec4(main_color,1);
}
