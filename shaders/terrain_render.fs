#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;



void main(void) 
{
  fragColor = vec4(ex_Tex,1);
}
