#version 330

in vec4 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;


void main(void) 
{
  fragColor = ex_Tex;
}
