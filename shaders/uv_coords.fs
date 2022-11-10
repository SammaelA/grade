#version 330

in vec3 ex_Tex;
in vec2 ex_Pos;

uniform sampler2D mask;

out vec4 fragUV;

void main(void) 
{
  fragUV.b = 1;
  fragUV.rg = ex_Tex.xy;
  fragUV.a = 1;
}
