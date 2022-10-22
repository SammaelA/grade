#version 330

in vec3 ex_Tex;

uniform sampler2D mask;

layout (location = 0) out vec4 fragUV;

void main(void) 
{
  fragUV.b = texture(mask,ex_Tex.xy).r;
  fragUV.rg = ex_Tex.xy;
  fragUV.a = 1
}
