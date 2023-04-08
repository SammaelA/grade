#version 330

in vec3 ex_Tex;
in vec2 ex_Pos;
out vec4 fragUV;
uniform sampler2D tex;
void main(void) 
{
  fragUV = texture(tex,ex_Tex.xy);
}
