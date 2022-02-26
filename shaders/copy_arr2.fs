#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
void main(void) 
{
  fragColor = vec4(textureLod(tex,vec3(ex_Tex,layer),0).xyz,1);
}