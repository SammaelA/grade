#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
void main(void) 
{
  fragColor = (ex_Tex.x >= 0 && ex_Tex.x <= 1 && ex_Tex.y >= 0 && ex_Tex.y <= 1) ? 
               vec4(textureLod(tex,vec3(ex_Tex,layer),0)) : vec4(0,0,0,0);
}