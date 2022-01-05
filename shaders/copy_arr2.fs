#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
void main(void) 
{
  fragColor = vec4(texture(tex,vec3(ex_Tex,layer)).xyz,1);
  //fragColor = texture(tex,vec3(ex_Tex,layer));
  //fragColor = vec4(1,0,0,1);
  //fragColor = vec4(vec3(ex_Tex,layer),1);
}