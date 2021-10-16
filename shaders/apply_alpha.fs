#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
void main(void) 
{
  vec4 col = texture(tex,vec3(ex_Tex,layer));
  fragColor = vec4(col.xyz*col.a,1);
}
