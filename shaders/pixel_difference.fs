#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex1;
uniform sampler2D tex2;
void main(void) 
{
  vec4 c1 = textureLod(tex1, ex_Tex, 0);
  vec4 c2 = textureLod(tex2, ex_Tex, 0);
  float diff = length(c1.a*c1.rgb - c2.a*c2.rgb);
  fragColor = vec4(diff,diff,diff,1);
}