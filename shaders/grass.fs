#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex1;
uniform sampler2D tex2;
void main(void) 
{
  fragColor = texture(tex1,ex_Tex.xy);
  if (fragColor.a<0.33)
    discard;
 fragColor.rgb/fragColor.a;
}