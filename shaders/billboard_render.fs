#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;


void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  //fragColor.r += 0.4*(ex_Tex.z);
  //fragColor.b += 0.4*(1 - ex_Tex.z);
  if (fragColor.a<0.33)
    discard;
}
