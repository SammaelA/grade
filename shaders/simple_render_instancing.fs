#version 330

in vec2 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;


void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (fragColor.a<0.33)
    discard;
}
