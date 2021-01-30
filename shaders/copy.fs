#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform vec4 multiplier;

void main(void) 
{
  fragColor = multiplier*texture(tex,ex_Tex);
}
