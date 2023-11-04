#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;

void main(void) 
{
  fragColor = vec4(texture(tex,ex_Tex).x > 0.5);
}
