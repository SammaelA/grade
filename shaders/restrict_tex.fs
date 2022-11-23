#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform float x_sh;
uniform float y_sh;
uniform float x_sz;
uniform float y_sz;

void main(void) 
{
  fragColor = texture(tex,vec2(x_sh, y_sh) + ex_Tex * vec2(x_sz, y_sz));
}
