#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform sampler2D mask;
uniform float x_sh;
uniform float y_sh;
uniform float x_sz;
uniform float y_sz;

void main(void) 
{
  if (texture(mask,ex_Tex).r > 0.8)
    fragColor = vec4(texture(mask,ex_Tex).r*texture(tex,ex_Tex).rgb, 1);
  else
    fragColor = vec4(0,0,0,1);
}
