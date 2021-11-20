#version 330

in vec2 ex_Tex;

out vec4 fragColor;

uniform sampler2D target;
uniform sampler2D prevTarget;
uniform float weight;
void main(void) 
{
  fragColor = weight*texture(prevTarget,ex_Tex) + (1 - weight)*texture(target,ex_Tex);
}
