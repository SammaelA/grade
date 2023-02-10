#version 330

in vec2 ex_Tex;
uniform sampler2D tex_1;
uniform sampler2D tex_2;
uniform float a;
uniform float b;
uniform int op;

out vec4 fragColor;

#define ADD 0
#define MUL 1
#define DIV 2

void main(void) 
{
  vec4 c1 = texture(tex_1,ex_Tex);
  vec4 c2 = texture(tex_2,ex_Tex);
  if (op == ADD)
    fragColor = a*c1 + b*c2;
  else if (op == MUL)
    fragColor = a*b*c1*c2;
  else 
    fragColor = a*b*c1/max(c2, vec4(1.0/255, 1.0/255, 1.0/255, 1.0/255));
}