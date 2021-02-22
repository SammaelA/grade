#version 430

in vec3 tc_a;
in vec3 tc_b;
in vec3 tc_t;
in vec3 q_abt;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec2 a_mult;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform sampler2D noise;
uniform vec4 screen_size;

void main(void) 
{
  vec2 noise_pos = vec2(fract(3*gl_FragCoord.x*screen_size.z), fract(3*gl_FragCoord.y*screen_size.w));
  float ns = texture(noise,noise_pos).x;
  fragColor = ns <= q_abt.x ? texture(tex,tc_a) : (ns <= q_abt.x + q_abt.y ? texture(tex,tc_b) : texture(tex,tc_t));
  if (fragColor.a<0.33)
    discard;
  
  if (a_mult.x < 5 && a_mult.x < abs(a_mult.y)*(-0.5*(sign(a_mult.y) - 1) + sign(a_mult.y)*ns))
    discard;
  
  fragColor.rgb /= fragColor.a;
}
