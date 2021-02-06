#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec2 a_mult;
out vec4 fragColor;

uniform sampler2DArray tex;
uniform sampler2D noise;
uniform vec4 screen_size;
void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (fragColor.a<0.33)
    discard;
  vec2 noise_pos = vec2(fract(3*gl_FragCoord.x*screen_size.z), fract(3*gl_FragCoord.y*screen_size.w));
  float ns = texture(noise,noise_pos).x;
  if (a_mult.x < 5 && a_mult.x < abs(a_mult.y)*(-0.5*(sign(a_mult.y) - 1) + sign(a_mult.y)*ns))
    discard;
}
