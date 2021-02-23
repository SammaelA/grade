#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec2 a_mult;
out vec4 fragColor;

uniform sampler2DArray tex;
uniform sampler2D noise;
uniform vec4 screen_size;
float gradientNoise(float x, float y)
{
  float f = 0.06711056f * x + 0.00583715f * y;
  return fract(52.9829189f * fract(f));
}
void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (fragColor.a<0.33)
    discard;
  fragColor.rgb /= fragColor.a;
  vec2 noise_pos = vec2(fract(25*gl_FragCoord.x*screen_size.z), fract(25*gl_FragCoord.y*screen_size.w));
  float ns = texture(noise,noise_pos).x;
  ns = gradientNoise(float(gl_FragCoord.x),float(gl_FragCoord.y));
  if ((a_mult.y > 0.1 && ns > a_mult.x) || (a_mult.y < -0.1 && ns <  1  - a_mult.x))
    discard;
}
