#version 430

in vec3 tc_a;
in vec3 tc_b;
in vec3 tc_t;
in vec3 q_abt;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec2 a_mult;
in mat4 rot_m;
flat in uint model_id;
out vec4 fragColor;

uniform sampler2DArray tex;
uniform sampler2D noise;
uniform vec4 screen_size;
uniform int debug_model_id;
float gradientNoise(float x, float y)
{
  float f = 0.06711056f * x + 0.00583715f * y;
  return fract(52.9829189f * fract(f));
}
void main(void) 
{
  vec2 noise_pos = vec2(fract(25*gl_FragCoord.x*screen_size.z), fract(25*gl_FragCoord.y*screen_size.w));
  float ns_imp = texture(noise,noise_pos).x;
  float ns = gradientNoise(float(gl_FragCoord.x),float(gl_FragCoord.y));
  fragColor = ns_imp <= q_abt.x ? texture(tex,tc_a) : (ns_imp <= q_abt.x + q_abt.y ? texture(tex,tc_b) : texture(tex,tc_t));
  if (fragColor.a<0.33)
    discard;

  if ((a_mult.y > 0.1 && ns > a_mult.x) || (a_mult.y < -0.1 && ns <  1  - a_mult.x))
    discard;
  fragColor.rgb /= fragColor.a;
  if (int(model_id) == debug_model_id)
    fragColor = vec4(1,0,1,1);
}
