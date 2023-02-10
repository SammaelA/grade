#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
uniform int radius;
uniform float alpha_thr;
out vec4 fragColor;

void main(void) 
{
  ivec2 tc = ivec2(tex_size*ex_Tex);
  float a_sum = 0;
  int a_cnt = 0;
  vec3 col_sum = vec3(0, 0, 0);
  for (int i=-radius; i<=radius; i++)
  {
    for (int j=-radius; j<=radius; j++)
    {
      vec4 p = texelFetch(tex, ivec2(tc.x + i, tc.y + j), 0);
      if (p.a > 0)
        a_cnt++;
      float a = p.a/max(i*i+j*j, 1);
      a_sum += a;
      col_sum += a*p.xyz;
    }
  }
  vec4 self_val = texelFetch(tex, tc, 0);
  if (self_val.a == 0)
  {
    col_sum /= a_sum;
    if (a_cnt < alpha_thr*(2*radius + 1)*(2*radius + 1))
      fragColor = vec4(0,0,0,0);
    else
      fragColor = vec4(col_sum,a_sum);
  }
  else
  {
    if (a_cnt <= 2*(2*radius + 1))
      fragColor = vec4(0,0,0,0);
    else
    fragColor = vec4(self_val.xyz, 1);
  }
}