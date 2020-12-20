#version 330

out vec2 ex_Tex;

out vec4 fragColor;
uniform sampler2D tex;
uniform vec4 screen_size;
void main(void) 
{
  ivec2 pos = ivec2(ex_Tex.x*screen_size.x,ex_Tex.y*screen_size.y);
  vec4 fc0 = texelFetch(tex,pos + ivec2(0,0),0);
  vec4 fc1 = texelFetch(tex,pos + ivec2(0,1),0);
  vec4 fc2 = texelFetch(tex,pos + ivec2(1,0),0);
  vec4 fc3 = texelFetch(tex,pos + ivec2(1,1),0);
  fragColor.rgb = (fc0.rgb*fc0.a + fc1.rgb*fc1.a + fc2.rgb*fc2.a + fc3.rgb*fc3.a)/max(1,fc0.a + fc1.a + fc2.a +fc3.a);
  fragColor.a = max(max(fc0.a,fc1.a),max(fc2.a,fc3.a));
  fragColor.rgb = fragColor.aaa;
}