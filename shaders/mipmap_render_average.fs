#version 330

out vec2 ex_Tex;

out vec4 fragColor;
uniform sampler2D tex;
uniform vec4 screen_size;
void main(void) 
{
  ivec2 pos = ivec2(2*gl_FragCoord.x,2*gl_FragCoord.y);
  vec4 fc0 = texelFetch(tex,pos - ivec2(0,0),0);
  vec4 fc1 = texelFetch(tex,pos - ivec2(0,1),0);
  vec4 fc2 = texelFetch(tex,pos - ivec2(1,0),0);
  vec4 fc3 = texelFetch(tex,pos - ivec2(1,1),0);

  fragColor = 0.25*(fc0 + fc1 + fc2 + fc3);
}