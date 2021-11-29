#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;

void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  fragColor.xyz = fragColor.a > 0 ? fragColor.xyz : vec3(0,0,0);
  fragColor.a = fragColor.a > 0 ? 1 : 0;
}
