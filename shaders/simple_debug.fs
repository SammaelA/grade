#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec4 ex_FragPosView;

uniform sampler2D tex;
uniform int type = 6;

layout (location = 0) out vec4 fragColor;
layout (location = 1) out vec4 fragNormal;
layout (location = 2) out vec4 fragViewPos;
layout (location = 3) out vec4 fragWorldPos;

void main(void) 
{
  fragColor.xyz = ex_Tex.xyz;
  fragColor.a = 1;

  fragNormal = vec4(ex_Normal.xyz,type);
  fragViewPos = vec4(ex_FragPosView.xyz,1);
  //fragWorldPos = vec4(ex_FragPos,1);

}