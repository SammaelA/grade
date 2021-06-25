#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec4 ex_FragPosView;
in vec4 FragPosLightSpace;

uniform sampler2D terrain;

layout (location = 0) out vec4 fragColor;
layout (location = 1) out vec4 fragNormal;
layout (location = 2) out vec4 fragViewPos;
layout (location = 3) out vec4 fragWorldPos;

void main(void) 
{
  fragColor = texture(terrain,ex_Tex.xy);
  fragNormal = vec4(ex_Normal.xyz,0);
  fragViewPos = vec4(ex_FragPosView.xyz,1);
  fragWorldPos = vec4(ex_FragPos,1);
}
