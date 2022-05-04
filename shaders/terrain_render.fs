#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec4 ex_FragPosView;
in vec4 FragPosLightSpace;

uniform sampler2D grass1;
uniform sampler2D grass2;
uniform sampler2D rock;
uniform sampler2D perlin;

layout (location = 0) out vec4 fragColor;
layout (location = 1) out vec4 fragNormal;
layout (location = 2) out vec4 fragViewPos;
layout (location = 3) out vec4 fragWorldPos;

void main(void) 
{
  vec2 uvm1 = abs(fract(0.0015*ex_FragPos.xz)-0.5);
  vec2 uv = 2*abs(fract(ex_Tex.xy)-0.5);
  vec4 t1 = texture(grass1,uv*(1+uvm1));
  vec4 t2 = texture(grass2,uv*(2-uvm1.yx));
  vec4 t3 = texture(rock,uv);

  float p = texture(perlin,abs(fract(0.001*ex_FragPos.xz)-0.5)).x;
  float rmul = clamp(2*(abs(dFdx(ex_FragPos.y)) + abs(dFdy(ex_FragPos.y)))/(length(dFdx(ex_FragPos.xz)) + length(dFdy(ex_FragPos.xz))),0,1);
  rmul = rmul*rmul;
  fragColor = rmul*t3 + (1-rmul)*(p*t1 + (1-p)*t2);
  fragNormal = vec4(ex_Normal.xyz,0);
  fragViewPos = vec4(ex_FragPosView.xyz,1);
  fragWorldPos = vec4(ex_FragPos,1);
}
