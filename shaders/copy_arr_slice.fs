#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
uniform vec4 in_slice_tr;
uniform vec4 slice_tr;
void main(void) 
{
  vec2 tc_slice = in_slice_tr.xy + ex_Tex.xy*in_slice_tr.zw;
  fragColor = (tc_slice.x >= 0 && tc_slice.x <= 1 && tc_slice.y >= 0 && tc_slice.y <= 1) ? 
               vec4(textureLod(tex,vec3(slice_tr.xy + tc_slice*slice_tr.zw,layer),0)) : vec4(0,0,0,0);
}