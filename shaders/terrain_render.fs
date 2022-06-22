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
uniform sampler2D debug_tex;

uniform int debug_render = 1;//0 - no debug, 1 - show grid, 2 - show debug tex, 3 - show both
uniform vec4 grid_params = vec4(0,0,50,50);//(start_x, start_y, size_x, size_y)
uniform vec4 grid_render_params = vec4(0.7,0.85,1,0.05);//(r,g,b,thickness as part of grid size)
uniform vec4 debug_tex_scale = vec4(0,0,0.01, 0.01);//tc = (worldPos.xz - scale.xy)*scale.zw;
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
  fragNormal = vec4(ex_Normal.xyz,1);
  fragViewPos = vec4(ex_FragPosView.xyz,1);
  fragWorldPos = vec4(ex_FragPos,1);

  if (debug_render % 2 == 1)
  {
    vec2 grd_pos = (ex_FragPos.xz - grid_params.xy)/grid_params.zw;
    grd_pos = abs(grd_pos - vec2(ivec2(grd_pos)));
    if ((min(min(grd_pos.x, grd_pos.y), min(1 - grd_pos.x, 1 - grd_pos.y)) < 0.5*grid_render_params.w) &&
        length(fragViewPos) < 1000)
    {
      fragColor.xyz += grid_render_params.xyz;
      fragColor.w = 1;
    }
  }
  if (debug_render / 2 % 2 == 1)
  {
    vec2 dbg_uv = (ex_FragPos.xz - debug_tex_scale.xy)*debug_tex_scale.zw;
    fragColor.xyz = (0.1 + 0.8*texture(debug_tex,dbg_uv).xyz);
  }
}
