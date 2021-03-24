#version 330

in vec4 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2D tex;
uniform sampler2DArray tex_arr;
uniform bool need_tex;
uniform bool need_arr_tex;
uniform bool need_coord;
uniform bool need_color;
uniform vec3 main_color;
uniform int slice;
void main(void) 
{
  vec4 color = vec4(0,0,0,1);
  if (need_tex)
    color += texture(tex,ex_Tex.xy);
  else if (need_arr_tex)
    color += texture(tex_arr,vec3(ex_Tex.xy,slice));
  else if (need_coord)
    color += ex_Tex;
    
  if (need_color)
  {
    
    fragColor = vec4(main_color,1);
  }
  else
  {
    fragColor = color;
  }
}
