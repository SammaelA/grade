#version 330

in vec2 ex_Tex;
in vec3 ex_Normal;
in float projection_error_packed;
out vec4 fragColor;

uniform sampler2D tex;
uniform int state;
uniform vec4 fixed_color;
void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (fragColor.a<0.05)
    discard;
  if (state == 1)
    fragColor = vec4(ex_Normal,projection_error_packed);
  else if (state == -1)
    fragColor.xyz = fixed_color.xyz;
  else if (state == -2)
    fragColor.xyz = vec3(fixed_color.x, fixed_color.y,  gl_FragCoord.z / gl_FragCoord.w);
}
