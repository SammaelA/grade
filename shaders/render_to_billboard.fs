#version 330

in vec2 ex_Tex;
in vec3 ex_Normal;
in float projection_error_packed;
out vec4 fragColor;

uniform sampler2D tex;
uniform int state;

void main(void) 
{
  fragColor = texture(tex,ex_Tex);
  if (fragColor.a<0.33)
    discard;
  if (state == 1)
    fragColor = vec4(ex_Normal,projection_error_packed);
}
