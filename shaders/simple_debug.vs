#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec4 ex_FragPosView;


void main(void) 
{ 
  ex_Normal = (transpose(inverse(model))*vec4(in_Normal,0)).xyz;
  ex_FragPos = (model*vec4(in_Position, 1.0f)).xyz;
  ex_FragPosView = view * vec4(ex_FragPos, 1.0f);
  gl_Position = projection * ex_FragPosView;
  ex_Tex = in_Tex.xyz;
}
