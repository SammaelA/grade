#version 460


in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projection;
uniform mat4 view;

out vec3 ex_Tex;
out vec2 ex_Pos;
void main(void) 
{
    gl_Position = projection * view * vec4(in_Position, 1.0f);
    ex_Pos = 0.5*(gl_Position.xy+1);
    ex_Tex.xy = in_Tex.xy;
    ex_Tex.z = 0;
}
