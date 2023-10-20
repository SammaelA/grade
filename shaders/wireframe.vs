#version 460


in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 viewProj;
uniform mat4 view;

void main(void) 
{
    gl_Position = viewProj * vec4(in_Position, 1.0f);
    gl_Position.y *= -1;
}