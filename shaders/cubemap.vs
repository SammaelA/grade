#version 330 core
layout (location = 0) in vec3 in_Position;

out vec3 TexCoords;

uniform mat4 projection;
uniform mat4 view;

void main()
{
    TexCoords = in_Position;
    gl_Position = projection * view * vec4(100*in_Position, 1.0);
}  