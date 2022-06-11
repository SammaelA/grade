#version 330 core

layout (location = 0) in vec3 position;

uniform mat4 view_proj_matrix;
uniform mat4 skybox_rotation_matrix;
uniform float zFar;

out vec3 texCoords;

void main()
{ 
    gl_Position = view_proj_matrix * vec4(position, 1);
    gl_Position.z = gl_Position.w;
    texCoords = (skybox_rotation_matrix * vec4(position, 1)).xyz;
}