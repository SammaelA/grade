#version 330 core

layout (location = 0) in vec3 position;

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;

void main()
{ 
    vec4 global_position = transform_matrix * vec4(position, 1.0f);
    gl_Position = view_proj_matrix * global_position;
}