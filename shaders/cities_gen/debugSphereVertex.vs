#version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 TexCoord;

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;
uniform vec3 sphere_center;

out vec3 frag_normal;

void main()
{ 
    vec4 global_position = transform_matrix * vec4(position, 1.0f);
    gl_Position = view_proj_matrix * global_position;
    frag_normal = normalize(global_position.xyz - sphere_center);
}