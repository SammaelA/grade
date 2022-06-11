#version 330 core

layout (location = 0) in vec3 position;
// layout (location = 1) in vec3 normal;
// layout (location = 2) in vec2 TexCoord;

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;

out vec3 frag_global_pos;

void main()
{ 
    vec4 globalPosition = transform_matrix * vec4(position, 1.0f);
    gl_Position = view_proj_matrix * globalPosition;
    frag_global_pos = globalPosition.xyz;
}