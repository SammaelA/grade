#version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 TexCoord;

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;
uniform mat4 sun_view_proj_matrix;

out vec3 frag_position;
out vec3 frag_normal;
out vec2 frag_tex_coord;
out vec4 frag_sun_coordinates;

void main()
{ 
    vec4 global_position = transform_matrix * vec4(position, 1.0f);
    gl_Position = view_proj_matrix * global_position;

    frag_position = global_position.xyz;
    frag_normal = normalize((transform_matrix * vec4(normal, 0.0f)).xyz);
    frag_tex_coord = TexCoord;
    frag_sun_coordinates = sun_view_proj_matrix * global_position;
}