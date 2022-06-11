#version 330 core

layout (location = 0) in vec3 _position;
layout (location = 1) in vec3 _normal;
layout (location = 2) in vec3 _tex_coords;

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;
// uniform mat4 sun_view_proj_matrix;

// out vec3 frag_global_position;
out vec3 frag_tex_coord;
// out float frag_marking_tex_coord;
// flat out vec4 frag_anchors_coords;
// flat out vec4 frag_handles_coords;
out vec3 frag_normal;
out vec3 frag_position;
// out vec4 frag_sun_coordinates;


void main()
{ 
    vec4 global_position = transform_matrix * vec4(_position, 1.0f);
    gl_Position = view_proj_matrix * global_position;

    frag_normal = normalize((transform_matrix * vec4(_normal, 0.0f)).xyz);
    frag_position = global_position.xyz;
    frag_tex_coord = _tex_coords;
    frag_tex_coord.y *= -1;
}