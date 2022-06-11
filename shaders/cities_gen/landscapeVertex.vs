#version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in float _slope_ratio;
layout (location = 3) in float _water_ratio;
// layout (location = 4) in vec2 texture_coords;

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;
uniform mat4 sun_view_proj_matrix;

//CLIPPING UNIFORMS
uniform vec4 clipping_plane;

// 0 - default
// 1 - clipping
uniform int render_mode;

out vec3 frag_position;
out vec3 frag_normal;
// out vec2 frag_tex_coord;
out vec4 frag_sun_coordinates;
out float slope_ratio;
out float water_ratio;

void main()
{ 
    vec4 global_position = transform_matrix * vec4(position, 1.0f);
    gl_Position = view_proj_matrix * global_position;

    frag_position = global_position.xyz;
    frag_normal = normal;//normalize((transform_matrix * vec4(normal, 0.0f)).xyz);
    // frag_tex_coord = texture_coords;
    frag_sun_coordinates = sun_view_proj_matrix * global_position;

    slope_ratio = _slope_ratio;
    water_ratio = _water_ratio;

    // fragCol = col;
    // vec4 global_position = transform_matrix * vec4(position, 1.0f);
    // gl_Position = view_proj_matrix * global_position;

    //CLIPPING
    if (render_mode == 1)
        gl_ClipDistance[0] = dot(global_position, clipping_plane);
}