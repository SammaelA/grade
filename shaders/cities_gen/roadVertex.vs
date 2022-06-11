#version 330 core

layout (location = 0) in vec3 _position;
layout (location = 1) in vec3 _normal;
layout (location = 2) in vec2 _asphalt_tex_coords;
//STRAIGHT SECTION: [0, 1] -> marking texCoord
//CROSSROADS: (1, 2) -> crossroads texcoords clamp
//ROUND TURNS: [2, 3] -> round turns markings, 2 for first segment, 3 for second segment
layout (location = 3) in float _marking_tex_coords;

//CROSSROADS: side points of road triangle, meaning centers
//ROUND TURNS: center position, small radiuses
layout (location = 4) in vec4 _anchors_coords;

//CROSSROADS: handle , then point of triangle, meaning cross and sec_edge
//ROUND TURNS: big radiuses
layout (location = 5) in vec4 _handles_coords; 

uniform mat4 transform_matrix;
uniform mat4 view_proj_matrix;
// uniform mat4 sun_view_proj_matrix;

out vec3 frag_global_position;
out vec2 frag_asphalt_tex_coord;
out float frag_marking_tex_coord;
flat out vec4 frag_anchors_coords;
flat out vec4 frag_handles_coords;
out vec3 frag_normal;
// out vec4 frag_sun_coordinates;


void main()
{ 
    vec4 global_position = transform_matrix * vec4(_position, 1.0f);
    gl_Position = view_proj_matrix * global_position;

    frag_global_position = global_position.xyz;
    frag_normal = normalize((transform_matrix * vec4(_normal, 0.0f)).xyz);
    frag_asphalt_tex_coord = _asphalt_tex_coords;
    frag_marking_tex_coord = _marking_tex_coords;

    frag_anchors_coords = _anchors_coords;
    frag_handles_coords = _handles_coords;

    // frag_position = global_position.xyz;
    
    // frag_tex_coord = TexCoord;
    // frag_sun_coordinates = sun_view_proj_matrix * global_position;

    // if (_marking_tex_coords < 1.99 || _marking_tex_coords > 3.01)
    // {
    //     frag_marking_tex_coord = 1;
    // }
    // else
    // {
    //     frag_marking_tex_coord = 2.5;
    // }
}