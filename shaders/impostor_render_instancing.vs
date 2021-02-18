#version 430

struct SliceVertexData 
{
    vec4 position;
	vec4 tcs;
};
layout(std140, binding=2) buffer Slices 
{
    SliceVertexData sliceVertexes[];
};
struct LOD_info
{
    uint offset;
    float min_dist;
    float max_dist;
    float pad;
};
struct InstanceData
{
    vec4 center_self;
    vec4 center_par;
    mat4 projection_camera;
};
struct IndexData
{
    uint index;
    //uint pad;
    float mn;
    float mx;
};
layout(std140, binding=1) buffer _models_intervals
{
    readonly uvec2 models_intervals[];
};
layout(std140, binding=3) buffer _instances
{
    readonly InstanceData instances[];
};
layout(std140, binding=4) buffer _instance_indexes//same size as instances
{
    readonly IndexData instance_indexes[];
};


in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projectionCamera;
uniform vec3 imp_center;
uniform vec3 camera_pos;
uniform int slice_count;
uniform int slice_offset;
uniform uint id;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec2 a_mult;

#define PI 3.1415926535897932384626433832795

mat4 translate(vec3 d)
{
	return mat4(1, 0, 0, d.x,
	            0, 1, 0, d.y,
	            0, 0, 1, d.z,
	            0, 0, 0, 1);
}
mat4 rotate(vec3 axis, float angle) 
{
  axis = normalize(axis);
  float s = sin(angle);
  float c = cos(angle);
  float oc = 1.0 - c;

  return mat4(
		oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
    oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
    oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
		0.0,                                0.0,                                0.0,                                1.0
	);
}
void main(void) {
	//Fragment Position in Model Space
	uint offset = models_intervals[id].x + gl_InstanceID;
	mat4 inst_mat = instances[instance_indexes[offset].index].projection_camera;
	a_mult = vec2(instance_indexes[offset].mn, instance_indexes[offset].mx);
	vec3 center = (inst_mat * vec4(imp_center,1.0f)).xyz;
    vec3 viewdir = center - camera_pos;
    viewdir.y = 0;
    viewdir = normalize(viewdir);
    float phi = acos(dot(vec3(0,0,-1),viewdir));
    if (viewdir.x > 0)
        phi = 2*PI - phi;
    int slice_n = int(round(slice_count*(phi/(2.0f*PI)))) % slice_count; 
    phi = - phi + (2.0f*PI*float(slice_n))/float(slice_count);
	mat4 rot = rotate(vec3(0,1,0),phi);
	slice_n += 1;

	vec4 pos = vec4(sliceVertexes[slice_offset + slice_n*4 + gl_VertexID].position.xyz,1);
	pos = vec4(center,0) + rot * (inst_mat * pos - vec4(center,0));
	ex_FragPos = pos.xyz;
	ex_Normal = in_Normal;	//Pass Normal

	//Fragment in Screen Space
	gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);
	ex_Tex = sliceVertexes[slice_offset + slice_n*4 + gl_VertexID].tcs.xyz;
}
