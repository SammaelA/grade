#version 430

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
    readonly uvec4 models_intervals[];
};
layout(std140, binding=3) buffer _instances
{
    readonly InstanceData instances[];
};
layout(std140, binding=4) buffer _instance_indexes//same size as instances
{
    readonly IndexData instance_indexes[];
};
layout(std140, binding=5) buffer _models_instance_count//same size as models_intervals
{
    readonly uint models_instance_count[];
};
in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projectionCamera;
uniform vec2 LOD_dist_min_max;
uniform vec3 camera_pos;
uniform uint id;

out vec2 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec2 a_mult;

void main(void) 
{
	uint offset = models_intervals[id].x + gl_InstanceID;
	mat4 inst_mat = instances[instance_indexes[offset].index].projection_camera;
	ex_FragPos = (inst_mat * vec4(in_Position, 1.0f)).xyz;
	a_mult = vec2(instance_indexes[offset].mn, instance_indexes[offset].mx);
	
	ex_Normal = in_Normal;
	gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);
	ex_Tex = in_Tex.xy;
	ex_Tex.x += 0.01*float(models_instance_count[id]);
}