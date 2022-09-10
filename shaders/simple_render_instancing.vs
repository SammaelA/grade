#version 460

struct TypeData
{
    uint offset;
    uint pad1,pad2,pad3;
};
layout(std140, binding=1) readonly buffer _TypeData
{
    TypeData typeData[];
};
struct InstanceData
{
    vec4 center_self;
    vec4 center_par;
    mat4 projection_camera;
};
layout(std140, binding=3) readonly buffer _InstanceData
{
    InstanceData instances[];
};
struct ModelData
{
    uint LOD;
    uint type;
    uint vertexes;
    uint first_index;
    uvec2 interval;
    uint culling;
    uint pad;
    vec4 x_s;
    vec4 y_s;
    vec4 z_s;
};
layout(std140, binding=4) readonly buffer _ModelData
{
    ModelData modelData[];
};
struct currentInstancesData
{
    uint index;
    uint pad;
    float mn;
    float mx;
};
layout(std140, binding=5) readonly buffer _curInsts
{
    currentInstancesData curInsts[];
};

layout(std140, binding=6) readonly buffer _curModels
{
    uvec4 curModels[];
};

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projection;
uniform mat4 view;
uniform vec2 LOD_dist_min_max;
uniform vec3 camera_pos;
uniform uint type_id;
uniform mat4 lightSpaceMatrix;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec4 ex_FragPosView;
out vec2 a_mult;
flat out uint model_id;
out vec4 FragPosLightSpace;
out vec3 colorMult;
void main(void) 
{
    uint id = typeData[type_id].offset + gl_DrawID;
	uint offset = curModels[id].x + gl_InstanceID;
	mat4 inst_mat = instances[curInsts[offset].index].projection_camera;
	ex_FragPos = (inst_mat * vec4(in_Position, 1.0f)).xyz;
	a_mult = vec2(curInsts[offset].mn, curInsts[offset].mx);
	
	ex_Normal = normalize(ex_Normal = (transpose(inverse(inst_mat))*vec4(in_Normal,0)).xyz);
	ex_FragPosView = view * vec4(ex_FragPos, 1.0f);
    gl_Position = projection * ex_FragPosView;
	ex_Tex = in_Tex.xyz;
    vec4 cp = instances[curInsts[offset].index].center_par;
    vec4 cs = instances[curInsts[offset].index].center_self;
    ex_Tex.z = cp.w;
    model_id = curModels[id].y;
	FragPosLightSpace = lightSpaceMatrix * vec4(ex_FragPos, 1.0);

    vec3 rnd1 = 2*fract(cp.xyz) - 1;
    rnd1.z = 0;
    vec3 rnd2 = 2*fract(cs.xyz) - 1;
    rnd2.z = 0;
    float rnd3 = 2*fract(10*length(cp.xyz)) - 1;
    float rnd4 = 2*fract(10*length(cs.xyz)) - 1;
    colorMult = (1 + 0.1*rnd1 + 0.05*rnd2)*(1 + 0.2*rnd3)*(1 + 0.1*rnd4);
}