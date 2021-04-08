#version 460

struct SliceVertexData 
{
    vec4 position;
	vec4 tcs;
};
layout(std140, binding=2) readonly buffer Slices 
{
    SliceVertexData sliceVertexes[];
};
struct ImpostorData
{
    int slice_offset;
    int slice_verts;
    int slice_count;
    int pad1;
    vec4 imp_center;
};
layout(std140, binding=9) readonly buffer _ImpostorData 
{
    ImpostorData impostorData[];
};
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

uniform mat4 projectionCamera;
//uniform vec3 imp_center;
uniform vec3 camera_pos;
//uniform int slice_count;
//uniform int slice_offset;
uniform float delta;
//uniform float angle_step;
uniform vec2 hor_vert_transition_thr;
uniform uint type_id;
uniform int vertex_id_offset;

out vec3 tc_a;
out vec3 tc_b;
out vec3 tc_t;
out vec3 q_abt;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec2 a_mult;
flat out uint model_id;
out mat4 rot_m;
#define PI 3.1415926535897932384626433832795
uniform mat4 lightSpaceMatrix;
out vec4 FragPosLightSpace;
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
	uint id = typeData[type_id].offset + gl_DrawID;
	uint offset = curModels[id].x + gl_InstanceID;
	mat4 inst_mat = instances[curInsts[offset].index].projection_camera;
	a_mult = vec2(curInsts[offset].mn, curInsts[offset].mx);
    
    uint model_rel_n = curModels[id].y - typeData[type_id].offset;
    int vertex_id = (gl_VertexID - vertex_id_offset) % 4;
    int slice_offset = impostorData[model_rel_n].slice_offset;
    int slice_verts = impostorData[model_rel_n].slice_verts;
    int slice_count = impostorData[model_rel_n].slice_count;
    vec3 imp_center = impostorData[model_rel_n].imp_center.xyz;
    float angle_step = (2*PI)/float(slice_count);
    
	vec3 center = (inst_mat * vec4(imp_center,1.0f)).xyz;
    vec3 viewdir = normalize(center - camera_pos);
    mat3 MVI = mat3(transpose(inverse(inst_mat)));
    vec3 front = normalize(MVI * vec3(0,0,1));
    vec3 plane = normalize(MVI * vec3(1,0,0));
    vec3 up = normalize(MVI * vec3(0,1,0));
    float psi = -asin(clamp(dot(up, viewdir),-1,0));
    viewdir.y = 0;
    viewdir = normalize(viewdir);
    float phi = acos(-dot(front,viewdir));
    if (dot(viewdir,plane) > 0)
        phi = 2*PI - phi;
    float a_phi = phi;
    
    int slice_n = int(floor(slice_count*(phi/(2.0f*PI)))) % slice_count;

    phi = - phi + (2.0f*PI*float(slice_n))/float(slice_count);
    
    int second_slice_n = (slice_n - int(sign(phi)) + slice_count) % slice_count; 

    slice_n += 1;
    second_slice_n += 1;

	mat4 rot = rotate(up,phi);
    rot_m = inst_mat;
    mat4 rot_top = rotate(cross(up, viewdir),psi);
	
	vec4 pos = vec4(sliceVertexes[slice_offset + slice_n*4 + vertex_id].position.xyz,1);

    //ex_FragPos.z = float(vertex_id + 0.0);
    //ex_FragPos.y = float(slice_offset + slice_n*4 + vertex_id + 0.0);
	pos = vec4(center,0) + rot * (inst_mat * pos - vec4(center,0));

    float s_phi =  sign(phi)*(angle_step - abs(phi));
    mat4 rot_s = rotate(up, - s_phi);
    vec4 pos_s = vec4(sliceVertexes[slice_offset + second_slice_n*4 + vertex_id].position.xyz,1);
	pos_s = vec4(center,0) + rot_s * (inst_mat * pos_s - vec4(center,0));
    float hth = clamp(hor_vert_transition_thr.x,0.0f,0.5f);
    float rel_delta = delta*angle_step;
    q_abt.y = smoothstep(hth, 1 - hth, 0.5*(1 + (abs(phi) - 0.5*angle_step)/rel_delta));
    q_abt.x = 1 - q_abt.y;
    pos = q_abt.x * pos + q_abt.y * pos_s;

    pos = vec4(center,0) + rot_top * (pos - vec4(center,0));
    mat4 rot_top_s = rotate(cross(up, viewdir),-0.5*PI + psi);
    mat4 rot_top_s2 = rotate(up, -a_phi);
    vec4 top_pos = inst_mat * vec4(sliceVertexes[slice_offset + vertex_id].position.xyz,1);
    top_pos = vec4(center,0) + rot_top_s2 * (top_pos - vec4(center,0));
    
    float vth = clamp(hor_vert_transition_thr.y,0.0f,0.5f);
    q_abt.z = smoothstep(vth, 1-vth, abs(psi)/(0.5*PI));
    q_abt.xy *= (1 - q_abt.z);
    pos = (1 - q_abt.z)*pos + q_abt.z*top_pos;
	ex_FragPos = pos.xyz;
    ex_Normal = (transpose(inverse(inst_mat))*vec4(normalize(camera_pos - pos.xyz),0)).xyz;

	//Fragment in Screen Space
	gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);
	tc_a = sliceVertexes[slice_offset + slice_n*4 + vertex_id].tcs.xyz;
    tc_b = sliceVertexes[slice_offset + second_slice_n*4 + vertex_id].tcs.xyz;
    tc_t = sliceVertexes[slice_offset + vertex_id].tcs.xyz;

    ex_FragPos = pos_s.xyz;
    ex_FragPos.x = curInsts[offset].index;
    ex_FragPos.y = offset;
    //ex_FragPos.y = slice_n;
    //ex_FragPos.z = gl_DrawID;
    model_id = curModels[id].y;
    FragPosLightSpace = lightSpaceMatrix * vec4(ex_FragPos, 1.0);
}
