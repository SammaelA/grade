#version 460

layout(std140, binding=11) readonly buffer _InstanceData
{
    mat4 instances[];
};

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projection;
uniform mat4 view;
uniform int inst_buf_offset;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec4 ex_FragPosView;


void main(void) 
{ 
    int instance_n = inst_buf_offset + int(gl_InstanceID);
    mat4 inst_mat = instances[instance_n];
    ex_Normal = normalize(ex_Normal = (transpose(inverse(inst_mat))*vec4(in_Normal,0)).xyz);
    ex_FragPos = (inst_mat*vec4(in_Position, 1.0f)).xyz;
    ex_FragPosView = view * vec4(ex_FragPos, 1.0f);
    gl_Position = projection * ex_FragPosView;
    ex_Tex.xy = in_Tex.xy;
    ex_Tex.z = 0;
}
