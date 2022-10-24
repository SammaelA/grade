#version 460

layout(std140, binding=12) readonly buffer _InstanceData
{
    mat4 instances[];
};

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

out vec3 ex_Tex;

void main(void) 
{
    ex_Tex.xy = in_Tex.xy;
    ex_Tex.z = 0;
}
