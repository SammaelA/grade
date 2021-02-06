#version 330

in vec3 in_Position;
in vec4 in_Tex;

out vec3 ex_Tex;
uniform float layer;

void main(void) 
{
	gl_Position = vec4(2*in_Position - vec3(1,1,0), 1.0f);
	ex_Tex = vec3(in_Position.xy,layer);
}