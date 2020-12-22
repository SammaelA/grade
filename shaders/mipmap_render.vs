#version 330

in vec3 in_Position;
in vec4 in_Tex;

out vec2 ex_Tex;


void main(void) 
{
	gl_Position = vec4(2*in_Position - vec3(1,1,0), 1.0f);
	ex_Tex = in_Position.xy;
}