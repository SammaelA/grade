#version 330

in vec3 in_Position;
in vec4 in_Tex;

out vec2 ex_Tex;
uniform vec4 tex_transform = vec4(0,0,1,1);
void main(void) 
{
	gl_Position = vec4(2*in_Position - vec3(1,1,0), 1.0f);
	ex_Tex = tex_transform.xy + tex_transform.zw*in_Position.xy;
}