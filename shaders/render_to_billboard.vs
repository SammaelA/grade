#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 model;
uniform mat4 projectionCamera;

out vec2 ex_Tex;
out vec3 ex_Normal;


void main(void) {

	mat4 tr = projectionCamera * model;
	gl_Position = tr * vec4(in_Position, 1.0f);
	ex_Normal = (transpose(inverse(tr))*vec4(in_Normal,0)).xyz;
	ex_Tex = in_Tex.xy;
}
