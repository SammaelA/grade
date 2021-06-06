#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 model;
uniform mat4 projectionCamera;
uniform float projection_zero;
out vec2 ex_Tex;
out vec3 ex_Normal;
out float projection_error_packed;

void main(void) {

	mat4 tr = projectionCamera * model;
	vec4 pos = model * vec4(in_Position, 1.0f);
	gl_Position = projectionCamera * pos;
	ex_Normal = 0.5*(normalize((transpose(inverse(tr))*vec4(in_Normal,0)).xyz) + 1);
	ex_Tex = in_Tex.xy;
	projection_error_packed = 0.5*(0.01*(pos.z - projection_zero) + 1);
}
