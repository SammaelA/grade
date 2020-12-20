#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 model;
uniform mat4 projectionCamera;

out vec2 ex_Tex;
out vec3 ex_Normal;


void main(void) {
	//Fragment Position in Model Space
	ex_Normal = in_Normal;	//Pass Normal

	//Fragment in Screen Space
	gl_Position = projectionCamera * model * vec4(in_Position, 1.0f);

	//Color from Normal Vector
	ex_Tex = in_Tex.xy;
}
