#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;
in mat4 in_Model;

uniform mat4 projectionCamera;

out vec2 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;


void main(void) {
	//Fragment Position in Model Space
	ex_FragPos = (in_Model * vec4(in_Position, 1.0f)).xyz;
	ex_Normal = in_Normal;	//Pass Normal

	//Fragment in Screen Space
	gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);

	//Color from Normal Vector
	ex_Tex = in_Tex.xy;
}
