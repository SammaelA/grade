#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 model;
uniform mat4 projectionCamera;
uniform mat4 lightSpaceMatrix;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec4 FragPosLightSpace;

void main(void) {
	//Fragment Position in Model Space
	ex_FragPos = (model * vec4(in_Position, 1.0f)).xyz;
	ex_Normal = (transpose(inverse(model))*vec4(in_Normal,0)).xyz;

	//Fragment in Screen Space
	gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);

	//Color from Normal Vector
	ex_Tex = in_Tex.xyz;
	FragPosLightSpace = lightSpaceMatrix * vec4(ex_FragPos, 1.0);
}
