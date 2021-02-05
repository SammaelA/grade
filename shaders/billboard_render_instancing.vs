#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;
in vec4 in_Center_par;
in vec4 in_Center_self;
in mat4 in_Model;

uniform mat4 projectionCamera;
uniform vec2 LOD_dist_min_max;
uniform vec3 camera_pos;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;


void main(void) {
	//Fragment Position in Model Space
	ex_FragPos = (in_Model * vec4(in_Position, 1.0f)).xyz;

	if (length(in_Center_par.xyz - camera_pos) > LOD_dist_min_max.y)
	{
		gl_Position = vec4(-1000,-1000,1,1);
	}
	else if (length(in_Center_self.xyz - camera_pos) < LOD_dist_min_max.x)
	{
		gl_Position = vec4(-1000,-1000,1,1);
	}
	else
	{
		ex_Normal = in_Normal;	//Pass Normal

		//Fragment in Screen Space
		gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);

		//Color from Normal Vector
		ex_Tex = in_Tex.xyz;
	}
}