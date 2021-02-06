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

out vec2 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec2 a_mult;

void main(void) {
	//Fragment Position in Model Space
	ex_FragPos = (in_Model * vec4(in_Position, 1.0f)).xyz;
	float mx_dist = LOD_dist_min_max.y - length(in_Center_par.xyz - camera_pos);
	float mn_dist = length(in_Center_self.xyz - camera_pos) - LOD_dist_min_max.x;
	float trans = 6;
	if (mx_dist < -trans)
	{
		gl_Position = vec4(-1000,-1000,1,1);
	}
	else if (mn_dist < -trans)
	{
		gl_Position = vec4(-1000,-1000,1,1);
	}
	else
	{
		if (mn_dist < trans)
		{
			a_mult = vec2(0.5*(mn_dist/trans + 1), 1);
		}
		else if (mx_dist < trans)
		{
			a_mult = vec2(0.5*(mx_dist/trans + 1), -1);
		}
		else
		{
			a_mult = vec2(10,0);
		}
		ex_Normal = in_Normal;	//Pass Normal

		//Fragment in Screen Space
		gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);

		//Color from Normal Vector
		ex_Tex = in_Tex.xy;
	}
}