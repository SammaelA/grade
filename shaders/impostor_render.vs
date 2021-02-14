#version 430

struct SliceVertexData 
{
    vec4 position;
};
layout(std140, binding=2) buffer Slices 
{
    readonly SliceVertexData sliceVertexes[];
};


in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projectionCamera;
uniform mat4 rot;
uniform int slice_n;
uniform int slice_offset;
uniform int slice_verts;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;


void main(void) {
	//Fragment Position in Model Space
	vec4 pos = sliceVertexes[slice_offset + slice_n*4 + gl_VertexID].position;
	pos = rot * pos;
	ex_FragPos = pos.xyz;
	ex_Normal = in_Normal;	//Pass Normal

	//Fragment in Screen Space
	gl_Position = projectionCamera * vec4(ex_FragPos, 1.0f);

	//Color from Normal Vector
	int sl_x = slice_n % 3;
	int sl_y = slice_n / 3;

	ex_Tex = in_Tex.xyz;

	ex_Tex.x = 0.333*(ex_Tex.x + sl_x);
	ex_Tex.y = 0.333*(ex_Tex.y + sl_y);
}
