#version 330

in vec2 ex_Tex;
uniform sampler2D trans;
uniform sampler2D depth;

out vec4 fragColor;

void main(void)
{
    float tr = texture(trans,ex_Tex).a;
    vec3 d = texture(depth,ex_Tex).xyz;
    fragColor = vec4(d,tr > 0 ? 0.5*(2 - 1/(tr+1)) : 0);
}