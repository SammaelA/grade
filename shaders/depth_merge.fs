#version 330

in vec2 ex_Tex;
uniform sampler2D trans;
uniform sampler2D depth;

out vec4 fragColor;

void main(void)
{
    vec4 tr = texture(trans,ex_Tex);
    vec4 d = texture(depth,ex_Tex);
    fragColor = vec4(1,0,0,1);
    if (d.w > 0.95)
      fragColor.xy = d.xy;//z and z*z
    if (tr.w > 0.05)
      fragColor.zw = tr.wx;//z and a
    //fragColor = vec4(d,tr > 0 ? 0.5*(2 - 1/(tr+1)) : 0);
}