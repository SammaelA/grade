#version 330 core
in vec3 ex_Tex;
out vec4 fragColor;
uniform sampler2D tex;
uniform float opaqueness = 1e4;
void main() 
{
    if (texture(tex, ex_Tex.xy).a < 0.95)
        discard;
    float z = (gl_FragCoord.z);
    fragColor = vec4 ( z, z*z, z, opaqueness);
}
