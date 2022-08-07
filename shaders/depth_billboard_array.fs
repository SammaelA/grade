#version 330 core
in vec3 ex_Tex;
out vec4 fragColor;
uniform sampler2DArray tex;
uniform float opaqueness = 1e4;
void main() 
{
  float alpha = texture(tex, ex_Tex).a;
    if (alpha < 0.05)
        discard;
    float z = (gl_FragCoord.z);
    fragColor = vec4 ( z, z*z, z, alpha*opaqueness);
}
