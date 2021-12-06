#version 330 core
in vec3 ex_Tex;
out vec4 fragColor;

void main() 
{
    float z = (gl_FragCoord.z);
    fragColor = vec4 ( z, z*z, z, 1);
}
