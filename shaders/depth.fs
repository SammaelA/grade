#version 330 core
out vec4 fragColor;
void main() 
{
    float z = gl_FragCoord.z;
    fragColor = vec4 ( z, z*z, 0, 1.0 );
}
