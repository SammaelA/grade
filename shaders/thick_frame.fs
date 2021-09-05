#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;
uniform float thickness;
uniform vec4 color;
void main(void) 
{
    if (1-abs(2*ex_Tex.x-1) < thickness || 1-abs(2*ex_Tex.y-1) < thickness)
        fragColor = vec4(color.xyz,1);
    else
        discard;
}