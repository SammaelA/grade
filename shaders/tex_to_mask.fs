#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
out vec4 fragColor;

void main(void) 
{
    fragColor = vec4(0, 0, 0, 0);
    if (texture(tex, ex_Tex).a > 0)
    {
        fragColor = vec4(1, 1, 1, 1);
    }
}
