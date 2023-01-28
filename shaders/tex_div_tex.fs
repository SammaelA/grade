#version 330

in vec2 ex_Tex;
uniform sampler2D tex1;
uniform sampler2D tex2;
out vec4 fragColor;

void main(void) 
{
    fragColor = texture(tex1, ex_Tex) - texture(tex2, ex_Tex);
}
