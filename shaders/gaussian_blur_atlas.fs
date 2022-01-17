#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
uniform int pass;
uniform vec2 tex_size_inv;

void main(void) 
{
    vec4 original = texture(tex, vec3(ex_Tex, layer));
    float weights[15] = float[15](0.1061154, 0.1028506, 0.1028506, 0.09364651, 
    0.09364651, 0.0801001, 0.0801001, 0.06436224, 0.06436224, 
    0.04858317, 0.04858317, 0.03445063, 0.03445063, 0.02294906, 
    0.02294906);
    float offsets[15] = float[15](0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7);
    float blur_step = 0.5;
    vec4 res = vec4(0, 0, 0, 1);
    if (pass == 0)
    {
        for (int i = 0; i < 15; i++)
            res += texture(tex, vec3(ex_Tex, layer) + vec3(blur_step*offsets[i]*tex_size_inv.x, 0, 0)) * weights[i];
    }
    else
    {
        for (int i = 0; i < 15; i++)
            res += texture(tex, vec3(ex_Tex, layer) + vec3(0, blur_step*offsets[i]*tex_size_inv.y, 0)) * weights[i];
    }
    fragColor = vec4(res.xy,original.zw);
    //fragColor = vec4(ex_Tex, layer, 1);
}
