#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
uniform float blur_step = 0.5;
uniform int pass;
uniform vec2 tex_size_inv;
uniform vec2 slice_size;
void main(void) 
{
    vec4 original = texture(tex, vec3(ex_Tex, layer));
    /*
    float weights[15] = float[15](0.1061154, 0.1028506, 0.1028506, 0.09364651, 
    0.09364651, 0.0801001, 0.0801001, 0.06436224, 0.06436224, 
    0.04858317, 0.04858317, 0.03445063, 0.03445063, 0.02294906, 
    0.02294906);*/
    float weights[15] = float[15](0.02294906, 0.03445063, 0.04858317, 0.06436224, 0.0801001, 0.09364651, 0.1028506,
                                  0.1061154,  
                                  0.1028506, 0.09364651,  0.0801001,  0.06436224, 0.04858317, 0.03445063, 0.02294906);
    float offsets[15] = float[15](-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7);
    //float blur_step = 0.5;
    vec4 res = vec4(0, 0, 0, 1);
    vec2 pixel_coords = ex_Tex/tex_size_inv;
    ivec2 slice_n = ivec2(pixel_coords/slice_size);
    vec4 borders = vec4(pixel_coords - vec2(slice_n)*slice_size, vec2(slice_n.x+1, slice_n.y+1)*slice_size - pixel_coords);
    float w = 1e-6;
    if (pass == 0)
    {
        for (int i = max(-7, int(-borders.x/blur_step)); i <= min(7, int(borders.z/blur_step)); i++)
        {
            res += texture(tex, vec3(ex_Tex, layer) + vec3(blur_step*offsets[i+7]*tex_size_inv.x, 0, 0)) * weights[i+7];
            w += weights[i+7];
        }
    }
    else
    {
        for (int i = max(-7, int(-borders.y/blur_step)); i <= min(7, int(borders.w/blur_step)); i++)
        {
            res += texture(tex, vec3(ex_Tex, layer) + vec3(0, blur_step*offsets[i+7]*tex_size_inv.y, 0)) * weights[i+7];
            w += weights[i+7];
        }
    }
    res.xy /= w;
    fragColor = vec4(res.xy,original.zw);
    //fragColor = vec4(ex_Tex, layer, 1);
}
