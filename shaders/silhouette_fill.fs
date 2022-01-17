#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
uniform int radius;
uniform int dir_threshold;
uniform vec2 tex_size_inv;
uniform vec2 slice_size;
uniform float threshold;

ivec2 ray(float step_x, float step_y, vec4 borders)
{
    ivec2 r = ivec2(0,0);
    for (int i=1;i<=radius;i++)
    {
        if (-i*step_x < borders.x && i*step_x < borders.z && -i*step_y < borders.y && i*step_y < borders.w)
        {
            vec2 pixel = texelFetch(tex, ivec3(ex_Tex/tex_size_inv + vec2(i*step_x, i*step_y),layer), 0).xy;
            r = max(r, ivec2(pixel.x > threshold, pixel.y > threshold));
        }
    }
    return clamp(r, ivec2(0,0), ivec2(1,1));
}
void main(void) 
{
    vec2 pixel_coords = ex_Tex/tex_size_inv;
    ivec2 slice_n = ivec2(pixel_coords/slice_size);
    vec4 borders = vec4(pixel_coords - vec2(slice_n)*slice_size, vec2(slice_n.x+1, slice_n.y+1)*slice_size - pixel_coords);
    vec4 original = texture(tex, vec3(ex_Tex, layer));
    if (original.x > threshold || original.y > threshold)
    {
        fragColor = original;
    }
    else
    {
        ivec2 r = ivec2(0,0);
        r += ray(-1,-1, borders);
        r += ray(-1,0, borders);
        r += ray(-1,1, borders);
        r += ray(0,-1, borders);
        r += ray(0,1, borders);
        r += ray(1,-1, borders);
        r += ray(1,0, borders);
        r += ray(1,1, borders);

        fragColor = vec4(0,0,0,0);
        fragColor.x = r.x > dir_threshold ? 1 : 0;
        fragColor.y = r.y > dir_threshold ? 1 : 0;
        fragColor.w = fragColor.x + fragColor.y > 0 ? 1 : 0;
    }
}