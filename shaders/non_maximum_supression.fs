#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
out vec4 fragColor;

void main(void) 
{
    vec4 res = vec4(0, 0, 0, 1);
    ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
    vec3 val = texelFetch(tex, pixel_coords, 0).xyz;
    float v1 = 0;
    float v2 = 0;
    int angle = int(round(4*val.x));
    if (angle == 0)
    {
      v1 = texelFetch(tex, pixel_coords + ivec2(-1,0), 0).y;
      v2 = texelFetch(tex, pixel_coords + ivec2(1,0), 0).y;
    }
    else if (angle == 1)
    {
      v1 = texelFetch(tex, pixel_coords + ivec2(-1,-1), 0).y;
      v2 = texelFetch(tex, pixel_coords + ivec2(1,1), 0).y;
    }
    else if (angle == 2)
    {
      v1 = texelFetch(tex, pixel_coords + ivec2(0,-1), 0).y;
      v2 = texelFetch(tex, pixel_coords + ivec2(0,1), 0).y;
    }
    else
    {
      v1 = texelFetch(tex, pixel_coords + ivec2(-1,1), 0).y;
      v2 = texelFetch(tex, pixel_coords + ivec2(1,-1), 0).y;
    }

    if (val.y >= v1 && val.y >= v2)
      res.x = val.y;
    fragColor = res;
}