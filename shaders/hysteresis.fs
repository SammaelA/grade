#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
uniform vec2 thresholds;
out vec4 fragColor;

void main(void) 
{
    vec4 res = vec4(0, 0, 0, 1);
    ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
    float val = texelFetch(tex, pixel_coords, 0).x;
    float lt = 0.25*thresholds.x;
    float ht = 0.25*thresholds.y;
    if (pixel_coords.x > 0 && pixel_coords.y > 0 && pixel_coords.x < int(tex_size.x-1) && pixel_coords.y < int(tex_size.y-1))
    {
      if (val > ht)
        res.x = 1;
      else if (val < lt)
        res.x = 0;
      else
      {
        #define GET(a,b) texelFetch(tex, pixel_coords + ivec2(a, b), 0).x
        float mx = max(max(max(GET(-1,-1), GET(-1, 0)), max(GET(-1,-1), GET(0, -1))),
                      max(max(GET(0,-1), GET(1, -1)), max(GET(1,0), GET(1, 1))));
        if (mx > ht)
          res.x = 1;
        else
          res.x = 0;
      }
    }
    fragColor = res;
}