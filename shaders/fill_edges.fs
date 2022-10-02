#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
uniform int search_radius;
out vec4 fragColor;

void main(void) 
{
  //if the pixel is between two edge pixels, it should be inside a body, that these edges are borders of
    vec4 res = vec4(0, 0, 0, 1);
    ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
    float val = texelFetch(tex, pixel_coords, 0).x;
    if (val == 0)
    {
      #define GET(a,b) texelFetch(tex, pixel_coords + ivec2(a, b), 0).x

      float x_steps[4] = float[4](1,1,0,-1);
      float y_steps[4] = float[4](0,1,1,1);
      int possible_dir = -1;
      for (int i = 0; i < 4; i++)
      {
        int dist = -1;
        for (int j=1;j<search_radius;j++)
        {
          float v = GET(j*x_steps[i], j*y_steps[j]);
          if (v > 0)
          {
            dist = j;
            break;
          }
        }
        if (dist >= 0)
        {
          dist = -1;
          for (int j=1;j<search_radius;j++)
          {
            float v = GET(-j*x_steps[i], -j*y_steps[j]);
            if (v > 0)
            {
              dist = j;
              break;
            }
          }
        }
        if (dist >= 0)
        {
          possible_dir = i;
          break;
        }
      }
      if (possible_dir >= 0)
      {
        res.x = 1;
      }
    }
    else
    {
      res.x = 1;
    }
    fragColor = res;
}