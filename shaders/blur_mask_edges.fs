#version 330

in vec2 ex_Tex;
uniform sampler2D tex_mask;
uniform int search_radius;
uniform vec2 tex_size;
out vec4 fragColor;

float get_dist(ivec2 pixel_coords, int max_dist)
{
  float c = texelFetch(tex_mask, pixel_coords, 0).x;
  float sr_sq = max_dist*max_dist;
  float min_dist = sr_sq;
  for (int i=-max_dist;i<=max_dist;i++)
  {
   for (int j=-max_dist;j<=max_dist;j++)
    {
      float r_q = i*i + j*j;
      if (r_q < sr_sq)
      {
        float c1 = texelFetch(tex_mask, pixel_coords + ivec2(i,j), 0).x;
        if ((c1 > 0) != (c > 0))
        {
          min_dist = min(min_dist, r_q);
        }
      }
    }   
  }

  return sqrt(min_dist);
}

void main(void) 
{
  vec4 res = vec4(0, 0, 0, 1);
  ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
  vec3 color = texelFetch(tex_mask, pixel_coords, 0).xyz;
  res.xyz = color;
  if (search_radius < 0)
  {
    //we blur _inside_ the mask border, i.e 0000055555 -> 0000012345
    if (color.x > 0)
    {
      float dist = get_dist(pixel_coords, -search_radius);
      res.xyz = -(dist/search_radius)*color;
    }
    else
      res.xyz = color;
  }
  else if (search_radius > 0)
  {
    //we blur _outside_ the mask border, i.e 0000055555 -> 0123455555
    if (color.x > 0)
      res.xyz = color;
    else
    {
      float dist = get_dist(pixel_coords, search_radius);
      res.xyz = (dist/search_radius)*color;
    }
  }

  fragColor = res;
}