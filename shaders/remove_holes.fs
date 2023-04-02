#version 330

in vec2 ex_Tex;
uniform sampler2D tex_mask;
uniform int search_radius;
uniform vec2 tex_size;
out vec4 fragColor;

void main(void) 
{
  ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
  vec3 color = texelFetch(tex_mask, pixel_coords, 0).xyz;
  if (color.x < 0.5)
  {
    vec2 tex_size_inv = vec2(1/tex_size.x, 1/tex_size.y);
    int rays = search_radius + 4;
    int intersects = 0;
    #define PI 3.14159
    float phi = 2*PI/rays;
    for (int i=0;i<rays;i++)
    {
      vec2 dir = tex_size_inv*vec2(cos(i*phi), sin(i*phi));
      for (int r=2;r<=search_radius;r++)
      {
        if (texture(tex_mask, ex_Tex + r*dir).x > 0.5)
        {
          intersects += 1;
          break;
        }
      }
    }
    if (intersects == rays)
      fragColor = vec4(1,1,1,1);
    else
      fragColor = vec4(0,0,0,1);
  }
  else
    fragColor = vec4(1,1,1,1);
}