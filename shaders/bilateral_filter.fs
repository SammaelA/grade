#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
out vec4 fragColor;
uniform float sigma_d;
uniform float sigma_r;
void main(void) 
{
  ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
  vec4 center = texelFetch(tex, pixel_coords + ivec2(0, 0), 0);
  int r = int(3*sigma_d);
  float total_w = 0;
  vec4 total_color = vec4(0,0,0,0);
  for (int i=max(pixel_coords.x - r, 0); i <= min(pixel_coords.x + r, int(tex_size.x) - 1); i++)
  {
    for (int j=max(pixel_coords.y - r, 0); j <= min(pixel_coords.y + r, int(tex_size.y) - 1); j++)
    {
      vec4 c = texelFetch(tex, ivec2(i, j), 0);
      float w0 = ((i - pixel_coords.x)*(i - pixel_coords.x) + (j - pixel_coords.y)*(j - pixel_coords.y))/(2.0f*sigma_d*sigma_d);
      float diff = dot(center - c, vec4(0.299, 0.587, 0.114, 1));//set bigger weight to prevent mixing with different alpha values
      float w1 = diff*diff/(2*sigma_r*sigma_r);
      float w = exp(-(w0+w1));
      total_color += w*c;
      total_w += w;
    }
  }
  fragColor = total_color/total_w;
}