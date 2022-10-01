#version 330

in vec2 ex_Tex;
uniform sampler2D tex;
uniform vec2 tex_size;
out vec4 fragColor;

void main(void) 
{
    vec4 res = vec4(0, 0, 0, 1);
    ivec2 pixel_coords = ivec2(tex_size*ex_Tex);
    if (pixel_coords.x > 0 && pixel_coords.y > 0 && pixel_coords.x < int(tex_size.x-1) && pixel_coords.y < int(tex_size.y-1))
    {
      #define GET(a,b) dot(texelFetch(tex, pixel_coords + ivec2(a, b), 0).xyz, vec3(0.2126, 0.7152, 0.0722))
      #define PI 3.14159
      float grad_x = -GET(-1,-1) + GET(1, -1) - 2*GET(-1, 0) + 2*GET(1,0) - GET(-1,1) + GET(1,1);
      float grad_y = -GET(-1,-1) + GET(-1, 1) - 2*GET(0, -1) + 2*GET(0,1) - GET(1,-1) + GET(1,1);

      float angle = (grad_x == 0.0 ? sign(grad_y)*PI/2 : atan(grad_y, grad_x));
      if (angle < 0)
        angle = PI + angle;
      int d_angle = int(round(4*angle/PI)) % 4;
      res.x = 0.25*d_angle;
      res.y = 0.25*length(vec2(grad_x, grad_y));
    }
    fragColor = res;
}
