#version 330

out vec4 result;
uniform vec3 points_color;

void main(void) 
{
  result = vec4(points_color,1);
}
