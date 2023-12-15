#version 460

in vec3 in_Position;

uniform mat4 viewProj;
uniform float points_base_size;

void main(void) 
{
  gl_Position = viewProj * vec4(in_Position, 1.0f);
  gl_PointSize = points_base_size*gl_Position.z; 
}
