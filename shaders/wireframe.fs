#version 330

out vec4 fragColor;
uniform vec3 wireframe_color;
void main(void) 
{
  fragColor = vec4(wireframe_color,1);
}
