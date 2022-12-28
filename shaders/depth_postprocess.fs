#version 330

in vec2 ex_Tex;
uniform sampler2D depth_tex;
uniform vec2 tex_size;
uniform float zNear;
uniform float zFar;
out vec4 fragColor;

void main(void) 
{
  ivec2 tc = ivec2(tex_size*ex_Tex);
  float depth = texelFetch(depth_tex, tc, 0).x;
  float z = depth * 2.0 - 1.0;
  fragColor = vec4(0, 0, 0, clamp(100*(1 - depth), 0, 1));
  fragColor.x = 0.1 * (2.0 * zNear * zFar) / (zFar + zNear - z * (zFar - zNear));	
}