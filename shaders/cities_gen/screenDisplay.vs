#version 330 core

layout (location = 0) in vec2 position;

// If you want f. e. render 1x1 pic in 2x1 screen, use scaleRatio = (2, 1)
uniform vec2 scaleRatio;

out vec2 TexCoord; 

void main()
{
    gl_Position = vec4(position, 0, 1); //Квадрат [-1;1] x [-1;1]
    TexCoord = (position * scaleRatio) * 0.5 + 0.5;
}