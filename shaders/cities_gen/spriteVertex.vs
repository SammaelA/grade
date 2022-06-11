#version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec2 texCoord;

uniform vec4 Transform;
out vec2 TexCoord;

void main()
{
    vec4 screenTransform = Transform * vec4(2.0);
    vec4 out_pos = vec4(screenTransform[0], screenTransform[1], 1.0, 1.0) * vec4(position.x, position.y, position.z, 1.0);
    out_pos = out_pos + vec4(screenTransform[2], screenTransform[3], 0.0, 0.0);
    out_pos = out_pos + vec4(-1, -1, 0.0, 0.0);
    gl_Position = out_pos;
    TexCoord = vec2(texCoord.x, 1 - texCoord.y); 
}