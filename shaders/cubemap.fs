#version 330 core

in vec3 TexCoords;

uniform samplerCube skybox;

layout (location = 0) out vec4 fragColor;
layout (location = 1) out vec4 fragNormal;
layout (location = 2) out vec4 fragViewPos;
layout (location = 3) out vec4 fragWorldPos;

void main()
{    
    fragColor = texture(skybox, TexCoords);
    fragNormal = vec4(0,0,0,-1000);
    fragViewPos = vec4(0);
    fragWorldPos = vec4(0);
}