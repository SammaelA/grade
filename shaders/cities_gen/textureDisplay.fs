#version 330 core

in vec2 TexCoord;

uniform sampler2D Texture;
uniform float scaleFrom;
uniform float scaleTo;

out vec4 color;

void main()
{
    color = texture(Texture, TexCoord);
    if (TexCoord[0] < 0 || TexCoord[1] < 0 || TexCoord[0] > 1 || TexCoord[1] > 1)
        color.rgb = vec3(1.0, 1.0, 1.0);

    color.rgb = max(color.rgb, scaleFrom);
    color.rgb = min(color.rgb, scaleTo);
    color.rgb = (color.rgb - scaleFrom) / (scaleTo - scaleFrom);
}