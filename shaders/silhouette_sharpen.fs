#version 330

in vec2 ex_Tex;
in vec3 ex_FragPos;

out vec4 fragColor;

uniform sampler2DArray tex;
uniform float layer;
uniform float threshold;

void main(void) 
{
    vec4 original = texture(tex, vec3(ex_Tex, layer));
    fragColor = vec4(0,0,0,0);
    fragColor.x = original.x > threshold ? 1 : 0;
    fragColor.y = original.y > threshold ? 1 : 0;
    fragColor.z = max(fragColor.x, fragColor.y);
    fragColor.w = fragColor.x + fragColor.y > 0 ? 1 : 0;
}