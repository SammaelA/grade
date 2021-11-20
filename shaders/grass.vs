#version 330

in vec3 in_Position;
in vec3 in_Normal;
in vec4 in_Tex;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 lightSpaceMatrix;

uniform sampler2D hmap;
uniform sampler2D perlin;
uniform sampler2D noise;

out vec3 ex_Tex;
out vec3 ex_Normal;
out vec3 ex_FragPos;
out vec4 ex_FragPosView;
out vec4 FragPosLightSpace;

void main(void) 
{
    vec3 hmap_scales_inv = vec3(1.0/1024,1e3,1.0/1024);
    vec3 hmap_center = vec3(0,0,0);
    vec3 shift = 5.12*vec3(gl_InstanceID / 200 - 100, 0, gl_InstanceID % 200 - 100);
    vec2 tc = 0.5*((shift*hmap_scales_inv).zx + 1);
    float per = texture(perlin,tc).x;
    vec3 noise = 7*(2*texture(noise,tc).xyz - 1);
    if (per > 0.25)
    {
        float h = hmap_scales_inv.y*texture(hmap,tc).x;
        shift.y = h;
        float scale = 7*(per + 1);
        shift.xz += noise.xz;
        vec2 sin_cos = vec2(sin(noise.y),cos(noise.y));
        vec3 pos = vec3(sin_cos.y*in_Position.x - sin_cos.x*in_Position.z,in_Position.y,sin_cos.x*in_Position.x + sin_cos.y*in_Position.z);
        vec3 n = vec3(sin_cos.y*in_Normal.x + sin_cos.x*in_Normal.z, in_Normal.y, -sin_cos.x*in_Normal.x + sin_cos.y*in_Normal.z);
        ex_FragPos = (vec4(scale*pos + shift, 1.0f)).xyz;
        ex_Normal = normalize(n.xyz);
        ex_FragPosView = view * vec4(ex_FragPos, 1.0f);
        gl_Position = projection * ex_FragPosView;
        FragPosLightSpace = lightSpaceMatrix * vec4(ex_FragPos, 1.0);
        ex_Tex.xy = in_Tex.xy;
        ex_Tex.z = 0;
    }
    else
    {
        ex_FragPos = vec3(0,0,0);
        ex_Normal = vec3(0,0,0);
        ex_Tex = vec3(0,0,0);
        gl_Position = vec4(0,0,-1000,1);
    }
}