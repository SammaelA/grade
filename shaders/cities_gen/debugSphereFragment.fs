#version 330 core

uniform vec3 sun_pos;
uniform vec3 sun_color;
uniform vec4 diffuse_color;

in vec3 frag_normal;

out vec4 color;

const float AMBIENT_FORCE = 0.5;


void main()
{
    vec3 frag_normal1 = normalize(frag_normal);
    vec3 directionToLight = normalize(sun_pos);
    float normalLightAngleCos = dot(frag_normal1, directionToLight);
    float normalLightAngleCosRatio = max(0.0, normalLightAngleCos);

    vec4 fragmentColor = diffuse_color;
    vec3 ambientComp = fragmentColor.xyz * AMBIENT_FORCE;
    vec3 diffuseComp = fragmentColor.xyz * normalLightAngleCosRatio;
    
    color = vec4((ambientComp + diffuseComp) * sun_color,1.0);
}
