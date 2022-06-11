#version 330 core

uniform sampler2D building_texture;
uniform sampler2D roof_texture;
uniform sampler2D specular_texture;

uniform vec3 sun_pos;
uniform vec3 sun_color;
uniform vec3 camera_position;

in vec3 frag_normal;
in vec3 frag_position;
in vec3 frag_tex_coord;

out vec4 color;


const float AMBIENT_FORCE = 0.4;
const float SPECULAR_RATIO = 0.8;

void main()
{
    // vec4 fragmentColor = vec4(1, 0, 1, 1);
    vec4 fragmentColor;
    float specularColor;
    if (frag_tex_coord.z < 0)
    {
        fragmentColor = texture(building_texture, frag_tex_coord.xy);
        specularColor = texture(specular_texture, frag_tex_coord.xy). r;
    }
    else
    {
        fragmentColor = texture(roof_texture, frag_tex_coord.xy);
        specularColor = 0.f;
    }

    vec3 directionToLight = normalize(sun_pos);
    vec3 reflectDir = reflect(-directionToLight, frag_normal);
    vec3 viewDir = normalize(frag_position - camera_position);
    float normalLightAngleCos = dot(frag_normal, directionToLight);
    float normalLightAngleCosRatio = max(0.0, normalLightAngleCos);

    vec3 ambientComp = fragmentColor.xyz * AMBIENT_FORCE;
    vec3 diffuseComp = fragmentColor.xyz * normalLightAngleCosRatio;
    float specularComp = specularColor * pow(max(dot(-viewDir, reflectDir), 0.0), 60.0) * SPECULAR_RATIO;
    
    color = vec4((ambientComp + diffuseComp + specularComp) * sun_color,1.0);
}
