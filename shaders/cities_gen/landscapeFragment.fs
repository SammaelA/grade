#version 330 core

uniform vec3 sun_pos;
uniform vec3 sun_color;
uniform vec3 camera_position;

uniform float grass_texture_scale;
uniform float soil_texture_scale;
uniform float sand_texture_scale;
uniform sampler2D grass_texture;
uniform sampler2D soil_texture;
uniform sampler2D sand_texture;

// uniform vec4 diffuse_color;
// uniform bool is_isotropic_color;
// uniform bool is_depth_only;


in vec3 frag_normal;
in vec3 frag_position;
in vec4 frag_sun_coordinates;
// in vec2 frag_tex_coord;
in float slope_ratio;
in float water_ratio;

out vec4 color;

const float AMBIENT_FORCE = 0.6;
const float SPECULAR_RATIO = 0.0;

void main()
{
    // vec4 grassColor = texture(grass_texture, frag_tex_coord * grass_texture_scale);
    // vec4 rockColor = texture(soil_texture, frag_tex_coord * soil_texture_scale);
    // vec4 sandColor = texture(sand_texture, frag_tex_coord * sand_texture_scale);
    vec4 grassColor = texture(grass_texture, frag_position.xz * grass_texture_scale);
    vec4 soilColor = texture(soil_texture, frag_position.xz * soil_texture_scale);
    vec4 sandColor = texture(sand_texture, frag_position.xz * sand_texture_scale);
    float beachFactor = clamp(water_ratio, 0, 1);
    vec4 fragmentColor = mix(mix(grassColor, soilColor, beachFactor), soilColor, slope_ratio);
     
    // vec3 directionToLight = normalize(sun_pos - frag_position);
    //ASSUMING ORTHOGRAPHIC SUN ANGLE
    vec3 directionToLight = normalize(sun_pos);

    vec3 reflectDir = reflect(-directionToLight, frag_normal);
    vec3 viewDir = normalize(frag_position - camera_position);
    float normalLightAngleCos = dot(frag_normal, directionToLight);
    float normalLightAngleCosRatio = max(0.0, normalLightAngleCos);

    vec3 ambientComp = fragmentColor.xyz * AMBIENT_FORCE;
    vec3 diffuseComp = fragmentColor.xyz * normalLightAngleCosRatio;
    float specularComp = pow(max(dot(-viewDir, reflectDir), 0.0), 10.0) * SPECULAR_RATIO;


    float shadowRatio = 0;


    color = vec4((ambientComp + (specularComp + diffuseComp) * (1 - shadowRatio)) * sun_color,1.0);
}
