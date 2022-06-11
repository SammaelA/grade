#version 330 core

uniform sampler2D diffuse_texture;
uniform sampler2D shadow_texture;

uniform vec3 sun_pos;
uniform vec3 sun_color;
uniform vec4 diffuse_color;
uniform float specular_ratio;
uniform bool is_isotropic_color;
uniform bool is_depth_only;
uniform vec3 camera_position;
uniform float roughness_factor;

in vec2 frag_tex_coord;
in vec3 frag_normal;
in vec3 frag_position;
in vec4 frag_sun_coordinates;

out vec4 color;

const float AMBIENT_FORCE = 0.4;
const float BASE_SPECULAR_RATIO = 0.6;

float isShaded(vec4 sunCoords, sampler2D shadowTex, float lightAngleCos)
{
    const float BIAS_BASE = 0.001;
    // float BIAS = max(BIAS_BASE * (1.0 - lightAngleCos), BIAS_BASE * 0.1);
    float tanAlpha = sqrt(1 - lightAngleCos * lightAngleCos) / lightAngleCos;
    float BIAS = BIAS_BASE * tanAlpha + BIAS_BASE * 2;

    vec3 resultSunCoordinates = sunCoords.rgb / sunCoords.w * 0.5 + 0.5;
    float sunDepthFromMap = texture(shadowTex, resultSunCoordinates.xy).r;
    float sunDepthReal = resultSunCoordinates.z;
    return (sunDepthFromMap + BIAS < sunDepthReal) ? 1.0 : 0.0;
}

void main()
{
    if (is_depth_only)
    {
        // float zDiff = zFar - zNear; 
        // float interpolatedDepth = 
        //     (test.w/test.z) * zFar * zNear/zDiff + 0.5 * (zFar + zNear)/zDiff + 0.5;    
        color = vec4(gl_FragCoord.z, 0, 0, 1);
    }
    else
    {
        vec4 fragmentColor;
        if (!is_isotropic_color)
            fragmentColor = texture(diffuse_texture, frag_tex_coord);
        else
            fragmentColor = diffuse_color;
        
        // vec3 directionToLight = normalize(sun_pos - frag_position);
        //ASSUMING ORTHOGRAPHIC SUN ANGLE
        vec3 directionToLight = normalize(sun_pos);

        vec3 reflectDir = reflect(-directionToLight, frag_normal);
        vec3 viewDir = normalize(frag_position - camera_position);
        float normalLightAngleCos = dot(frag_normal, directionToLight);
        float normalLightAngleCosRatio = max(0.0, normalLightAngleCos);

        vec3 ambientComp = fragmentColor.xyz * AMBIENT_FORCE;
        vec3 diffuseComp = fragmentColor.xyz * normalLightAngleCosRatio;
        float SPECULAR_RATIO = BASE_SPECULAR_RATIO * (1 - roughness_factor);
        float specularComp = pow(max(dot(-viewDir, reflectDir), 0.0), 14.0) * SPECULAR_RATIO;
        float shadowRatio = isShaded(frag_sun_coordinates, shadow_texture, normalLightAngleCos);
        
        color = vec4((ambientComp + (specularComp + diffuseComp) * (1 - shadowRatio)) * sun_color,1.0);
    }
}
