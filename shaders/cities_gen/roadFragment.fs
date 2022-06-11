#version 330 core

uniform sampler2D asphalt_texture;
uniform sampler2D marking_texture;
uniform sampler2D cross_marking_texture;
// uniform sampler2D shadow_texture;

uniform float asphalt_scale;
uniform vec3 sun_pos;
// uniform vec3 sun_color;
// uniform vec4 diffuse_color;
// uniform bool is_isotropic_color;
// uniform bool is_depth_only;
// uniform vec3 camera_position;
// uniform float roughness_factor;

in vec2 frag_asphalt_tex_coord;
in float frag_marking_tex_coord;
in vec3 frag_global_position;
flat in vec4 frag_anchors_coords;
flat in vec4 frag_handles_coords;
in vec3 frag_normal;
// in vec4 frag_sun_coordinates;

out vec4 color;

const float AMBIENT_FORCE = 0.4;
// const float BASE_SPECULAR_RATIO = 0.6;

float clampCoef(float t)
{
    return (sign(sign(t + 0.01) * sign(1.01 - t) + 0.001) + 1) * 0.5;
}

void main()
{
    // vec4 fragmentColor;
    // if (!is_isotropic_color)
    //     fragmentColor = texture(diffuse_texture, frag_tex_coord);
    // else
    //     fragmentColor = diffuse_color;

    // vec3 directionToLight = normalize(sun_pos);

    // vec3 reflectDir = reflect(-directionToLight, frag_normal);
    // vec3 viewDir = normalize(frag_position - camera_position);
    // float normalLightAngleCos = dot(frag_normal, directionToLight);
    // float normalLightAngleCosRatio = max(0.0, normalLightAngleCos);

    // vec3 ambientComp = fragmentColor.xyz * AMBIENT_FORCE;
    // vec3 diffuseComp = fragmentColor.xyz * normalLightAngleCosRatio;
    // float SPECULAR_RATIO = BASE_SPECULAR_RATIO * (1 - roughness_factor);
    // float specularComp = pow(max(dot(-viewDir, reflectDir), 0.0), 14.0) * SPECULAR_RATIO;
    // float shadowRatio = isShaded(frag_sun_coordinates, shadow_texture, normalLightAngleCos);

    // color = vec4((ambientComp + (specularComp + diffuseComp) * (1 - shadowRatio)) * sun_color,1.0);

    // color = texture(asphalt_texture, frag_global_position.xz * asphalt_scale);
    vec4 asphaltColor = texture(asphalt_texture, frag_asphalt_tex_coord * asphalt_scale);
    vec4 markingColor;
    if (frag_marking_tex_coord <= 1.1)
    {
        markingColor = texture(marking_texture, vec2(frag_marking_tex_coord, 0.5));
    }
    else if (frag_marking_tex_coord < 1.9)
    {
        vec2 startCoordinates = frag_handles_coords.zw;
        vec2 p0 = frag_anchors_coords.xy - startCoordinates;
        vec2 p1 = frag_handles_coords.xy - startCoordinates;
        vec2 p2 = frag_anchors_coords.zw - startCoordinates;
        vec2 localPos = frag_global_position.xz - startCoordinates;

        // localPos.y = sign(localPos.y) * max(abs(localPos.y), 0.00001);
        float k = localPos.x / localPos.y;
        float a = (p0.y * k - p0.x) + (p2.y * k - p2.x) - 2 * (p1.y * k - p1.x);
        float b = 2 * (p1.y * k - p1.x) - 2 * (p0.y * k - p0.x);
        float c = p0.y * k - p0.x;
        float D = pow(b, 2) - 4 * a * c;
        float t1 = (-b + sqrt(D)) / (2 * a);
        float t2 = (-b - sqrt(D)) / (2 * a);
        float t1Ratio = clampCoef(t1);
        float t = clamp(t1Ratio * t1 + (1 - t1Ratio) * t2, 0, 1);
        vec2 bezierPos = pow((1 - t), 2) * p0 + 2 * t * (1 - t) * p1 + pow(t, 2) * p2;
        float r = (length(localPos) / length(bezierPos)) * 0.5;
        r = clamp(r, 0, frag_marking_tex_coord - 1);
        markingColor = texture(marking_texture, vec2(r, 0.5));
    }
    else
    {
        float turnRatio = frag_marking_tex_coord - 2;
        vec2 centerPosition = frag_anchors_coords.xy;
        float smallRadius = mix(frag_anchors_coords.z, frag_anchors_coords.w, turnRatio);
        float bigRadius = mix(frag_handles_coords.x, frag_handles_coords.y, turnRatio);
        float realMarkingTexCoords = clamp(
            (length(frag_global_position.xz - centerPosition) - smallRadius) / (bigRadius - smallRadius),
            0,
            1
        );
        markingColor = texture(cross_marking_texture, vec2(realMarkingTexCoords, 0.5));
    }
    vec4 diffuseColor = vec4(mix(asphaltColor.rgb, markingColor.rgb, markingColor.a), 1);


    vec3 ambientComp = diffuseColor.rgb * AMBIENT_FORCE;

    vec3 directionToLight = normalize(sun_pos);
    float normalLightAngleCos = dot(frag_normal, directionToLight);
    float normalLightAngleCosRatio = max(0.0, normalLightAngleCos);
    vec3 diffuseComp = diffuseColor.rgb * normalLightAngleCosRatio;

    color = vec4(ambientComp + diffuseComp, 1.0);
}
