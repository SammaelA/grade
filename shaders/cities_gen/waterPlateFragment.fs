#version 330 core

uniform vec2 window_size;
uniform float time_seconds;
uniform vec3 sun_pos;
uniform vec3 sun_color;
uniform vec3 camera_pos;
uniform vec2 z_near_far;
uniform mat4 inv_proj_matrix;


uniform sampler2D under_water_texture;
uniform sampler2D above_water_texture;
uniform sampler2D dudv_texture;
uniform sampler2D normal_texture;
uniform sampler2D bloom_texture;
uniform sampler2D depth_under_water_texture;

// 0 - default
// 2 - bloom
uniform int render_mode;


in vec3 frag_global_pos;

out vec4 color;


const float NORMAL_MAP_SCALE = 0.1;
const float DISTORTION_FORCE = 0.08;
const float WATER_SPEED = 0.5 * NORMAL_MAP_SCALE;
const float SPECULAR_DEGREE = 56;
const float NORMAL_STRAIGHTEN_FORCE = 0.16;
const float REFLECTIVE_BASE = 0.1;
const float SHLICK_DEGREE = 2;
const float SPECULAR_FROM = 0.4;
const int BLOOM_RADIUS = 3;
const float GAUSS_SIGMA = 2.0; // less means narrower bloom
const float BLOOM_DEFORM_FACTOR = 1;
const float BLOOM_DISTANCE = 35;
const float MAX_REFRACTION_WATER_DEPTH = 8;
const vec4 WATER_COLOR = vec4(0, 31, 26, 257) / 257;
const float ARTIFACT_SUPPRESS_MAX_DEPTH = 0.51;
const float MAX_DISTANCE_SUPPRESS = 50;

//derivative
const float GAUSS_MULTIPLIER_OUT = 1 / (GAUSS_SIGMA * 2.50662827463);
const float GAUSS_MULTIPLIER_IN = 1 / (GAUSS_SIGMA * GAUSS_SIGMA * 2);
float zNear = z_near_far.x;
float zFar  = z_near_far.y;

float Random (vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

float gaussCoef(float x)
{
    return GAUSS_MULTIPLIER_OUT * exp(-(x*x) * GAUSS_MULTIPLIER_IN);
}

float bloomWeight (vec2 shift)
{
    return gaussCoef(length(shift));
}

float DistanceToCamera(float depth)
{
    float z_clip = depth * 2 - 1;
    float z_eye_space = - (2.0 * zNear * zFar) / (zFar + zNear - z_clip * (zFar - zNear));
    vec4 vertexOutput = vec4(0);
    vertexOutput.w = -z_eye_space;
    vertexOutput.xyz = vec3(vec2(gl_FragCoord.xy / window_size), depth);
    vertexOutput.xyz = vertexOutput.xyz * 2 - 1;
    vertexOutput.xyz *= vertexOutput.w;
    vec3 eyeCoords = (inv_proj_matrix * vertexOutput).xyz;
    return length(eyeCoords);
}

void main()
{
    //COORDINATE MOVE AND MOVE
    vec2 textureCoords = frag_global_pos.xz * NORMAL_MAP_SCALE;
    textureCoords += normalize(vec2(1,1)) * time_seconds * WATER_SPEED;

    //NORMAL CALCULATION
    vec3 normalMapSample = texture(normal_texture, textureCoords).xyz;
    vec3 normal = vec3(normalMapSample.r * 2 - 1, normalMapSample.b, normalMapSample.g * 2 - 1);
    normal = normalize(normal);
    normal = mix(normal, vec3(0,1,0), NORMAL_STRAIGHTEN_FORCE);
    vec3 viewDir = normalize(frag_global_pos - camera_pos);
    vec3 pseudoNormal = vec3(0,1,0);
    float normalAngleCos = dot(pseudoNormal, -viewDir);

    if (render_mode == 2)
    {
        //SUN SPECULAR
        vec3 directionToLight = normalize(sun_pos);
        vec3 reflectDir = reflect(-directionToLight, normal);
        float specular = max(0.0, dot(reflectDir, -viewDir));
        specular = pow(specular, SPECULAR_DEGREE);
        specular *= 1 - normalAngleCos;
        specular *= clamp(length(camera_pos - frag_global_pos) / BLOOM_DISTANCE, 0, 1);
        specular = (specular > SPECULAR_FROM) ? 1 : 0;
        color = vec4(specular, 0, 0, 1);
    }

    if (render_mode == 0)
    {
        vec2 screenCoords = gl_FragCoord.xy / window_size;
        float distanceToSurface = DistanceToCamera(gl_FragCoord.z);
        float waterDepth = DistanceToCamera(texture(depth_under_water_texture, screenCoords).r);
        waterDepth -= distanceToSurface;

        //COORDINATE DISTORTION
        vec2 dudvDistortion = texture(dudv_texture, textureCoords).rg;
        dudvDistortion -= vec2(0.5);
        dudvDistortion *= DISTORTION_FORCE;
        
        vec2 distortedCoords = screenCoords + dudvDistortion;

        //COORDINATES FOR REFLECTION AND REFRACTION
        vec2 refractionCoords = distortedCoords;
        vec2 reflectionCoords = refractionCoords;
        reflectionCoords.y = 1 - reflectionCoords.y;

        //REFLECTION AND REFRACTION COLORS
        vec4 refractionColor = texture(under_water_texture, refractionCoords);
        vec4 reflectionColor = texture(above_water_texture, reflectionCoords);

        //WATER DEPTH EFFECT
        float depthEffectFactor = clamp(waterDepth / MAX_REFRACTION_WATER_DEPTH, 0, 1);
        refractionColor = mix(refractionColor, WATER_COLOR, depthEffectFactor);

        //FRESNEL EFFECT
        float reflectiveFactor = 
            REFLECTIVE_BASE + (1 - REFLECTIVE_BASE) * pow((1 - normalAngleCos), SHLICK_DEGREE);

        //BLOOM
        vec2 texOffset = 1.0 / window_size;
        float specular = 0;
        for (int i = -BLOOM_RADIUS; i <= BLOOM_RADIUS; i++)
        {
            for (int j = -BLOOM_RADIUS; j <= BLOOM_RADIUS; j++)
            {
                float sample = texture(bloom_texture, screenCoords + texOffset * vec2(i, j)).r;
                float x = length(vec2(i, j)) + abs(i * j) * BLOOM_DEFORM_FACTOR;
                specular += sample * gaussCoef(x);
            }
        }
        specular = clamp(specular, 0, 1);
        vec4 specularColor = vec4(sun_color, 1);

        //MIXING
        vec4 unspecularColor = mix(refractionColor, reflectionColor, reflectiveFactor);
        color = mix(unspecularColor, specularColor, specular);

        //SUPPRESS EDGE ARTIFACT
        float distanceScale = mix(1, MAX_DISTANCE_SUPPRESS, distanceToSurface / zFar);
        color.a = clamp(pow(waterDepth / ARTIFACT_SUPPRESS_MAX_DEPTH, 1) / distanceScale, 0, 1);
    }
}
