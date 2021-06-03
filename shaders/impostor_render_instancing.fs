#version 430

in vec3 tc_a;
in vec3 tc_b;
in vec3 tc_t;
in vec3 q_abt;
in vec3 ex_Normal;
in vec3 ex_FragPos;
in vec2 a_mult;
in mat4 rot_m;
in mat4 normalTr;
in vec4 FragPosLightSpace;
flat in uint model_id;
out vec4 fragColor;
uniform vec3 dir_to_sun;
uniform vec3 camera_pos;
uniform vec3 ambient_diffuse_specular;
uniform vec3 light_color;
uniform sampler2DArray color_tex;
uniform sampler2DArray normal_tex;
uniform sampler2D noise;
uniform vec4 screen_size;
uniform int debug_model_id;
uniform sampler2D shadowMap;
uniform bool need_shadow;
float ShadowCalculation(vec4 fragPosLightSpace)
{
        // Выполняем деление перспективы
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
 
    // Преобразуем в диапазон [0,1]
    projCoords = projCoords * 0.5 + 0.5;
 
    // Получаем наиболее близкое значение глубины исходя из перспективы глазами источника света (используя диапазон [0,1] fragPosLight в качестве координат)
    float closestDepth = texture(shadowMap, projCoords.xy).r; 
 
    // Получаем глубину текущего фрагмента исходя из перспективы глазами источника света
    float currentDepth = projCoords.z;
 
    // Проверяем, находится ли текущий фрагмент в тени
  float bias = 1e-5;

  float shadow = currentDepth - bias > closestDepth ? 1.0 : 0.0;
  if(projCoords.z > 1.0)
        shadow = 0.0;
    return shadow;
}
float gradientNoise(float x, float y)
{
  float f = 0.06711056f * x + 0.00583715f * y;
  return fract(52.9829189f * fract(f));
}
void main(void) 
{
  vec2 noise_pos = vec2(fract(25*gl_FragCoord.x*screen_size.z), fract(25*gl_FragCoord.y*screen_size.w));
  float ns_imp = texture(noise,noise_pos).x;
  float ns = gradientNoise(float(gl_FragCoord.x),float(gl_FragCoord.y));
  vec3 tc = ns_imp <= q_abt.x ? tc_a : (ns_imp <= q_abt.x + q_abt.y ? tc_b : tc_t);
  fragColor = texture(color_tex,tc);
  if (fragColor.a<0.33)
    discard;

  if ((a_mult.y > 0.1 && ns > a_mult.x) || (a_mult.y < -0.1 && ns <  1  - a_mult.x))
    discard;
  fragColor.rgb /= fragColor.a;
  if (int(model_id) == debug_model_id)
    fragColor = vec4(1,0,1,1);
  if (need_shadow)
  {
    float shadow = need_shadow ? ShadowCalculation(FragPosLightSpace) : 0;
    vec3 n = (normalTr*vec4(texture(normal_tex,tc).xyz,0)).xyz;
    vec3 ads = ambient_diffuse_specular;
    vec3 dir = normalize(camera_pos - ex_FragPos);
    n = dot(n,dir_to_sun) > 0 ? n : -n;
    vec3 h = reflect(-dir_to_sun,n);
    vec3 color = vec3(1,1,1)*ads.x + (1-shadow)*light_color*(ads.y*dot(n,dir_to_sun)+0*pow(dot(n,h),16));
    fragColor.xyz *= color;
  }
}
