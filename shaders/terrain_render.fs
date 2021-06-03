#version 330

in vec3 ex_Tex;
in vec3 ex_Normal;
in vec3 ex_FragPos;

uniform vec3 dir_to_sun;
uniform vec3 camera_pos;
uniform vec3 ambient_diffuse_specular;
uniform vec3 light_color;

in vec4 FragPosLightSpace;
out vec4 fragColor;
uniform sampler2D shadowMap;
uniform bool need_shadow;
#define VSM 1

float ShadowCalculation(vec4 fragPosLightSpace)
{
  vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
  projCoords = projCoords * 0.5 + 0.5;

  #if VSM == 1
    vec2 vsm = texture(shadowMap, projCoords.xy).rg; 
    if (projCoords.z >= vsm.x)
    {
      float mu = vsm.x;
      float s2 = vsm.y - mu*mu;
      float	pmax = s2 / ( s2 + (projCoords.z - mu)*(projCoords.z - mu) );
      return pmax;
    }
    else
      return 0;
  #else
    float closestDepth = texture(shadowMap, projCoords.xy).r; 
    float currentDepth = projCoords.z;
    float bias = 1e-5;

    float shadow = currentDepth - bias > closestDepth ? 1.0 : 0.0;
    if(projCoords.z > 1.0)
        shadow = 0.0;
    return shadow;
  #endif
}

void main(void) 
{
  if (need_shadow)
  {
    fragColor = 0.67*vec4(0.4,0.7,0.1,1) + 0.33*vec4(ex_Tex,1);
    float shadow = need_shadow ? ShadowCalculation(FragPosLightSpace) : 0;
    vec3 n = ex_Normal;
    vec3 ads = ambient_diffuse_specular;
    vec3 dir = normalize(camera_pos - ex_FragPos);
    n = dot(n,dir_to_sun) > 0 ? n : -n;
    vec3 h = reflect(-dir_to_sun,n);
    vec3 color = vec3(1,1,1)*ads.x + light_color*(ads.y*dot(n,dir_to_sun)+0*pow(dot(n,h),16))*(1-shadow);
    fragColor.xyz *= color;
  }
  else
    fragColor = vec4(ex_Tex,1);
}
