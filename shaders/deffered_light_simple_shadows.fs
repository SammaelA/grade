#version 330

in vec2 ex_Tex;

uniform sampler2D colorTex;
uniform sampler2D normalsTex;
uniform sampler2D viewPosTex;
uniform sampler2D worldPosTex;
uniform sampler2D aoTex;
uniform sampler2D cubeTex;
uniform int mode;
uniform vec3 dir_to_sun;
uniform vec3 camera_pos;
uniform vec3 ambient_diffuse_specular;
uniform vec3 light_color;
uniform vec2 sts_inv;
uniform sampler2D shadowMap;
uniform bool need_shadow;
uniform mat4 shadow_mat;

out vec4 fragColor;
#define VSM 1

//different types of pixel
//should match with scene_generation_helper.cpp
const int PIXEL_TYPE_NONE = 0;
const int PIXEL_TYPE_TERRAIN = 1;
const int PIXEL_TYPE_TREES = 2;
const int PIXEL_TYPE_GRASS = 3;
const int PIXEL_TYPE_MODELS = 4;
const int PIXEL_TYPE_DEBUG_LIGHT = 5;
const int PIXEL_TYPE_DEBUG_NO_LIGHT = 6;

float ShadowCalculation(vec4 fragPosLightSpace, vec3 world_pos, float bias, int samples)
{
  vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
  projCoords = projCoords * 0.5 + 0.5;

  float res = 0;
  if(projCoords.z > 1.0)
    return 0.0;

    vec2 cD_tD = texture(shadowMap, projCoords.xy).wz;
    float trans_mult = sqrt(sqrt(cD_tD.y));
    float currentDepth = projCoords.z;
    res = currentDepth - bias > cD_tD.x ? trans_mult : 0.0;
  
  #if VSM == 1
    vec2 vsm = texture(shadowMap, projCoords.xy).xy; 
    if (projCoords.z >= vsm.x + bias)
    {
      res += vsm.x;
    }
    res = clamp(res,0,1);
  #endif

  return res;
}

void main(void) 
{
      vec4 color_alpha = texture(colorTex,ex_Tex);
      vec4 normal_type = texture(normalsTex,ex_Tex);
      vec3 cube_color = pow(texture(cubeTex, ex_Tex).xyz,vec3(1/2.2));//remove pow is cube tex is in sRGB
      int type = int(normal_type.a);
      
      if (color_alpha.a < -0.01)
      {
        fragColor = vec4(cube_color,1);
      }
      else if (type != PIXEL_TYPE_DEBUG_NO_LIGHT)
      {
        vec3 world_pos = texture(worldPosTex,ex_Tex).xyz;
        vec3 view_pos = texture(viewPosTex,ex_Tex).xyz;
        float ao = texture(aoTex,ex_Tex).x;
        fragColor = vec4(color_alpha.xyz/color_alpha.a,1);
        fragColor = vec4(pow(color_alpha.xyz,vec3(2.2)),1);
        float shadow = 0;
        if (need_shadow)
        {
            float bias = 0.0001;
            int samples = 1;

            if (type == PIXEL_TYPE_TERRAIN)
                samples = 4;
            else 
                bias = 0.015;
            vec4 FragPosLightSpace = shadow_mat * vec4(world_pos,1);
            shadow = ShadowCalculation(FragPosLightSpace, world_pos, bias, samples);
        }
        vec3 n = normalize(normal_type.xyz);
        vec3 ads = ambient_diffuse_specular;
        n = dot(n,dir_to_sun) > 0 ? n : -n;
        float lambertian = max(dot(dir_to_sun, n), 0.0);
        float specular = 0.0;
        float shininess = 16.0;
        if (lambertian > 0.0) 
        {
          vec3 viewDir = normalize(world_pos - camera_pos);
          vec3 halfDir = normalize(dir_to_sun + viewDir);
          float specAngle = max(dot(halfDir, n), 0.0);
          specular = pow(specAngle, shininess);
        }
        vec3 h = reflect(-dir_to_sun,n);
        vec3 color = vec3(1,1,1)*ads.x + (1-shadow)*(light_color*ads.y*lambertian + vec3(1,1,1)*ads.z*specular);
        fragColor.xyz *= color;
        fragColor.xyz = pow(fragColor.xyz,vec3(1.0/2.2));
        fragColor.xyz = fragColor.xyz*color_alpha.a + cube_color*(1 - color_alpha.a);
      }
      else
      {
        fragColor = color_alpha;
      }
}