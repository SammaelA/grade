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

float ShadowCalculation(vec4 fragPosLightSpace, vec3 world_pos, float bias, int samples)
{
vec2 poissonDisk[64];
poissonDisk[0] = vec2(-0.613392, 0.617481);
poissonDisk[1] = vec2(0.170019, -0.040254);
poissonDisk[2] = vec2(-0.299417, 0.791925);
poissonDisk[3] = vec2(0.645680, 0.493210);
poissonDisk[4] = vec2(-0.651784, 0.717887);
poissonDisk[5] = vec2(0.421003, 0.027070);
poissonDisk[6] = vec2(-0.817194, -0.271096);
poissonDisk[7] = vec2(-0.705374, -0.668203);
poissonDisk[8] = vec2(0.977050, -0.108615);
poissonDisk[9] = vec2(0.063326, 0.142369);
poissonDisk[10] = vec2(0.203528, 0.214331);
poissonDisk[11] = vec2(-0.667531, 0.326090);
poissonDisk[12] = vec2(-0.098422, -0.295755);
poissonDisk[13] = vec2(-0.885922, 0.215369);
poissonDisk[14] = vec2(0.566637, 0.605213);
poissonDisk[15] = vec2(0.039766, -0.396100);
poissonDisk[16] = vec2(0.751946, 0.453352);
poissonDisk[17] = vec2(0.078707, -0.715323);
poissonDisk[18] = vec2(-0.075838, -0.529344);
poissonDisk[19] = vec2(0.724479, -0.580798);
poissonDisk[20] = vec2(0.222999, -0.215125);
poissonDisk[21] = vec2(-0.467574, -0.405438);
poissonDisk[22] = vec2(-0.248268, -0.814753);
poissonDisk[23] = vec2(0.354411, -0.887570);
poissonDisk[24] = vec2(0.175817, 0.382366);
poissonDisk[25] = vec2(0.487472, -0.063082);
poissonDisk[26] = vec2(-0.084078, 0.898312);
poissonDisk[27] = vec2(0.488876, -0.783441);
poissonDisk[28] = vec2(0.470016, 0.217933);
poissonDisk[29] = vec2(-0.696890, -0.549791);
poissonDisk[30] = vec2(-0.149693, 0.605762);
poissonDisk[31] = vec2(0.034211, 0.979980);
poissonDisk[32] = vec2(0.503098, -0.308878);
poissonDisk[33] = vec2(-0.016205, -0.872921);
poissonDisk[34] = vec2(0.385784, -0.393902);
poissonDisk[35] = vec2(-0.146886, -0.859249);
poissonDisk[36] = vec2(0.643361, 0.164098);
poissonDisk[37] = vec2(0.634388, -0.049471);
poissonDisk[38] = vec2(-0.688894, 0.007843);
poissonDisk[39] = vec2(0.464034, -0.188818);
poissonDisk[40] = vec2(-0.440840, 0.137486);
poissonDisk[41] = vec2(0.364483, 0.511704);
poissonDisk[42] = vec2(0.034028, 0.325968);
poissonDisk[43] = vec2(0.099094, -0.308023);
poissonDisk[44] = vec2(0.693960, -0.366253);
poissonDisk[45] = vec2(0.678884, -0.204688);
poissonDisk[46] = vec2(0.001801, 0.780328);
poissonDisk[47] = vec2(0.145177, -0.898984);
poissonDisk[48] = vec2(0.062655, -0.611866);
poissonDisk[49] = vec2(0.315226, -0.604297);
poissonDisk[50] = vec2(-0.780145, 0.486251);
poissonDisk[51] = vec2(-0.371868, 0.882138);
poissonDisk[52] = vec2(0.200476, 0.494430);
poissonDisk[53] = vec2(-0.494552, -0.711051);
poissonDisk[54] = vec2(0.612476, 0.705252);
poissonDisk[55] = vec2(-0.578845, -0.768792);
poissonDisk[56] = vec2(-0.772454, -0.090976);
poissonDisk[57] = vec2(0.504440, 0.372295);
poissonDisk[58] = vec2(0.155736, 0.065157);
poissonDisk[59] = vec2(0.391522, 0.849605);
poissonDisk[60] = vec2(-0.620106, -0.328104);
poissonDisk[61] = vec2(0.789239, -0.419965);
poissonDisk[62] = vec2(-0.545396, 0.538133);
poissonDisk[63] = vec2(-0.178564, -0.596057);


  vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
  projCoords = projCoords * 0.5 + 0.5;

  float res = 0;
  if(projCoords.z > 1.0)
    return 0.0;

    float base_depth = texture(shadowMap, projCoords.xy).b;
    float currentDepth = projCoords.z;
    float delta = 1 + 30*smoothstep(0.0,0.006,currentDepth - base_depth);
    int start = int(13*fract(40*world_pos.x + 41*world_pos.y + 51*world_pos.z));

    for (int i = 0;i< samples;i++)
    {
      vec2 cD_tD = texture(shadowMap, projCoords.xy + delta*sts_inv*poissonDisk[(start + 13*i) % 64]).zw; 
      float trans_mult = sqrt(sqrt(cD_tD.y));
      float shadow = currentDepth - bias > cD_tD.x ? trans_mult : 0.0;
      res += shadow/samples;
    }
  
  #if VSM == 1
    vec3 vsm = texture(shadowMap, projCoords.xy).rgb; 
    if (projCoords.z >= vsm.x + bias)
    {
      float mu = vsm.x;
      float s2 = vsm.y - mu*mu;
      float	pmax = s2 / ( s2 + (projCoords.z - mu)*(projCoords.z - mu) );
      res += 0.25*clamp(pmax,0,1);
    }
    res = clamp((1/1.25)*res,0,1);
  #endif

  return res;
}

void main(void) 
{
    if (mode == 0)
    {
      vec4 color_none = texture(colorTex,ex_Tex);
      vec4 normal_type = texture(normalsTex,ex_Tex);
      vec3 cube_color = texture(cubeTex, ex_Tex).xyz;
      if (color_none.a < 0.01)
      {
        fragColor = vec4(cube_color,1);
      }
      else
      {
        vec3 world_pos = texture(worldPosTex,ex_Tex).xyz;
        vec3 view_pos = texture(viewPosTex,ex_Tex).xyz;
        float ao = texture(aoTex,ex_Tex).x;
        fragColor = vec4(color_none.xyz/color_none.a,1);
        fragColor = vec4(pow(color_none.xyz,vec3(2.2)),1);
        float shadow = 0;
        if (need_shadow)
        {
            float bias = 1e-6;
            int samples = 2;
            int type = int(normal_type.a);

            if (type == 0)
                samples = 16;
            else 
                bias = 3*1e-6;
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
        fragColor.xyz = fragColor.xyz*color_none.a + cube_color*(1 - color_none.a);
      }
    }
    else
    {
        fragColor = texture(colorTex,ex_Tex);
    }
}