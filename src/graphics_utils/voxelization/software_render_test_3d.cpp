#include "voxelization.h"
#include "tinyEngine/camera.h"
#include "third_party/stb_image_write.h"
#include "third_party/stb_image.h"
#include "voxel_array.h"
#include "common_utils/distribution.h"
#include <cstdio>
#include <string>

using glm::degrees;
using glm::mat4;
using glm::vec2;
using glm::vec3;
using glm::vec4;

namespace voxelization
{
  float max(float x, float y)
  {
    return x > y ? x : y;
  }
  float min(float x, float y)
  {
    return x < y ? x : y;
  }
  float clamp(float x, float x_min, float x_max)
  {
    return min(max(x, x_min), x_max);
  }

  struct Image
  {
    static constexpr int channels = 3;
    Image(std::string name, float min_val = 0, float max_val = 1)
    {
      int image_w = 0, image_h = 0, image_channels = 0;
      unsigned char *image_data = stbi_load(name.c_str(), &image_w, &image_h, &image_channels, channels);

      w = image_w;
      h = image_h;
      data = new float[w * h * channels];
      for (int i = 0; i < w * h * channels; i++)
      {
        data[i] = min_val + (image_data[i] / 255.0) * (max_val - min_val);
      }

      stbi_image_free(image_data);
    }
    Image(const Image &image)
    {
      w = image.w;
      h = image.h;
      data = new float[w * h * channels];
      memcpy(data, image.data, sizeof(float) * w * h * channels);
    }
    Image(int _w, int _h) : w(_w), h(_h)
    {
      data = new float[w * h * channels];
    }
    ~Image()
    {
      if (data)
        delete[] data;
    }
    vec3 texelFetch(int x, int y) const
    {
      unsigned pos = channels * (w * clamp(y, 0, h - 1) + clamp(x, 0, w - 1));
      vec3 res = vec3(data[pos], data[pos + 1], data[pos + 2]);
      return res;
    }
    vec3 texture(float x, float y) const
    {
      x = x * w;
      y = y * h;
      int ix = floor(x);
      int iy = floor(y);
      float dx = x - ix;
      float dy = y - iy;
      vec3 p00 = texelFetch(ix, iy);
      vec3 p01 = texelFetch(ix, iy + 1);
      vec3 p10 = texelFetch(ix + 1, iy);
      vec3 p11 = texelFetch(ix + 1, iy + 1);
      vec3 res = (1 - dx) * (1 - dy) * p00 + (1 - dx) * dy * p01 + dx * (1 - dy) * p10 + dx * dy * p11;
      return res;
    }

    void set_pixel(int x, int y, const vec3 &pixel)
    {
      if (x >= 0 && y >= 0 && x < w && y < h)
      {
        unsigned pos = channels * (w * y + x);
        data[pos] = pixel.x;
        data[pos + 1] = pixel.y;
        data[pos + 2] = pixel.z;
      }
    }
    int w, h;
    float *data;
  };

  void save_image(const Image &image, std::string name, float min = 0, float max = 1)
  {
    unsigned char *data = new unsigned char[image.w * image.h * Image::channels];
    for (int i = 0; i < image.w * image.h * Image::channels; i++)
      data[i] = clamp((unsigned)(255 * clamp(image.data[i], min, max)), 0, 255);
    stbi_write_png(name.c_str(), image.w, image.h, Image::channels, data, 0);
    delete[] data;
  }

  void render_3d_scene(VoxelArray<vec4> &voxel_array, CameraSettings &camera,
                       Image &out_image,
                       float image_w, float image_h, float max_distance, int max_steps,
                       int spp)
  {
    mat4 projection = glm::perspective(camera.fov_rad, (float)image_w / image_h, camera.z_near, camera.z_far);
    mat4 view = glm::lookAt(camera.origin, camera.target, camera.up);
    mat4 view_proj = projection * view;
    mat4 view_proj_inv = glm::inverse(view_proj);

    for (int i=0;i<image_w;i++)
    {
      for (int j=0;j<image_h;j++)
      {
        vec3 full_color(0,0,0);
        float full_a = 0;
        for (int sample = 0; sample < spp; sample++)
        {
          vec4 t_pos = vec4(2*(i+urand())/image_w - 1, 2*(j+urand())/image_h - 1, 1, 1);
          vec4 p1 = view_proj_inv * t_pos;
          vec3 ray = glm::normalize(vec3(p1.x/p1.w, p1.y/p1.w, p1.z/p1.w));

          vec3 color(0,0,0);
          float a = 1;
          float dist = (max_distance/max_steps);
          vec3 step = dist*ray;
          for (float k=0;k<max_steps;k++)
          {
            vec3 p = camera.origin + k*step;
            vec4 voxel = voxel_array.get(p);
            float density = voxel.w;
            vec3 c = vec3(voxel.r, voxel.g, voxel.b);
            float a_i = expf(-density*dist);
            color += a*(1-a_i)*c;
            a *= a_i;
            if (a < 1e-4)
              break;
          }
          full_color += color;
          full_a += a;
        }
        out_image.set_pixel(i,j,full_color/(float)spp);
      }
    }
  }
  void software_render_test_3d()
  {
    VoxelArray<vec4> voxel_array(glm::vec3(-15,-15,-15), glm::vec3(15,15,15), glm::ivec3(64,64,64), vec4(0,0,0,0));
    voxel_array.set_circle(vec3(0,0,0), 5, vec4(1,0,0,0.5));

    voxel_array.set_circle(vec3(10,10,10), 5, vec4(0,0,1,0.5));
    voxel_array.set_circle(vec3(10,10,-10), 5, vec4(0,0,1,0.5));
    voxel_array.set_circle(vec3(-10,10,10), 5, vec4(0,0,1,0.5));
    voxel_array.set_circle(vec3(-10,10,-10), 5, vec4(0,0,1,0.5));

    voxel_array.set_circle(vec3(10,-10,10), 5, vec4(0,1,0,0.5));
    voxel_array.set_circle(vec3(10,-10,-10), 5, vec4(0,1,0,0.5));
    voxel_array.set_circle(vec3(-10,-10,10), 5, vec4(0,1,0,0.5));
    voxel_array.set_circle(vec3(-10,-10,-10), 5, vec4(0,1,0,0.5));

    int n = 0;
    CameraSettings camera;
    camera.target = vec3(0, 0, 0);
    camera.up = vec3(0, 1, 0);

    float image_w = 512;
    float image_h = 512;

    Image test_image(image_w, image_h);
    
    for (float t=0;t<=1;t+=0.02)
    {
      camera.origin = vec3(50*sin(2*M_PI*t), 0, 50*cos(2*M_PI*t));
      render_3d_scene(voxel_array, camera, test_image, image_w, image_h, camera.z_far, 256, 1);
      save_image(test_image, "saves/3d_render/test_image.png");
      //save_image(test_image, "saves/3d_render/test_image_"+std::to_string(n)+".png");
      n++;
    }
  }
};