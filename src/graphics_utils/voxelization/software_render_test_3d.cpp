#include "voxelization.h"
#include "tinyEngine/camera.h"
#include "third_party/stb_image_write.h"
#include "third_party/stb_image.h"
#include "voxel_array.h"
#include "common_utils/distribution.h"
#include <cppad/cppad.hpp>
#include "diff_generators/vectors.h"
#include "common_utils/optimization/optimization.h"
#include "common_utils/blk.h"
#include <cstdio>
#include <string>
#include <chrono>

using glm::degrees;
using glm::mat4;
using glm::vec2;
using glm::vec3;
using glm::vec4;

using dgen::dfloat;
using dgen::dvec4;
using dgen::dvec3;

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
    Image& operator=(const Image& image)
    {
      if (this == &image)
        return *this;
 
      if (data)
        delete[] data;
        
      w = image.w;
      h = image.h;
      data = new float[w * h * channels];
      memcpy(data, image.data, sizeof(float) * w * h * channels);

      return *this;
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


  void render_reference_image(CameraSettings &camera, Image &out_image,
                              float image_w, float image_h)
  {
    mat4 projection = glm::perspective(camera.fov_rad, (float)image_w / image_h, camera.z_near, camera.z_far);
    mat4 view = glm::lookAt(camera.origin, camera.target, camera.up);
    mat4 view_proj = projection * view;
    mat4 view_proj_inv = glm::inverse(view_proj);

    float max_distance = camera.z_far;
    float max_steps = 512;

    for (int i=0;i<image_w;i++)
    {
      for (int j=0;j<image_h;j++)
      {
        vec3 full_color(0,0,0);
        float full_a = 0;
        for (int sample = 0; sample < 1; sample++)
        {
          vec4 t_pos = vec4(2*(i+0.5)/image_w - 1, 2*(j+0.5)/image_h - 1, 1, 1);
          vec4 p1 = view_proj_inv * t_pos;
          vec3 ray = glm::normalize(vec3(p1.x/p1.w, p1.y/p1.w, p1.z/p1.w));

          vec3 color(0,0,0);
          float a = 1;
          float dist = (max_distance/max_steps);
          vec3 step = dist*ray;
          for (float k=0;k<max_steps;k++)
          {
            vec3 p = camera.origin + k*step;
            float density = 0;
            vec3 c = vec3(0,0,0);
            if (abs(p.x - 0) < 3 && abs(p.y - 0) < 3 && abs(p.z - 0) < 3)
            {
              //printf("[%d %d]\n", i, j);
              density = 1;
              c = vec3(0,1,0);
            }
            float a_i = expf(-density*dist);
            color += a*(1-a_i)*c;
            a *= a_i;
            if (a < 1e-4)
              break;
          }
          full_color += color;
          full_a += a;
        }
        out_image.set_pixel(i,j,full_color);
      }
    }
  }

  void render_3d_scene(VoxelArray<vec4> &voxel_array, const CameraSettings &camera,
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
          vec2 spd = spp == 1 ? vec2(0.5, 0.5) : vec2(urand(), urand());
          vec4 t_pos = vec4(2*(i+spd.x)/image_w - 1, 2*(j+spd.y)/image_h - 1, 1, 1);
          vec4 p1 = view_proj_inv * t_pos;
          vec3 ray = glm::normalize(vec3(p1.x/p1.w, p1.y/p1.w, p1.z/p1.w));

          vec3 color(0,0,0);
          float a = 1;
          float dist = (max_distance/max_steps);
          vec3 step = dist*ray;
          for (float k=0;k<max_steps;k++)
          {
            vec3 p = camera.origin + k*step;
            vec4 voxel = voxel_array.get_trilinear(p);
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

  void diff_render_3d_naive(vec3 p0, vec3 p1, glm::ivec3 voxel_array_size,
                            const std::vector<CameraSettings> &cameras,
                            const std::vector<Image> &reference_images,
                            float image_w, float image_h, float max_distance, int max_steps,
                            int spp, int iterations, std::string filename)
  {

    size_t x_n = 4 * voxel_array_size.x * voxel_array_size.y * voxel_array_size.z;
    std::vector<dfloat> X(x_n, 0);
    std::vector<dfloat> Y;

    // declare independent variables and start recording operation sequence
    CppAD::Independent(X);

    std::vector<dvec4> voxel_data(x_n/4, dvec4(0,0,0,0));

    VoxelArray<dvec4> voxel_array(p0, p1, voxel_array_size, dvec4(0,0,0,0), voxel_data.data());
    for (int i=0;i<x_n/4;i++)
      voxel_array.set_direct(i, dvec4(X[4*i], X[4*i+1], X[4*i+2], X[4*i+3]));

    Y.push_back(0);

    int ref_cnt = min(cameras.size(), reference_images.size());

    for (int ref_n = 0; ref_n < ref_cnt; ref_n++)
    {
      auto &camera = cameras[ref_n];
      auto &reference = reference_images[ref_n];

      mat4 projection = glm::perspective(camera.fov_rad, (float)image_w / image_h, camera.z_near, camera.z_far);
      mat4 view = glm::lookAt(camera.origin, camera.target, camera.up);
      mat4 view_proj = projection * view;
      mat4 view_proj_inv = glm::inverse(view_proj);

      for (int i=0;i<image_w;i++)
      {
        for (int j=0;j<image_h;j++)
        {
          dvec3 full_color(0,0,0);
          dfloat full_a = 0;
          for (int sample = 0; sample < spp; sample++)
          {
            vec4 t_pos = vec4(2*(i+0.5)/image_w - 1, 2*(j+0.5)/image_h - 1, 1, 1);
            vec4 p1 = view_proj_inv * t_pos;
            vec3 ray = glm::normalize(vec3(p1.x/p1.w, p1.y/p1.w, p1.z/p1.w));

            dvec3 color(0,0,0);
            dfloat a = 1;
            float dist = (max_distance/max_steps);
            vec3 step = dist*ray;
            for (float k=0;k<max_steps;k++)
            {
              vec3 p = camera.origin + k*step;
              dvec4 voxel = voxel_array.get_trilinear(p);
              dfloat density = voxel.w;
              dvec3 c = dvec3(voxel.x, voxel.y, voxel.z);
              dfloat a_i = exp(-density*dist);
              color += a*(1-a_i)*c;
              a *= a_i;
              if (a < 1e-4)
                break;
            }
            full_color += color;
            full_a += a;
          }
          vec3 ref_c = reference.texelFetch(i, j);
          dvec3 color_diff = full_color/spp - dvec3(ref_c.x, ref_c.y, ref_c.z);
          Y[0] += (color_diff.x*color_diff.x + color_diff.y*color_diff.y + color_diff.z*color_diff.z)/(image_w*image_h*ref_cnt);
        }
      }
    }
    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y);

    VoxelArray<vec4> res_voxel_array(p0, p1, voxel_array_size, vec4(0,0,0,0));
    Image res_image(image_w, image_h);

    int it = 0;

    opt::opt_func_with_grad F_to_optimize = [&](std::vector<float> &params) -> std::pair<float, std::vector<float>>
    {
      if (it % 10 == 0)
      {
        for (int i=0;i<x_n/4;i++)
          res_voxel_array.set_direct(i, vec4(params[4*i], params[4*i+1], params[4*i+2], params[4*i+3]));
        render_3d_scene(res_voxel_array, cameras[0], res_image, image_w, image_h, max_distance, max_steps, spp);
        save_image(res_image, "saves/3d_render/res_image.png");
      }
      it++;
      return {f.Forward(0, params)[0], f.Jacobian(params)};
    };

    std::vector<float> X0(x_n, 0);
    std::vector<float> params_min(x_n, 0);
    std::vector<float> params_max(x_n, 1);
    for (int i=0;i<x_n;i++)
      X0[i] = urand();

    Block optimizer_settings;
    optimizer_settings.add_arr("initial_params", X0);
    optimizer_settings.add_double("learning_rate", 0.01);
    optimizer_settings.add_int("iterations", iterations);
    optimizer_settings.add_bool("verbose", true);
    opt::Optimizer *optimizer = new opt::Adam();

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    optimizer->optimize(F_to_optimize, params_min, params_max, optimizer_settings);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float sec = 1e-4 * std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    debug("Optimization finished, %.2f s/iter\n", sec/iterations);

    {
      std::vector<float> X_best = optimizer->get_best_result();
      for (int i=0;i<x_n/4;i++)
        res_voxel_array.set_direct(i, vec4(X_best[4*i], X_best[4*i+1], X_best[4*i+2], X_best[4*i+3]));
      render_3d_scene(res_voxel_array, cameras[0], res_image, image_w, image_h, max_distance, max_steps, spp);
      save_image(res_image, "saves/3d_render/res_image.png");
      res_voxel_array.write_to_binary_file(filename);
      }

    delete optimizer;
  }

  void render_reference_1()
  {
    int n = 0;
    CameraSettings camera;
    camera.target = vec3(0, 0, 0);
    camera.up = vec3(0, 1, 0);

    float image_w = 512;
    float image_h = 512;

    Image test_image(image_w, image_h);

    int psi_cnt = 8;
    int phi_cnt = 8;
    float dist = 25;

    for (int psi = 1; psi < psi_cnt; psi++)
    {
      for (int phi = 0; phi < phi_cnt; phi++)
      {
        float phi_f = 2*M_PI*phi/(float)phi_cnt;
        float psi_f = M_PI*psi/(float)psi_cnt - 0.5*M_PI;
        camera.origin = dist*vec3(cos(psi_f)*sin(phi_f), sin(psi_f), cos(psi_f)*cos(phi_f));
        render_reference_image(camera, test_image, image_w, image_h);
        save_image(test_image, "saves/3d_render/reference_image_"+std::to_string(psi)+"_"+std::to_string(phi)+".png");
      }
    }
  }

  void render_test_3d()
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

  void diff_render_naive_test_1(std::string filename, glm::ivec3 voxel_size)
  {
    CameraSettings camera;
    camera.target = vec3(0, 0, 0);
    camera.up = vec3(0, 1, 0);
    camera.origin = vec3(0,0,25);

    float image_w = 256;
    float image_h = 256;

    Image test_image(image_w, image_h);

    render_reference_image(camera, test_image, image_w, image_h);

    save_image(test_image, "saves/3d_render/reference_image.png");

    diff_render_3d_naive(vec3(-5,-5,-5), vec3(5,5,5), voxel_size, {camera}, {test_image},
                         image_w, image_h, 100, 256, 1, 250, filename);
  }

  void diff_render_naive_test_2(std::string filename, glm::ivec3 voxel_size)
  {
    std::vector<CameraSettings> cameras;
    std::vector<Image> images;

    int psi_cnt = 2;
    int phi_cnt = 6;
    float dist = 25;

    float image_w = 128;
    float image_h = 128;

    for (int psi = 1; psi < psi_cnt; psi++)
    {
      for (int phi = 0; phi < phi_cnt; phi++)
      {
        CameraSettings camera;
        camera.target = vec3(0, 0, 0);
        camera.up = vec3(0, 1, 0);
        float phi_f = 2*M_PI*phi/(float)phi_cnt;
        float psi_f = M_PI*psi/(float)psi_cnt - 0.5*M_PI;
        camera.origin = dist*vec3(cos(psi_f)*sin(phi_f), sin(psi_f), cos(psi_f)*cos(phi_f));

        images.emplace_back(image_w, image_h);
        cameras.push_back(camera);

        render_reference_image(camera, images.back(), image_w, image_h);
        save_image(images.back(), "saves/3d_render/reference_image_"+std::to_string(psi)+"_"+std::to_string(phi)+".png");
      }
    }

    diff_render_3d_naive(vec3(-5,-5,-5), vec3(5,5,5), voxel_size, cameras, images,
                         image_w, image_h, 35, 256, 1, 200, filename);
  }

  void render_saved_voxel_array(std::string filename, glm::ivec3 vox_size, int image_size)
  {
    VoxelArray<vec4> voxel_array(vec3(-5,-5,-5), vec3(5,5,5), vox_size, vec4(0));
    voxel_array.read_from_binary_file(filename);
    int n = 0;
    CameraSettings camera;
    camera.target = vec3(0, 0, 0);
    camera.up = vec3(0, 1, 0);

    float image_w = image_size;
    float image_h = image_size;

    Image test_image(image_w, image_h);
    
    for (float t=0;t<=1;t+=0.02)
    {
      camera.origin = vec3(25*sin(2*M_PI*t), 0, 25*cos(2*M_PI*t));
      render_3d_scene(voxel_array, camera, test_image, image_w, image_h, 35, 256, 1);
      save_image(test_image, "saves/3d_render/test_image.png");
      save_image(test_image, "saves/3d_render/test_image_"+std::to_string(n)+".png");
      n++;
    }
  }

  void software_render_test_3d()
  {
    diff_render_naive_test_2("saves/3d_render/res_array_32.bin", glm::ivec3(32,32,32));
    render_saved_voxel_array("saves/3d_render/res_array_32.bin", glm::ivec3(32,32,32), 256);
  }
};