#include "graphics_utils/voxelization/voxelization.h"
#include "tinyEngine/camera.h"
#include "third_party/stb_image_write.h"
#include "third_party/stb_image.h"
#include "graphics_utils/voxelization/voxel_array.h"
#include "common_utils/distribution.h"
#include <cppad/cppad.hpp>
#include "diff_generators/vectors.h"
#include "common_utils/optimization/optimization.h"
#include "common_utils/blk.h"
#include "tinyEngine/resources.h"
#include "tinyEngine/engine.h"
#include "tree_utils/tree_modeling.h"
#include "tinyEngine/shader.h"
#include "tinyEngine/model.h"
#include "diff_generators/diff_geometry_generation.h"
#include "diff_generators/differentiable_generators.h"
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

  void render_reference_image_cup(CameraSettings &camera, float image_w, float image_h, 
                                  std::string save_path, Model *m)
  {
    Shader UV({"uv_coords.vs", "simple_render_diffuse.fs"}, {"in_Position", "in_Normal", "in_Tex"});
    GLuint fbo = create_framebuffer();
    Texture t = engine::textureManager->create_texture(image_w, image_h);
    Texture porcelain = engine::textureManager->load_unnamed_tex("resources/textures/porcelain_3.png", 1);

    int w = t.get_W();
    int h = t.get_H();
    //bind FBO
    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    float borderColorDepth[] = {1.0f, 1.0f, 1.0f, 1.0f};
    Texture photo_transformed = engine::textureManager->create_texture(w, h);
    Texture UV_tex = engine::textureManager->create_texture(w, h, GL_RGB32F);
    Texture depthTex = engine::textureManager->create_texture(w, h, GL_DEPTH_COMPONENT16, 1, NULL, GL_DEPTH_COMPONENT, GL_FLOAT);
    glBindTexture(GL_TEXTURE_2D, depthTex.texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColorDepth);
    glBindTexture(GL_TEXTURE_2D, 0);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTex.texture, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, UV_tex.texture, 0);
    glViewport(0, 0, w, h);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glm::mat4 projection = glm::perspective(camera.fov_rad, (float)image_w / image_h, camera.z_near, camera.z_far);
    glm::mat4 view = glm::lookAt(0.15f*camera.origin, camera.target, -camera.up); //to save to file without inversion

    UV.use();
    UV.uniform("projection", projection);
    UV.uniform("view", view);
    UV.texture("tex", porcelain);
    m->render();

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);

    engine::textureManager->save_png_directly(UV_tex, save_path);
    engine::view->next_frame();

    delete_framebuffer(fbo);
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

  void render_reference_image_circle(CameraSettings &camera, Image &out_image,
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
            if (length(p - vec3(0,0,0)) < 4)
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
                            int spp, int iterations, float lr, std::string filename)
  {

    size_t x_n = 4 * voxel_array_size.x * voxel_array_size.y * voxel_array_size.z;
    size_t params_cnt = x_n + 2;
    VoxelArray<vec4> res_voxel_array(p0, p1, voxel_array_size, vec4(0,0,0,0));
    Image res_image(image_w, image_h);
    std::vector<float> last_params;
    int it = 0;
    opt::opt_func_with_grad F_to_optimize = [&](std::vector<float> &params) -> std::pair<float, std::vector<float>>
    {
    std::vector<dfloat> X(params_cnt, 0);
    std::vector<dfloat> Y;
    // declare independent variables and start recording operation sequence
    CppAD::Independent(X);

    std::vector<dvec4> voxel_data(x_n/4, dvec4(0,0,0,0));

    VoxelArray<dvec4> voxel_array(p0, p1, voxel_array_size, dvec4(0,0,0,0), voxel_data.data());
    for (int i=0;i<x_n/4;i++)
      voxel_array.set_direct(i, dvec4(X[4*i], X[4*i+1], X[4*i+2], X[4*i+3]));

    Y.push_back(0);

    //MSE
    int ref_cnt = min(cameras.size(), reference_images.size());
    ref_cnt = 8;
    for (int y = 0; y < ref_cnt; y++)
    {
      int ref_n = urandi(0, cameras.size());
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
            float dx = CppAD::Value(CppAD::Var2Par(X[x_n]));
            float dy = CppAD::Value(CppAD::Var2Par(X[x_n+1]));
            vec4 t_pos = vec4(2*(i+dx)/image_w - 1, 2*(j+dy)/image_h - 1, 1, 1);
            vec4 p1 = view_proj_inv * t_pos;
            vec3 ray = glm::normalize(vec3(p1.x/p1.w, p1.y/p1.w, p1.z/p1.w));

            dvec3 color(0,0,0);
            dfloat a = 1;
            float dist = (max_distance/max_steps);
            vec3 step = dist*ray;
            for (float k=0;k<max_steps;k++)
            {
              vec3 p = camera.origin + k*step;
              dvec4 voxel = voxel_array.get(p);
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

    //regularization
    dfloat q = 0.001;
    for (int i=0;i<voxel_array_size.x;i++)
    {
      for (int j=0;j<voxel_array_size.y;j++)
      {
        Y[0] += q*voxel_array.get_direct(glm::ivec3(i,j,0)).w;
        Y[0] += q*voxel_array.get_direct(glm::ivec3(i,j,voxel_array_size.z - 1)).w;
      }
    }
    for (int i=0;i<voxel_array_size.y;i++)
    {
      for (int j=0;j<voxel_array_size.z;j++)
      {
        Y[0] += q*voxel_array.get_direct(glm::ivec3(0,i,j)).w;
        Y[0] += q*voxel_array.get_direct(glm::ivec3(voxel_array_size.x - 1,i,j)).w;
      }
    }
    for (int i=0;i<voxel_array_size.x;i++)
    {
      for (int j=0;j<voxel_array_size.z;j++)
      {
        Y[0] += q*voxel_array.get_direct(glm::ivec3(i,0,j)).w;
        Y[0] += q*voxel_array.get_direct(glm::ivec3(i,voxel_array_size.y - 1,j)).w;
      }
    }

    size_t y_n = Y.size();
    CppAD::ADFun<float> f(X, Y);

      if (it % 10 == 0)
      {
        for (int i=0;i<x_n/4;i++)
          res_voxel_array.set_direct(i, vec4(params[4*i], params[4*i+1], params[4*i+2], params[4*i+3]));
        render_3d_scene(res_voxel_array, cameras[0], res_image, image_w, image_h, max_distance, max_steps, spp);
        save_image(res_image, "saves/3d_render/res_image.png");
        save_image(res_image, "saves/3d_render/res_image_"+std::to_string(it)+".png");
      }
      if (it == iterations - 1)
        last_params = params;
      it++;
      params[x_n] = urand();
      params[x_n+1] = urand();
      return {f.Forward(0, params)[0], f.Jacobian(params)};
    };

    std::vector<float> X0(params_cnt, 0);
    std::vector<float> params_min(params_cnt, 0);
    std::vector<float> params_max(params_cnt, 10);
    for (int i=0;i<params_cnt;i++)
      X0[i] = urand();


    VoxelArray<vec4> voxel_array(vec3(-5,-5,-5), vec3(5,5,5), voxel_array_size, vec4(0));
    voxel_array.read_from_binary_file(filename);
    for (int i=0;i<x_n/4;i++)
    {
      vec4 r = voxel_array.get_direct(i);
      X0[4*i] = r.x;
      X0[4*i+1] = r.y;
      X0[4*i+2] = r.z;
      X0[4*i+3] = r.w;
    }

    Block optimizer_settings;
    optimizer_settings.add_arr("initial_params", X0);
    optimizer_settings.add_double("learning_rate", lr);
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

  void render_test_3d_2()
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

    float image_w = 2048;
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
                         image_w, image_h, 100, 256, 1, 250, 0.1, filename);
  }

  void diff_render_naive_test_2(std::string filename, glm::ivec3 voxel_size, int image_size)
  {
    std::vector<CameraSettings> cameras;
    std::vector<Image> images;

    int psi_cnt = 4;
    int phi_cnt = 16;
    float dist = 25;

    float image_w = image_size;
    float image_h = image_size;

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
                         image_w, image_h, 35, 256, 1, 250, 0.05, filename);
  }

  void diff_render_naive_test_3(std::string filename, glm::ivec3 voxel_size)
  {
    CameraSettings camera;
    camera.target = vec3(0, 0, 0);
    camera.up = vec3(0, 1, 0);
    camera.origin = vec3(0,0,17);

    Image test_image("saves/3d_render/cup_test.png");

    float image_w = test_image.w;
    float image_h = test_image.h;

    diff_render_3d_naive(vec3(-5,-5,-5), vec3(5,5,5), voxel_size, {camera}, {test_image},
                         image_w, image_h, 100, 256, 1, 250, 0.1, filename);
  }

  void diff_render_naive_test_4(std::string filename, glm::ivec3 voxel_size, int image_size)
  {
    dgen::GeneratorDescription gd = dgen::get_generator_by_name("dishes");
    Model *m = new Model();
    std::vector<float> cup_params{3.318, 3.634, 3.873, 4.091, 4.220, 4.309, 4.356, 4.4, 4.4, 1.0, 0.987, 0.041, 0.477, 0.223, 0.225, 0.277, 0.296, 0.310, 0.317, 0.310, 0.302, 0.296, 0.289, 0.287, 0.286, 0.291, 0.295, 0.306, 0.317, 0.326, 0.345, 0.364, 0.375, 1.991, 1.990, 1.209, 1.221, 1.190, 1.201, 1.153, 1.131, 1.088, 1.048, 0.998, 0.952, 0.899, 0.853, 0.826, 0.842, 0.831, 0.930, 1.098, 1.690}; 
    std::vector<float> scene_params{0, 0, 0, 0, 0, 0, 510.403, 6.013, 26.971, 1.000, 56.560, 0.9, 0.125};
    dgen::DFModel res;
    dgen::dgen_test("dishes", cup_params, res, false, dgen::ModelQuality(false, 1));
    dgen::transform_by_scene_parameters(scene_params, res.first);
    visualizer::simple_mesh_to_model_332(res.first, m);
    m->update();

    std::vector<CameraSettings> cameras;
    std::vector<Image> images;

    int psi_cnt = 1;
    int phi_cnt = 50;
    float dist = 17;

    float image_w = image_size;
    float image_h = image_size;

    for (int psi = 0; psi < psi_cnt; psi++)
    {
      for (int phi = 0; phi < phi_cnt; phi++)
      {
        CameraSettings camera;
        camera.target = vec3(0, 0, 0);
        camera.up = vec3(0, 1, 0);
        float phi_f = 2*M_PI*phi/(float)phi_cnt;
        float psi_f = 0.5*M_PI*psi/(float)psi_cnt;
        camera.origin = dist*vec3(cos(psi_f)*sin(phi_f), sin(psi_f), cos(psi_f)*cos(phi_f));

        render_reference_image_cup(camera, image_w, image_h, "saves/3d_render/cup_test.png", m);
        images.emplace_back("saves/3d_render/cup_test.png");
        cameras.push_back(camera);
        save_image(images.back(), "saves/3d_render/reference_image_"+std::to_string(psi)+"_"+std::to_string(phi)+".png");
      }
    }

    delete m;

    diff_render_3d_naive(vec3(-5,-5,-5), vec3(5,5,5), voxel_size, cameras, images,
                         image_w, image_h, 25, 100, 1, 100, 0.1, filename);
  }

  void render_saved_voxel_array(std::string filename, glm::ivec3 vox_size, int image_size, float dist,
                                float value_thr = 0)
  {
    VoxelArray<vec4> voxel_array(vec3(-5,-5,-5), vec3(5,5,5), vox_size, vec4(0));
    voxel_array.read_from_binary_file(filename);
    for (int i=0;i<voxel_array.get_total_vox_count();i++)
      if (voxel_array.get_direct(i).w < value_thr)
        voxel_array.set_direct(i, vec4(0,0,0,0));
    int n = 0;
    CameraSettings camera;
    camera.target = vec3(0, 0, 0);
    camera.up = vec3(0, 1, 0);

    float image_w = image_size;
    float image_h = image_size;

    Image test_image(image_w, image_h);
    
    for (float t=0;t<=1;t+=0.02)
    {
      camera.origin = vec3(dist*sin(2*M_PI*t), 0, dist*cos(2*M_PI*t));
      render_3d_scene(voxel_array, camera, test_image, image_w, image_h, 25, 100, 4);
      save_image(test_image, "saves/3d_render/test_image.png");
      save_image(test_image, "saves/3d_render/test_image_"+std::to_string(n)+".png");
      n++;
    }
  }

  void diff_render_naive_test_5(std::string filename, glm::ivec3 voxel_size, int image_size)
  {
    std::vector<CameraSettings> cameras;
    std::vector<Image> images;

    int psi_cnt = 4;
    int phi_cnt = 16;
    float dist = 25;

    float image_w = image_size;
    float image_h = image_size;

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

        render_reference_image_circle(camera, images.back(), image_w, image_h);
        save_image(images.back(), "saves/3d_render/reference_image_"+std::to_string(psi)+"_"+std::to_string(phi)+".png");
      }
    }

    diff_render_3d_naive(vec3(-5,-5,-5), vec3(5,5,5), voxel_size, cameras, images,
                         image_w, image_h, 35, 256, 1, 150, 0.1, filename);
  }

  void software_render_test_3d()
  {
    //CameraSettings camera;
    //camera.target = vec3(0, 0, 0);
    //camera.up = vec3(0, 1, 0);
    //camera.origin = vec3(0,0,25);
    //render_reference_image_cup(camera, 256, 256, "saves/3d_render/cup_test.png");
    //diff_render_naive_test_2("saves/3d_render/res_array_32.bin", glm::ivec3(32,32,32), 150);
    //render_saved_voxel_array("saves/3d_render/res_array_32.bin", glm::ivec3(32,32,32), 256);
    
    diff_render_naive_test_4("saves/3d_render/cup_array_128_2.bin", glm::ivec3(128,128,128), 256);
    render_saved_voxel_array("saves/3d_render/cup_array_128_2.bin", glm::ivec3(128,128,128), 256, 17, 0.5);

    //diff_render_naive_test_5("saves/3d_render/circle_array_32.bin", glm::ivec3(32,32,32), 100);
    //render_saved_voxel_array("saves/3d_render/circle_array_32.bin", glm::ivec3(32,32,32), 256, 25);
  }

  void render_test_3d(VoxelArray <glm::vec4> voxel_array)
  {
    //CameraSettings camera;
    //camera.target = vec3(0, 0, 0);
    //camera.up = vec3(0, 1, 0);
    //camera.origin = vec3(0,0,25);
    //render_reference_image_cup(camera, 256, 256, "saves/3d_render/cup_test.png");
    //diff_render_naive_test_2("saves/3d_render/res_array_32.bin", glm::ivec3(32,32,32), 150);
    //render_saved_voxel_array("saves/3d_render/res_array_32.bin", glm::ivec3(32,32,32), 256);
    
    //diff_render_naive_test_4("saves/3d_render/cup_array_128_2.bin", glm::ivec3(128,128,128), 256);
    //render_saved_voxel_array("saves/3d_render/cup_array_128_2.bin", glm::ivec3(128,128,128), 256, 17, 0.5);

    //diff_render_naive_test_5("saves/3d_render/circle_array_32.bin", glm::ivec3(32,32,32), 100);
    //render_saved_voxel_array("saves/3d_render/circle_array_32.bin", glm::ivec3(32,32,32), 256, 25);
    Image im(512, 512);
    CameraSettings camera = {glm::vec3{-3, 0, 0}, glm::vec3{0, 0, 0}, glm::vec3{0, 1, 0}};
    float dist = 2.5;
    int n = 0;
    for (float t=0;t<=1;t+=0.1)
    {
      camera.origin = vec3(dist*sin(2*M_PI*t), 0, dist*cos(2*M_PI*t));
      render_3d_scene(voxel_array, camera, im, 512, 512, 25, 256, 4);
      save_image(im, "saves/3d_render/test_image_"+std::to_string(n)+".png");
      n++;
    }
  }
};