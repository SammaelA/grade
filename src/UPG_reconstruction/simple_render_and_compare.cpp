#include "simple_render_and_compare.h"
#include "tinyEngine/engine.h"

NonDiffRender::NonDiffRender():
  render_silhouette({"render_silhouette.vs", "render_silhouette.fs"}, {"in_Position"}),
    diff_loss_sum({"diff_loss_sum.comp"},{}),
  diff_loss("diff_loss.fs")
  {
    fbo = create_framebuffer();

    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
    
    results_buf = create_buffer();
    vao = create_vertex_array();
    glBindVertexArray(vao);
    vbo = create_buffer();

    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
  }

  NonDiffRender::~NonDiffRender()
  {
    delete_framebuffer(fbo);
    delete_buffer(results_buf);
    glDisableVertexAttribArray(vao);
    delete_vertex_array(vao);
    delete_buffer(vbo);
  }
void NonDiffRender::init_optimization(const std::vector<Texture> &reference_images, 
                                      Settings render_settings, 
                                      bool save_intermediate_images)
{
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  assert(reference_images.size() >= 1);

  tex_w = render_settings.image_w;
  tex_h = render_settings.image_h;
  save_images = save_intermediate_images;
  reference_textures = reference_images;

  tex1 = engine::textureManager->create_texture(tex_w, tex_h, GL_R8);
  tex2 = engine::textureManager->create_texture(tex_w, tex_h, GL_RGB16F);

  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);  
}

float NonDiffRender::render_and_compare_silhouette(const std::vector<float> &positions, 
                                                   const std::vector<CameraSettings> &cameras)
{
  auto t0 = std::chrono::steady_clock::now();

    assert(cameras.size() == reference_textures.size());
    assert(positions.size()%9 == 0);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(float), positions.data(), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    int prev_FBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    float total_res = 0;

    for (int i=0;i<cameras.size();i++)
    {
      glBindVertexArray(vao);
      glBindBuffer(GL_ARRAY_BUFFER, vbo);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
      //render model silhouette
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex1.texture, 0);
      glViewport(0, 0, tex_w, tex_h);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glm::mat4 viewProj = cameras[i].get_viewProj();
      render_silhouette.use();
      render_silhouette.uniform("viewProj", viewProj);
      glDrawArrays(GL_TRIANGLES, 0, positions.size()/3);

      //create diff_loss texture (0.5*(1 + img - target_img))
      glMemoryBarrier(GL_TEXTURE_UPDATE_BARRIER_BIT);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex2.texture, 0);
      diff_loss.use();
      diff_loss.get_shader().uniform("tex_size", glm::vec2(tex_w, tex_h));
      diff_loss.get_shader().texture("tex", tex1);
      diff_loss.get_shader().texture("tex_reference", reference_textures[i]);
      diff_loss.render();

      //save rendered image if needed
      if (save_images)
      {
        glMemoryBarrier(GL_TEXTURE_UPDATE_BARRIER_BIT);
        std::string dir = std::string("NonDiffRender_")+std::to_string(iteration)+std::string("_cam_")+std::to_string(i);
        engine::textureManager->save_png(tex1, dir);
      }

      #define NUM_TILES 16

      //calculate MSE per tile
      glMemoryBarrier(GL_ALL_BARRIER_BITS);
      diff_loss_sum.use();
      diff_loss_sum.uniform("tex_size", glm::vec2(tex_w, tex_h));
      diff_loss_sum.texture("tex_diff", tex2);
      glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, results_buf);
      glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*NUM_TILES*NUM_TILES, nullptr, GL_STREAM_READ);
      glDispatchCompute(1, 1, 1);

      //accumulate per-tile MSE to one value
      glMemoryBarrier(GL_ALL_BARRIER_BITS);
      float res = 0;
      float results[NUM_TILES*NUM_TILES];
      glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
      GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
      memcpy(results,ptr,sizeof(float)*NUM_TILES*NUM_TILES);
      for (int j=0;j<NUM_TILES*NUM_TILES;j++)
        res += results[j];
      total_res += res;
      glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);  
    }

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

    auto t4 = std::chrono::steady_clock::now();
    int t_all = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t0).count();
    //logerr("Non-diff. render took = %.1f ms", 1e-3*t_all);

    iteration++;
    static float average_ms = 0.0;
    average_ms += 1e-3*t_all;
    //logerr("average val %f", average_ms/(iteration));

    return total_res/cameras.size();
}

const float *NonDiffRender::get_vertex_grad() const
{
  logerr("Trying to get gradint from NON-DIFFERENTIABLE renderer!");
  return nullptr;
}