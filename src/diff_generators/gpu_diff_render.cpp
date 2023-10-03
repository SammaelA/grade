#include "gpu_diff_render.h"
#include "tinyEngine/engine.h"
#include "graphics_utils/modeling.h"
#include "cities_generator/global.h"
#include <glm/glm.hpp>

struct EdgeGPU
{
  glm::vec3 p0;
  float local_pdf;
  glm::vec3 p1;
  float _pad;
};

DiffRenderGPU::DiffRenderGPU():
render_silhouette({"render_silhouette.vs", "render_silhouette.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
diff_loss("diff_loss.fs"),
diff_loss_sum({"diff_loss_sum.comp"},{}),
edges_get_pdf({"edges_get_pdf.comp"},{}),
edge_sampling({"edge_sampling.comp"},{})
{
  fbo = create_framebuffer();

  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
  
  results_buf = create_buffer();
  pdf_buf = create_buffer();
  edges_buf = create_buffer();
  grad_buf = create_buffer();

  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}

DiffRenderGPU::~DiffRenderGPU()
{
  delete_framebuffer(fbo);
  delete_buffer(results_buf);
  delete_buffer(pdf_buf);
  delete_buffer(edges_buf);
  delete_buffer(grad_buf);
}

void DiffRenderGPU::init_optimization(const std::vector<std::string> &reference_image_dir,
                                      MitsubaInterface::RenderSettings render_settings,
                                      bool save_intermediate_images) 
{
  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  assert(reference_image_dir.size() >= 1);
  if (render_settings.renderStyle != MitsubaInterface::SILHOUETTE)
    logerr("DiffRenderGPU supports only silhouette render style");

  tex_w = render_settings.image_w;
  tex_h = render_settings.image_h;
  save_images = save_intermediate_images;
  int edge_samples = tex_w*tex_h;

  for (auto &dir : reference_image_dir)
  {
    reference_textures.push_back(engine::textureManager->load_unnamed_tex(dir, 1));
    assert(reference_textures.back().get_W() == tex_w && reference_textures.back().get_H() == tex_h);
  }

  tex1 = engine::textureManager->create_texture(tex_w, tex_h, GL_R8);
  tex2 = engine::textureManager->create_texture(tex_w, tex_h, GL_R16F);

  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
}

float DiffRenderGPU::render_and_compare(const dgen::DFModel &model, const std::vector<CameraSettings> &cameras,
                                        const std::vector<float> &scene_params)
{
  assert(cameras.size() == reference_textures.size());
  auto t1 = std::chrono::steady_clock::now();

  pos_grad.resize(model.first.size(), 0);
  cam_grad.resize(scene_params.size(), 0);

  auto t2 = std::chrono::steady_clock::now();

  Model m;
  visualizer::simple_mesh_to_model_332(model.first, &m);

  auto t3 = std::chrono::steady_clock::now();

  m.update();

  auto t4 = std::chrono::steady_clock::now();

  int prev_FBO = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex1.texture, 0);
  glViewport(0, 0, tex_w, tex_h);

  glm::mat4 y_swap = glm::mat4(glm::vec4(1,0,0,0), glm::vec4(0,-1,0,0), glm::vec4(0,0,1,0),glm::vec4(0,0,0,1));

  float tr_z = scene_params[2] * tan(0.5*0.25) / tan(0.5*cameras[0].fov_rad);
  glm::mat4 transform = glm::rotate(glm::rotate(glm::rotate(glm::translate(glm::mat4(1.0f),
                                       glm::vec3(scene_params[0], scene_params[1], tr_z)), 
                                       scene_params[5], glm::vec3(0,0,1)), 
                                       scene_params[4], glm::vec3(0,1,0)), 
                                       scene_params[3], glm::vec3(1,0,0));

  auto t5 = std::chrono::steady_clock::now();
  float total_res = 0;

  for (int i=0;i<cameras.size();i++)
  {
    //render model silhouette
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glm::mat4 projection = glm::perspective(cameras[i].fov_rad, 1.0f, cameras[i].z_near, cameras[i].z_far);
    glm::mat4 view = glm::lookAt(cameras[i].origin, cameras[i].target, cameras[i].up);
    glm::mat4 viewProj = y_swap * projection * view * transform;
    render_silhouette.use();
    render_silhouette.uniform("viewProj", viewProj);
    m.render();

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
      std::string dir = std::string("DiffRenderGPU_")+std::to_string(iter)+std::string("_cam_")+std::to_string(i);
      engine::textureManager->save_png(tex1, dir);
    }

    auto t6 = std::chrono::steady_clock::now();
    #define NUM_TILES 16

    //calculate MSE per tile
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    diff_loss_sum.use();
    diff_loss_sum.uniform("tex_size", glm::vec2(tex_w, tex_h));
    diff_loss_sum.texture("tex_diff", tex2);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, results_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*NUM_TILES*NUM_TILES, nullptr, GL_STREAM_READ);
    glDispatchCompute(1, 1, 1);

    auto t7 = std::chrono::steady_clock::now();

    //put edges in buffer, then reproject them and calculate pdfs
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    std::vector<EdgeGPU> edges(m.get_size());
    for (int j=0;j<edges.size(); j+=3)
    {
      edges[j].p0 = glm::vec3(m.positions[3*j], m.positions[3*j+1], m.positions[3*j+2]);
      edges[j].p1 = glm::vec3(m.positions[3*(j+1)], m.positions[3*(j+1)+1], m.positions[3*(j+1)+2]);

      edges[j+1].p0 = glm::vec3(m.positions[3*(j+1)], m.positions[3*(j+1)+1], m.positions[3*(j+1)+2]);
      edges[j+1].p1 = glm::vec3(m.positions[3*(j+2)], m.positions[3*(j+2)+1], m.positions[3*(j+2)+2]);

      edges[j+2].p0 = glm::vec3(m.positions[3*j], m.positions[3*j+1], m.positions[3*j+2]);
      edges[j+2].p1 = glm::vec3(m.positions[3*(j+2)], m.positions[3*(j+2)+1], m.positions[3*(j+2)+2]);
    }
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, edges_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(EdgeGPU)*edges.size(), edges.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, pdf_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*NUM_TILES*NUM_TILES, nullptr, GL_STREAM_READ);
    edges_get_pdf.use();
    edges_get_pdf.uniform("viewProj", viewProj);
    edges_get_pdf.uniform("triangles_cnt", (int)(edges.size()/3));
    edges_get_pdf.uniform("accumulation_pass", false);
    glDispatchCompute(1, 1, 1);
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    edges_get_pdf.uniform("accumulation_pass", true);
    glDispatchCompute(1, 1, 1);

    //TEST: get edges
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, edges_buf);
    GLvoid* eptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(edges.data(),eptr,sizeof(EdgeGPU)*edges.size());
    for (int j=0;j<edges.size();j++)
    {
      //logerr("P %f %f %f %f %f %f", m.positions[3*j], m.positions[3*j+1], m.positions[3*j+2],
      //                              model.first[8*j], model.first[8*j+1],model.first[8*j+2] );
      //logerr("edge %d (%f %f %f) (%f %f %f) %f", j, edges[j].p0.x, edges[j].p0.y, edges[j].p0.z, edges[j].p1.x, edges[j].p1.y, edges[j].p1.z, edges[j].local_pdf);
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); 

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

    auto t8 = std::chrono::steady_clock::now();

    int dt1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    int dt2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    int dt3 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
    int dt4 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count();
    int dt5 = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
    int dt6 = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
    int dt7 = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();
    int dt_all = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t1).count();
    logerr("DiffRenderGPU took %d (%d + %d + %d + %d + %d + %d + %d) us", dt_all,dt1,dt2,dt3,dt4,dt5,dt6,dt7);
  }

  glMemoryBarrier(GL_ALL_BARRIER_BITS);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);
  iter++;

  return 10*log10(MAX(1e-9f,total_res/cameras.size()));
}

const float *DiffRenderGPU::get_vertex_grad() const
{
  return pos_grad.data();
}

const float *DiffRenderGPU::get_scene_params_grad() const
{
  return cam_grad.data();
}
