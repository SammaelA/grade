#include "world_renderer.h"
#define GLEW_EXPERIMENTAL
#include "tinyEngine/TinyEngine.h"
#include "core/body.h"
#include "cities_generator/global.h"

#define DEL_IT(a) if (a) {delete a;a = nullptr;}
const int HALTON_COUNT = 8;
const float JITTER_SCALE = 2.0;
  glm::vec2 haltonSequence[HALTON_COUNT] = {
    glm::vec2(1.0 / 2.0, 1.0 / 3.0),
    glm::vec2(1.0 / 4.0, 2.0 / 3.0),
    glm::vec2(3.0 / 4.0, 1.0 / 9.0),
    glm::vec2(1.0 / 8.0, 4.0 / 9.0),
    glm::vec2(5.0 / 8.0, 7.0 / 9.0),
    glm::vec2(3.0 / 8.0, 2.0 / 9.0),
    glm::vec2(7.0 / 8.0, 5.0 / 9.0),
    glm::vec2(1.0 / 16.0, 8.0 / 9.0)};
WorldRenderer::~WorldRenderer()
{
  regenerate_shadows = true;
  DEL_IT(hbaoRenderer)
  DEL_IT(cubemap)
  DEL_IT(defferedLight)
  DEL_IT(startScreenShader)
  DEL_IT(defaultShader)
  DEL_IT(debugShader)
  DEL_IT(groveRenderer)
  DEL_IT(heightmapTex)
  DEL_IT(grassRenderer)
  DEL_IT(terrainRenderer)
  DEL_IT(debugVisualizer)
  DEL_IT(grassRenderer2);
  DEL_IT(renderReadback);
}

void WorldRenderer::set_resolution(int w, int h)
{
  if (w == screen_w && h == screen_h)
    return;
  screen_w = w;
  screen_h = h;
  
  projection = glm::perspective(fov, (float)w / h, 1.0f, 3000.0f);
}

void WorldRenderer::init(int _h, int _w, Block &render_settings)
{
  if (inited)
    return;
  inited = true;
  float rt_overscale = render_settings.get_double("rt_overscale",1.5);
  int w = rt_overscale*_w;
  int h = rt_overscale*_h;
  target_w = w;
  target_h = h;

  set_resolution(w, h);

  light.dir = glm::normalize(glm::vec3(-0.2, 0.5, -0));
  light.color = glm::vec3(0.99, 0.9, 0.7);
  light.intensity = 2;
  light.ambient_q = 0.1;
  light.diffuse_q = 0.8;
  light.specular_q = 0.1;
  light.has_shadow_map = true;
  light.shadow_map_size = glm::vec2(4096, 4096);

  startScreenShader = new PostFx("simple_render.fs");
  shadowMap.create(light.shadow_map_size.x, light.shadow_map_size.y);
  defferedTarget.create(w, h);
  defferedTarget.set_clear_color(glm::vec4(0.0, 0.0, 0.0, 0.0));

  hbaoRenderer = new HBAORenderer();
  hbaoRenderer->create(w/2, h/2);
  cubemap = new Cubemap(w, h);
  defferedLight = new PostFx("deffered_light.fs");
  debugShader = new Shader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  defaultShader = new Shader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  simpleInstancingShader = new Shader({"simple_instancing.vs", "simple_instancing.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  simpleInstancingShaderShadow = new Shader({"simple_instancing.vs", "simple_instancing_shadow.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  debugVisualizer = new DebugVisualizer(textureManager.get("wood"), defaultShader);
  renderReadback = new RenderReadback();
  targets[0].create(w,h);
  targets[1].create(w,h);
  taa = new PostFx("taa.fs");

}

    void WorldRenderer::set_heightmap(Heightmap &heightmap)
    {
        bool need_simple_grass = false;
        glm::vec2 precision = glm::vec2(8, 8);
        remove_heightmap();
        terrainRenderer = new TerrainRenderer(heightmap, heightmap.get_pos(), heightmap.get_size(), precision);

        if (need_simple_grass)
        {
            heightmapTex = new HeightmapTex(heightmap, 1024, 1024);
            grassRenderer = new GrassRenderer();
        }
    }
    void WorldRenderer::remove_heightmap()
    {
        DEL_IT(terrainRenderer)
        DEL_IT(heightmapTex)
        DEL_IT(grassRenderer)
        on_scene_changed();
    }

    void WorldRenderer::set_grove(GrovePacked &source, GroveGenerationData &gen_data)
    {
        remove_grove();
        int p_int = render_settings.get_int("grove_renderer_precision", (int)GroveRenderer::MEDIUM);   
        GroveRenderer::Precision pres = p_int <= 0 ? GroveRenderer::LOW : (p_int == 1 ? GroveRenderer::MEDIUM : GroveRenderer::DEBUG);
        bool print_perf = render_settings.get_bool("print_perf",false);
        std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
        if (pres == GroveRenderer::Precision::LOW)
            LODs_dists.back() = -10;
        groveRenderer = new GroveRenderer(&source, &gen_data, 5, LODs_dists, print_perf, pres);
    }
    void WorldRenderer::remove_grove()
    {
        delete groveRenderer;
        on_scene_changed();
    }

    void WorldRenderer::set_voxels_debug(LightVoxelsCube &voxels)
    {
        remove_voxels_debug();
        debugVisualizer->visualize_light_voxels(&voxels);
    }
    void WorldRenderer::remove_voxels_debug()
    {
        on_scene_changed();
    }

    void WorldRenderer::add_body_debug(Body *body)
    {
        remove_body_debug();
        debugVisualizer->add_bodies(body, 1);
    }

    void WorldRenderer::add_aabb_debug(const AABB &box)
    {
      glm::vec3 sz = box.max_pos - box.min_pos;
      Box *box_ptr = new Box(box.min_pos, glm::vec3(sz.x,0,0), glm::vec3(0,sz.y,0), glm::vec3(0,0,sz.z));
      debugVisualizer->add_bodies(box_ptr,1);
      on_scene_changed();
      delete box_ptr;
    }

    void WorldRenderer::remove_body_debug()
    {
        on_scene_changed();
    }

    void WorldRenderer::set_grass(GrassPacked &grass_data)
    {
      grassRenderer2 = new GrassRenderer2(grass_data);
    }
    void WorldRenderer::remove_grass()
    {
      DEL_IT(grassRenderer2);
      on_scene_changed();
    }
    
    void WorldRenderer::add_instanced_models(std::vector<Scene::InstancedModel> &_models)
    {
      models = _models;
      int offset = 0;
      std::vector<glm::mat4> all_matrices;
      for (auto &im : models)
      {
        inst_offsets.push_back(offset);
        offset += im.instances.size();
        all_matrices.insert(all_matrices.end(), im.instances.begin(), im.instances.end());
      }
      glGenBuffers(1, &simple_instances_buffer);
      glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 12, simple_instances_buffer);
      glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::mat4)*offset, all_matrices.data(), GL_STATIC_DRAW);
      
      on_scene_changed();
    }

    void WorldRenderer::remove_all_instanced_models()
    {
      glDeleteBuffers(1, &simple_instances_buffer);
      on_scene_changed();
    }

    void WorldRenderer::on_scene_changed()
    {
        regenerate_shadows = true;
    }

    void WorldRenderer::set_render_mode(int _render_mode)
    {
        render_mode = _render_mode;
    }
    void WorldRenderer::set_forced_LOD(int _forced_LOD)
    {
        forced_LOD = _forced_LOD;
    }
    void WorldRenderer::set_groveRendererDebugParams(GroveRendererDebugParams _groveRendererDebugParams)
    {
        groveRendererDebugParams = _groveRendererDebugParams;
    }
    void WorldRenderer::clear_all()
    {
        remove_heightmap();
        remove_grove();
        remove_voxels_debug();
        remove_body_debug();
        render_mode = -1;
        forced_LOD = -1;
    }

void WorldRenderer::render(float dt, Camera &camera)
{
  frame++;
  float deltaWidth = 1.0 / target_w;
  float deltaHeight = 1.0 / target_h;
  uint index = frame % HALTON_COUNT;
 
  glm::vec2 jitter = glm::vec2(haltonSequence[index].x * deltaWidth, haltonSequence [index].y * deltaHeight);
  
  projectionNoJitter = projection;

  projection[3][0] += jitter.x * JITTER_SCALE;
  projection[3][1] += jitter.y * JITTER_SCALE;

  // SHADOW PASS //
  checkForGlErrors("render pre shadow", true);
  if (regenerate_shadows)
  {
    regenerate_shadows = false;
    shadowMap.use(light);
    glm::mat4 sh_viewproj = shadowMap.get_transform();
    if (!models.empty())
    {
      simpleInstancingShaderShadow->use();
      simpleInstancingShaderShadow->uniform("projection", shadowMap.get_projection());
      simpleInstancingShaderShadow->uniform("view", shadowMap.get_view());
      for (int i=0;i<models.size();i++)
      {
        simpleInstancingShaderShadow->uniform("inst_buf_offset", inst_offsets[i]);

        glBindVertexArray(models[i].model->vao);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, models[i].model->ibo);
        glDrawElementsInstanced(GL_TRIANGLES, models[i].model->SIZE, GL_UNSIGNED_INT, 0, models[i].instances.size());
      }
    }
    if (groveRenderer)
    {
      groveRenderer->render(groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), camera,
                                glm::vec2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
                                light, groveRendererDebugParams, sh_viewproj, 0, true);
    }
    if (terrainRenderer)
    {
      terrainRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(),
                                  0, camera.pos, light, 0, glm::vec4(0,0,1,1), true);
    }
    if (grassRenderer2)
    {
           grassRenderer2->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                camera.pos, *heightmapTex, light, true); 
    }
    else if (grassRenderer)
    {
      grassRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                camera.pos, *heightmapTex, light, true);
    }
    shadowMap.start_trans_pass();
    if (groveRenderer)
    {
      groveRenderer->render(groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), camera,
                                glm::vec2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
                                light, groveRendererDebugParams, sh_viewproj, 0, true);
    }
    if (grassRenderer2)
    {
      grassRenderer2->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                camera.pos, *heightmapTex, light, true);
    }
    else if (grassRenderer)
    {
      grassRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                camera.pos, *heightmapTex, light, true);
    }
    shadowMap.finish_trans_pass();

    shadowMap.blur();
  }
  checkForGlErrors("render shadow", true);
      /*
    debugShader->use();
    debugShader->uniform("projectionCamera", projection * camera.camera());
    if (render_mode == DEBUG_RENDER_MODE)
    {
      debugShader->texture("tex", textureManager.get(debug_tex));
      debugShader->texture("tex_arr", textureManager.get(debug_tex));
      debugShader->uniform("need_tex", true);
      debugShader->uniform("need_arr_tex", false);
      debugShader->uniform("need_coord", false);
      debugShader->uniform("slice", 0);
    }
    else if (render_mode == ARRAY_TEX_DEBUG_RENDER_MODE)
    {
      debugShader->texture("tex", textureManager.get_arr(debug_tex));
      debugShader->texture("tex_arr", textureManager.get_arr(debug_tex));
      debugShader->uniform("need_tex", false);
      debugShader->uniform("need_arr_tex", true);
      debugShader->uniform("need_coord", false);
      debugShader->uniform("slice", debug_layer);
    }
    */
  // FILL GBUFFER //

  defferedTarget.target();
  if (render_mode >= DEBUG_ONLY_RENDER_MODE)
  {
    if (terrainRenderer)
    {
      terrainRenderer->render(projection, camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                              camera.pos, light, 
                              debugInfo.get_bool("render_grid_debug", false) ? 1 : 0,
                              debugInfo.get_vec4("grid_params", glm::vec4(0,0,100,100)));
    }
  }
  checkForGlErrors("render terrain", true);
  if (render_mode == ALL_RENDER_MODE || render_mode == MODELS_ONLY_RENDER_MODE)
  {
    if (!models.empty())
    {
      simpleInstancingShader->use();
      simpleInstancingShader->uniform("projection", projection);
      simpleInstancingShader->uniform("view", camera.camera());
      for (int i=0;i<models.size();i++)
      {
        simpleInstancingShader->uniform("inst_buf_offset", inst_offsets[i]);
        simpleInstancingShader->texture("tex", models[i].tex.texture);

        glBindVertexArray(models[i].model->vao);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, models[i].model->ibo);
        glDrawElementsInstanced(GL_TRIANGLES, models[i].model->SIZE, GL_UNSIGNED_INT, 0, models[i].instances.size());
      }
    }
    checkForGlErrors("render models", true);
    if (grassRenderer2)
    {
      grassRenderer2->render(projection, camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                                camera.pos, *heightmapTex, light);
    }
    else if (grassRenderer)
    {
      grassRenderer->render(projection, camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                                camera.pos, *heightmapTex, light);
    }
    checkForGlErrors("render grass", true);
    if (render_mode != 2 && groveRenderer)
    {
      groveRenderer->render(forced_LOD, projection, camera.camera(), camera,
                                 glm::vec2(Tiny::view.WIDTH, Tiny::view.HEIGHT), light,
                                 groveRendererDebugParams, shadowMap.get_transform(), 0 * shadowMap.getTex());
    }
    checkForGlErrors("render trees", true);
  }
  if (render_mode == ALL_RENDER_MODE || render_mode == DEBUG_ONLY_RENDER_MODE)
  {
    if (debugVisualizer)
      debugVisualizer->render(projection * camera.camera(), render_mode);
  }
  checkForGlErrors("render debug", true);

  cubemap->render(projection, camera.camera(), camera);
  checkForGlErrors("render cubemap", true);

  // READBACK FROM GBUFFER //
  if (readback_required && renderReadback)
  {
    RRD.valid = true;
    glm::vec4 wp = renderReadback->get_world_pos(RRID.cursor_screen_pos*defferedTarget.size(), 
                                                 defferedTarget.get_world_pos(), defferedTarget.get_color());
    RRD.cursor_world_pos_type = wp;
    RRD.cursor_on_geometry = wp.w > 0;
  }
  checkForGlErrors("render readback", true);

  // RESOLVE //

  targets[current_target].target();

  glm::vec3 ads = glm::vec3(light.ambient_q, light.diffuse_q, light.specular_q);
  defferedLight->use();
  defferedLight->get_shader().texture("colorTex", defferedTarget.get_color());
  defferedLight->get_shader().texture("normalsTex", defferedTarget.get_normals());
  defferedLight->get_shader().texture("viewPosTex", defferedTarget.get_view_pos());
  defferedLight->get_shader().texture("worldPosTex", defferedTarget.get_world_pos());
  defferedLight->get_shader().texture("aoTex", hbaoRenderer->get_tex());
  defferedLight->get_shader().texture("shadowMap", shadowMap.getTex());
  defferedLight->get_shader().texture("cubeTex", cubemap->get_tex());
  defferedLight->get_shader().uniform("dir_to_sun", light.dir);
  defferedLight->get_shader().uniform("camera_pos", camera.pos);
  defferedLight->get_shader().uniform("ambient_diffuse_specular", ads);
  defferedLight->get_shader().uniform("light_color", light.color*light.intensity);
  defferedLight->get_shader().uniform("need_shadow", shadowMap.getTex() != 0);
  defferedLight->get_shader().uniform("shadow_mat", shadowMap.get_transform());
  defferedLight->get_shader().uniform("sts_inv", 1.0f / light.shadow_map_size);
  defferedLight->render();

  Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  
  taa->use();
  taa->get_shader().texture("target",targets[current_target].get_tex());
  taa->get_shader().texture("prevTarget",targets[(current_target + 1) % 2].get_tex());
  taa->get_shader().uniform("weight",0.95f);
  taa->render();
  checkForGlErrors("render resolve", true);
  
  current_target = (current_target + 1) % 2;
  projection = projectionNoJitter;
}