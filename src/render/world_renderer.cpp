#include "world_renderer.h"
#define GLEW_EXPERIMENTAL
#include "common_utils/body.h"
#include "tree_utils/tree_modeling.h"
#include "tinyEngine/engine.h"

#define DEL_IT(a) if (a) {delete a;a = nullptr;}
const int HALTON_COUNT = 8;
const float JITTER_SCALE = 2.0;
  float2 haltonSequence[HALTON_COUNT] = {
    float2(1.0 / 2.0, 1.0 / 3.0),
    float2(1.0 / 4.0, 2.0 / 3.0),
    float2(3.0 / 4.0, 1.0 / 9.0),
    float2(1.0 / 8.0, 4.0 / 9.0),
    float2(5.0 / 8.0, 7.0 / 9.0),
    float2(3.0 / 8.0, 2.0 / 9.0),
    float2(7.0 / 8.0, 5.0 / 9.0),
    float2(1.0 / 16.0, 8.0 / 9.0)};
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
  DEL_IT(grassRenderer2)
  DEL_IT(renderReadback)
  DEL_IT(simpleInstancingShader)
  DEL_IT(simpleInstancingShaderShadow)
  DEL_IT(taa)
}

void WorldRenderer::set_resolution(int w, int h)
{
  if (w == screen_w && h == screen_h)
    return;
  screen_w = w;
  screen_h = h;
  
  projection = LiteMath::perspective(fov, (float)w / h, 1.0f, 3000.0f);
  
  defferedTarget = DefferedTarget();
  defferedTarget.create(render_settings.RT_overscale *w, render_settings.RT_overscale*h);
  defferedTarget.set_clear_color(float4(0.0, 0.0, 0.0, 0.0));
}

void WorldRenderer::init(int _h, int _w, Block &render_preset)
{
  if (inited)
    return;

  std::string preset = render_preset.get_string("preset","high_quality");
  if (preset == "high_quality")
  {
    render_settings.RT_overscale = 2.0;
    render_settings.shadow_quality = RenderSettings::HIGH;
  }
  else if (preset == "medium_quality")
  {
    render_settings.RT_overscale = 1.5;
    render_settings.shadow_quality = RenderSettings::LOW;
  }
  else if (preset == "low_quality")
  {
    render_settings.RT_overscale = 1.0;
    render_settings.shadow_quality = RenderSettings::NONE;
  }
  else
  {
    logerr("Unknown render preset %s", preset.c_str());
  }

  inited = true;
  int w = _w;
  int h = _h;
  target_w = w;
  target_h = h;

  set_resolution(w, h);

  light.dir = normalize(float3(-0.5, 0.5, -0.15));
  light.color = float3(0.99, 0.9, 0.7);
  light.intensity = 2;
  light.ambient_q = 0.1;
  light.diffuse_q = 0.8;
  light.specular_q = 0.1;
  light.has_shadow_map = true;
  light.shadow_map_size = render_settings.shadow_quality == RenderSettings::HIGH ? float2(4096, 4096) : float2(2048, 2048);

  startScreenShader = new PostFx("simple_render.fs");
  shadowMap.create(light.shadow_map_size.x, light.shadow_map_size.y);

  hbaoRenderer = new HBAORenderer();
  hbaoRenderer->create(w/2, h/2);
  cubemap = new Cubemap(w, h);
  defaultShader = new Shader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  debugShader = new Shader({"simple_debug.vs", "simple_debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  simpleInstancingShader = new Shader({"simple_instancing.vs", "simple_instancing.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  simpleInstancingShaderShadow = new Shader({"simple_instancing.vs", "simple_instancing_shadow.fs"}, {"in_Position", "in_Normal", "in_Tex"});

  if (render_settings.shadow_quality == RenderSettings::HIGH)
    defferedLight = new PostFx("deffered_light.fs");
  else if (render_settings.shadow_quality == RenderSettings::LOW)
    defferedLight = new PostFx("deffered_light_simple_shadows.fs");
  else 
    defferedLight = new PostFx("deffered_light_no_shadows.fs");

  renderReadback = new RenderReadback();
  targets[0] = RenderTarget();
  targets[0].create(w,h);

  targets[1] = RenderTarget();
  targets[1].create(w,h);
  taa = new PostFx("taa.fs");

}

    void WorldRenderer::set_heightmap(Heightmap &heightmap)
    {
        bool need_simple_grass = false;
        float2 precision = float2(8, 8);
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

    void WorldRenderer::set_grove(const GrovePacked &source, AABB2D scene_bbox, const std::vector<TreeTypeData> &types)
    {
        remove_grove();
        int p_int = (int)GroveRenderer::MEDIUM;  
        GroveRenderer::Precision pres = p_int <= 0 ? GroveRenderer::LOW : (p_int == 1 ? GroveRenderer::MEDIUM : GroveRenderer::DEBUG);
        bool print_perf = false;
        std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
        if (pres == GroveRenderer::Precision::LOW)
            LODs_dists.back() = -10;
        groveRenderer = new GroveRenderer(&source, scene_bbox, types, 5, LODs_dists, print_perf, pres);
    }
    void WorldRenderer::remove_grove()
    {
      delete groveRenderer;
      on_scene_changed();
    }

    void WorldRenderer::set_grass(const GrassPacked &grass_data)
    {
      remove_grass();
      grassRenderer2 = new GrassRenderer2(grass_data);
    }
    void WorldRenderer::remove_grass()
    {
      DEL_IT(grassRenderer2);
      on_scene_changed();
    }
    
    void WorldRenderer::add_instanced_models(const std::vector<Scene::InstancedModel> &_models)
    {
      models = _models;
      int offset = 0;
      std::vector<float4x4> all_matrices;
      for (auto &im : models)
      {
        inst_offsets.push_back(offset);
        offset += im.instances.size();
        for (auto &in : im.instances)
          all_matrices.push_back(in);
      }
      simple_instances_buffer = create_buffer();
      glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 12, simple_instances_buffer);
      glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float4x4)*all_matrices.size(), all_matrices.data(), GL_STATIC_DRAW);
      
      on_scene_changed();
    }

    void WorldRenderer::remove_all_instanced_models()
    {
      inst_offsets.clear();
      if (simple_instances_buffer)
      {
        delete_buffer(simple_instances_buffer);
        simple_instances_buffer = 0;
      }
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
        render_mode = -1;
        forced_LOD = -1;
    }

void WorldRenderer::render(float dt, Camera &camera)
{
  frame++;
  float deltaWidth = 1.0 / target_w;
  float deltaHeight = 1.0 / target_h;
  uint index = frame % HALTON_COUNT;
 
  float2 jitter = float2(haltonSequence[index].x * deltaWidth, haltonSequence [index].y * deltaHeight);
  
  projectionNoJitter = projection;

  projection(0,3) += jitter.x * render_settings.RT_overscale;
  projection(1,3) += jitter.y * render_settings.RT_overscale;

  // SHADOW PASS //
  checkGLErrors("render pre shadow", true);
  if (regenerate_shadows && render_settings.shadow_quality != RenderSettings::NONE)
  {
    regenerate_shadows = false;
    shadowMap.use(light);
    float4x4 sh_viewproj = shadowMap.get_transform();
    if (!models.empty())
    {
      simpleInstancingShaderShadow->use();
      simpleInstancingShaderShadow->uniform("projection", shadowMap.get_projection());
      simpleInstancingShaderShadow->uniform("view", shadowMap.get_view());
      for (int i=0;i<models.size();i++)
      {
        simpleInstancingShaderShadow->uniform("inst_buf_offset", inst_offsets[i]);

        for (int j=0;j<models[i].model.models.size();j++)
        {
          Model &pm = *(models[i].model.models[j]);
          Material &mat = models[i].model.materials[j];

          glBindVertexArray(pm.vao);
          glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pm.ibo);
          glDrawElementsInstanced(GL_TRIANGLES, pm.SIZE, GL_UNSIGNED_INT, 0, models[i].instances.size());
        }
      }
    }
    if (terrainRenderer)
    {
      //terrainRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(),
      //                            0, camera.pos, light, true);
    }
    shadowMap.start_trans_pass();
    if (groveRenderer)
    {
      groveRenderer->render(groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), camera,
                                float2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
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
  checkGLErrors("render shadow", true);

  // FILL GBUFFER //

  defferedTarget.target();
  if (render_mode >= DEBUG_ONLY_RENDER_MODE)
  {
    if (terrainRenderer)
    {
      int debug_type = debugInfo.get_bool("render_grid_debug", false) ? 1 : 0;
      Texture *debug_tex = nullptr;
      if (debugInfo.get_bool("render_grove_mask_debug", false))
      {
        debug_type += 2;
        auto it = debug_textures.find("grove_mask");
        if (it != debug_textures.end())
          debug_tex = &(it->second);
      }
      else if (debugInfo.get_bool("render_biome_mask_debug", false))
      {
        debug_type += 2;
        auto it = debug_textures.find("biome_mask");
        if (it != debug_textures.end())
          debug_tex = &(it->second);
      }
      terrainRenderer->render(projection, camera.camera(), shadowMap.get_transform(), 0,
                              camera.pos, light, false, debug_type,
                              debugInfo.get_vec4("grid_params", float4(0,0,100,100)),
                              debugInfo.get_vec4("debug_tex_scale", float4(0,0,0.01,0.01)),
                              debug_tex);
    }
  }
  checkGLErrors("render terrain", true);
  if (render_mode == ALL_RENDER_MODE || render_mode == MODELS_ONLY_RENDER_MODE)
  {
    if (!models.empty())
    {
      simpleInstancingShader->use();
      simpleInstancingShader->uniform("projection", projection);
      simpleInstancingShader->uniform("view", camera.camera());
      simpleInstancingShader->uniform("type",(int)(RenderPixelTypes::PIXEL_TYPE_MODELS));
      for (int i=0;i<models.size();i++)
      {
        simpleInstancingShader->uniform("inst_buf_offset", inst_offsets[i]);
        for (int j=0;j<models[i].model.models.size();j++)
        {
          Model &pm = *(models[i].model.models[j]);
          Material &mat = models[i].model.materials[j];

          glActiveTexture(GL_TEXTURE0);
          glBindTexture(mat.map_Ka.type, mat.map_Ka.texture);
          simpleInstancingShader->uniform("tex", 0);

          glBindVertexArray(pm.vao);
          glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pm.ibo);
          glDrawElementsInstanced(GL_TRIANGLES, pm.SIZE, GL_UNSIGNED_INT, 0, models[i].instances.size());
        }
      }
    }
    checkGLErrors("render models", true);
    if (grassRenderer2)
    {
      grassRenderer2->render(projection, camera.camera(), shadowMap.get_transform(), 0,
                                camera.pos, *heightmapTex, light);
    }
    else if (grassRenderer)
    {
      grassRenderer->render(projection, camera.camera(), shadowMap.get_transform(), 0,
                                camera.pos, *heightmapTex, light);
    }
    checkGLErrors("render grass", true);
    if (render_mode != 2 && groveRenderer)
    {
      groveRenderer->render(groveRenderer->get_max_LOD(), projection, camera.camera(), camera,
                                 float2(engine::view->WIDTH, engine::view->HEIGHT), light,
                                 groveRendererDebugParams, shadowMap.get_transform(), 0);
    }
    checkGLErrors("render trees", true);
  }
  debugShader->use();
  debugShader->uniform("projection", projection);
  debugShader->uniform("view", camera.camera());
  for (int i=0;i<debug_models.size();i++)
  {
    debugShader->uniform("model", float4x4());
    if (debug_models[i].apply_light)
      debugShader->uniform("type",(int)(RenderPixelTypes::PIXEL_TYPE_DEBUG_LIGHT));
    else
      debugShader->uniform("type",(int)(RenderPixelTypes::PIXEL_TYPE_DEBUG_NO_LIGHT));
    if (debug_models[i].m)
      debug_models[i].m->render();
  }
  checkGLErrors("render debug", true);

  cubemap->render(projection, camera.camera(), camera);
  checkGLErrors("render cubemap", true);

  // READBACK FROM GBUFFER //
  if (readback_required && renderReadback)
  {
    RRD.valid = true;
    float4 wp = renderReadback->get_world_pos(RRID.cursor_screen_pos*defferedTarget.size(), 
                                                 defferedTarget.get_world_pos(), defferedTarget.get_color());
    RRD.cursor_world_pos_type = wp;
    RRD.cursor_on_geometry = wp.w > 0;
  }
  checkGLErrors("render readback", true);

  // RESOLVE //

  targets[current_target].target();

  float3 ads = float3(light.ambient_q, light.diffuse_q, light.specular_q);
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
  defferedLight->get_shader().uniform("need_shadow", shadowMap.getTex().texture != 0);
  defferedLight->get_shader().uniform("shadow_mat", shadowMap.get_transform());
  defferedLight->get_shader().uniform("sts_inv", 1.0f / light.shadow_map_size);
  defferedLight->render();

  engine::view->target(float3(0.6, 0.7, 1));
  
  taa->use();
  taa->get_shader().texture("target",targets[current_target].get_tex());
  taa->get_shader().texture("prevTarget",targets[(current_target + 1) % 2].get_tex());
  taa->get_shader().uniform("weight",0.95f);
  taa->render();
  checkGLErrors("render resolve", true);
  
  current_target = (current_target + 1) % 2;
  projection = projectionNoJitter;
}