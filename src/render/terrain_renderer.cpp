#include "terrain_renderer.h"
#include "tinyEngine/engine.h"
#include "tree_utils/tree_modeling.h"

TerrainRenderer::TerrainRenderer(Heightmap &h, float3 pos, float2 size, float2 step):
        terrain_tex1(engine::textureManager->get("terrain1")),
        terrain_tex2(engine::textureManager->get("terrain2")),
        terrain_tex3(engine::textureManager->get("terrain3")),
        perlin(engine::textureManager->get("noise")),
        terrain({"terrain_render.vs", "terrain_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
        terrainShadow({"terrain_render.vs", "depth.fs"}, {"in_Position", "in_Normal", "in_Tex"})
        {
            flat_terrain = new Model();
            visualizer::heightmap_to_model(h, flat_terrain, size, float2(8192, 8192), MIN(step.x, step.y),0);
            flat_terrain->update();
            Texture *texs[3] = {&terrain_tex1, &terrain_tex2, &terrain_tex3};
            for (auto *tex : texs)
            {
                float aniso = 0.0f;
                glBindTexture(tex->type, tex->texture);
                glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &aniso);
                glTexParameterf(tex->type, GL_TEXTURE_MAX_ANISOTROPY_EXT, aniso); 
                glTexParameteri(tex->type, GL_TEXTURE_WRAP_S, GL_REPEAT);
                glTexParameteri(tex->type, GL_TEXTURE_WRAP_T, GL_REPEAT);
                glTexParameteri(tex->type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
            }
        }
    TerrainRenderer::~TerrainRenderer()
    {
        if (flat_terrain)
            delete flat_terrain;
    }
    void TerrainRenderer::render(float4x4 projection, float4x4 view, float4x4 shadow_tr, GLuint shadow_tex, float3 camera_pos,
                                 DirectedLight &light, bool to_shadow, int debug_type, float4 grid_params, float4 debug_tex_scale,
                                 Texture *debug_tex_ptr)
    {
      Shader &shader = to_shadow ? terrainShadow : terrain;
      shader.use();
      shader.texture("grass1", terrain_tex1);
      shader.texture("grass2", terrain_tex2);
      shader.texture("rock", terrain_tex3);
      shader.texture("perlin", perlin);
      shader.uniform("projection", projection);
      shader.uniform("view", view);
      shader.uniform("model", flat_terrain->model);
      shader.uniform("debug_render", debug_type);
      shader.uniform("grid_params", grid_params);
      if (debug_tex_ptr)
        shader.texture("debug_tex", *debug_tex_ptr);
      shader.uniform("debug_tex_scale", debug_tex_scale);

      flat_terrain->render();
    }