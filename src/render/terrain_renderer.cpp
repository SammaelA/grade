#include "terrain_renderer.h"
#include "graphics_utils/texture_manager.h"
#include "graphics_utils/modeling.h"

TerrainRenderer::TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step):
        terrain({"terrain_render.vs", "terrain_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
        terrainShadow({"terrain_render.vs", "depth.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
        terrain_tex1(textureManager.get("terrain1")),
        terrain_tex2(textureManager.get("terrain2")),
        terrain_tex3(textureManager.get("terrain3")),
        perlin(textureManager.get("noise"))
        {
            flat_terrain = new Model();
            Visualizer vis;
            vis.heightmap_to_model(h, flat_terrain, size, glm::vec2(8192, 8192), MIN(step.x, step.y),0);
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
    void TerrainRenderer::render(glm::mat4 projection, glm::mat4 view, glm::mat4 shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                                 DirectedLight &light, bool to_shadow)
    {
        Shader &shader = to_shadow ? terrainShadow : terrain;
        shader.use();
        shader.texture("grass1",terrain_tex1);
        shader.texture("grass2",terrain_tex2);
        shader.texture("rock",terrain_tex3);
        shader.texture("perlin",perlin);
        shader.uniform("projection",projection);
        shader.uniform("view",view);
        shader.uniform("model",flat_terrain->model);

        flat_terrain->render();
    }