#include "terrain_renderer.h"
#include "graphics_utils/texture_manager.h"
#include "graphics_utils/modeling.h"

TerrainRenderer::TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step):
        terrain({"terrain_render.vs", "terrain_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
        terrainShadow({"terrain_render.vs", "depth.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
        terrain_tex(textureManager.get("terrain"))
        {
            flat_terrain = new Model();
            Visualizer vis;
            vis.heightmap_to_model(h, flat_terrain, size, glm::vec2(8192, 8192), MIN(step.x, step.y),0);

            float aniso = 0.0f;
            glBindTexture(terrain_tex.type, terrain_tex.texture);
            glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &aniso);
            glTexParameterf(terrain_tex.type, GL_TEXTURE_MAX_ANISOTROPY_EXT, aniso); 
            glTexParameteri(terrain_tex.type, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(terrain_tex.type, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(terrain_tex.type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
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
        shader.texture("terrain",terrain_tex);
        shader.uniform("projection",projection);
        shader.uniform("view",view);
        shader.uniform("model",flat_terrain->model);

        flat_terrain->render();
    }