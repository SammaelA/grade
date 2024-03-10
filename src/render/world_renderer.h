#pragma once
#include "grove_renderer.h"
#include "shadow.h"
#include "grass_renderer.h"
#include "ambient_occlusion.h"
#include "tinyEngine/cubemap.h"
#include "tinyEngine/deffered_target.h"
#include "terrain_renderer.h"
#include "tinyEngine/render_target.h"
#include "core/scene.h"
#include "render_readback.h"

class RenderCmdExecutor;
class WorldRenderer
{
public:
    friend class RenderCmdExecutor;
    WorldRenderer() {}
    ~WorldRenderer();
    void init(int h, int w, Block &render_settings);
    void render(float dt, Camera &camera);
    void set_readback(RenderReadbackInputData &rrid)
    {
        RRID = rrid;
        readback_required = true;
    }
    RenderReadbackData get_readback() {return RRD;}
    void set_heightmap(Heightmap &heightmap);
    void remove_heightmap();

    void set_grove(const GrovePacked &source, AABB2D scene_bbox, const std::vector<TreeTypeData> &types);
    void remove_grove();
    
    void add_instanced_models(const std::vector<Scene::InstancedModel> &models);
    void remove_all_instanced_models();

    void set_grass(const GrassPacked &grass_data);
    void remove_grass();

    void set_render_mode(int _render_mode);
    void set_forced_LOD(int _forced_LOD); //forced_LOD = -1 means no forced LOD
    void set_groveRendererDebugParams(GroveRendererDebugParams _groveRendererDebugParams);
    void clear_all();
    
    void set_resolution(int w, int h);
    bool is_inited() {return inited;}
private:
    void on_scene_changed();

    struct DebugModel
    {
      Model *m;
      uint64_t id;
      bool apply_light = false;
    };

    struct RenderSettings
    {
      float RT_overscale = 2.0;
      enum {
        NONE, LOW, HIGH
      } shadow_quality = HIGH;
    };

    RenderSettings render_settings;
    int render_mode = -1;
    int forced_LOD = -1;
    GroveRendererDebugParams groveRendererDebugParams;
    bool regenerate_shadows = true;
    glm::mat4 projection, projectionNoJitter;

    ShadowMap shadowMap;
    DefferedTarget defferedTarget;
    HBAORenderer *hbaoRenderer = nullptr;
    Cubemap *cubemap = nullptr;
    PostFx *defferedLight = nullptr;
    PostFx *startScreenShader = nullptr;
    Shader *defaultShader = nullptr;
    Shader *debugShader = nullptr;
    Shader *simpleInstancingShader = nullptr;
    Shader *simpleInstancingShaderShadow = nullptr;
    PostFx *taa = nullptr;

    GroveRenderer *groveRenderer = nullptr;
    HeightmapTex *heightmapTex = nullptr;
    GrassRenderer *grassRenderer = nullptr;
    GrassRenderer2 *grassRenderer2 = nullptr;
    TerrainRenderer *terrainRenderer = nullptr;
    RenderReadback *renderReadback = nullptr;
    DirectedLight light;
    RenderTarget targets[2];
    
    GLuint simple_instances_buffer = 0;
    std::vector<Scene::InstancedModel> models;
    std::vector<int> inst_offsets;//same size ^

    std::vector<DebugModel> debug_models;
    std::map<std::string, Texture> debug_textures;
    int current_target = 0;
    int target_w = 0, target_h = 0;
    int screen_w = 0, screen_h = 0;
    bool inited = false;
    unsigned long frame = 0;
    const int MODELS_ONLY_RENDER_MODE = 1;
    const int ALL_RENDER_MODE = 0;
    const int DEBUG_ONLY_RENDER_MODE = -1;
    const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
    const int MAX_RENDER_MODE = 2;
    const float fov = LiteMath::to_radians(90.0f);
    
    RenderReadbackInputData RRID;
    RenderReadbackData RRD;
    bool readback_required = false;
    Block debugInfo;
};