#pragma once
#include "grove_renderer.h"
#include "shadow.h"
#include "visualizer.h"
#include "grass_renderer.h"
#include "ambient_occlusion.h"
#include "tinyEngine/cubemap.h"
#include "tinyEngine/deffered_target.h"
#include "terrain_renderer.h"

class WorldRenderer
{
public:
    WorldRenderer() {}
    ~WorldRenderer();
    void init(int h, int w, Block &render_settings);
    void render(float dt, Camera &camera);

    void set_heightmap(Heightmap &heightmap);
    void remove_heightmap();

    void set_grove(GrovePacked &source, GroveGenerationData &gen_data);
    void remove_grove();

    void set_voxels_debug(LightVoxelsCube &voxels);
    void remove_voxels_debug();

    void add_body_debug(Body *body);
    void remove_body_debug();

    void set_render_mode(int _render_mode);
    void set_forced_LOD(int _forced_LOD); //forced_LOD = -1 means no forced LOD
    void set_groveRendererDebugParams(GroveRendererDebugParams _groveRendererDebugParams);
    void clear_all();
    
private:
    void on_scene_changed();

    Block render_settings;
    int render_mode = -1;
    int forced_LOD = -1;
    GroveRendererDebugParams groveRendererDebugParams;
    bool regenerate_shadows = true;
    glm::mat4 projection;

    ShadowMap shadowMap;
    DefferedTarget defferedTarget;
    HBAORenderer *hbaoRenderer = nullptr;
    Cubemap *cubemap = nullptr;
    PostFx *defferedLight = nullptr;
    PostFx *startScreenShader = nullptr;
    Shader *defaultShader = nullptr;
    Shader *debugShader = nullptr;

    GroveRenderer *groveRenderer = nullptr;
    HeightmapTex *heightmapTex = nullptr;
    GrassRenderer *grassRenderer = nullptr;
    TerrainRenderer *terrainRenderer = nullptr;
    DebugVisualizer *debugVisualizer = nullptr;
    DirectedLight light;

    const int DEBUG_RENDER_MODE = -2;
    const int ARRAY_TEX_DEBUG_RENDER_MODE = -3;
    const int MAX_RENDER_MODE = 2;
    const float fov = glm::radians(90.0f);
};