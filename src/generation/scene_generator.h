#pragma once
#include "grove_generator.h"
#include "core/scene.h"

class SceneGenerator
{
public:
    struct SceneGenerationContext
    {
        Scene *scene;
        Block settings;
        std::map<std::string,TreeTypeData> tree_types;
        GroveGenerationData global_ggd;
    };
    
    SceneGenerator(SceneGenerationContext &_ctx): ctx(_ctx){};
    void create_scene_auto();
private:
    LightVoxelsCube *create_grove_voxels(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                                         AABB &influence_box);
    AABB get_influence_AABB(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                            Heightmap &h);
    void generate_grove();
    SceneGenerationContext &ctx;
};