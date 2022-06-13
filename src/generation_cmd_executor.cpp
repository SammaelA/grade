#include "cmd_executors.h"
#include "common_utils/utility.h"

namespace scene_gen
{
    void create_heightmap(Scene *scene, Block &settings)
    {
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        scene->heightmap = new Heightmap(settings.get_vec3("pos", glm::vec3(0, 0, 0)),
                                         settings.get_vec2("size", glm::vec2(1000, 1000)),
                                         settings.get_double("pixel_size", 10));
        bool load_from_image = settings.get_bool("load_from_image", true);
        float base = settings.get_double("base_height", 0);
        float mn = settings.get_double("min_height", 0);
        float mx = settings.get_double("max_height", 100);
        if (load_from_image)
        {
            std::string hmap_tex_name = settings.get_string("tex_name", "heightmap1.jpg");
            scene->heightmap->load_from_image(base, mn, mx, hmap_tex_name);
        }
        else
        {
            scene->heightmap->random_generate(base, mn, mx);
        }
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        float ms = 1e-4 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        debug("created heightmap. Took %.2f ms\n", ms);
    }
}
void GenerationCmdExecutor::execute(int max_cmd_count)
{
    int cmd_left = max_cmd_count;
    while (!genCmdBuffer.empty() && cmd_left != 0)
    {
        auto cmd = genCmdBuffer.pop();
        switch (cmd.type)
        {
        case GC_GEN_HMAP:
            if (genCtx.scene->heightmap)
                delete genCtx.scene->heightmap;
            scene_gen::create_heightmap(genCtx.scene, cmd.args);
            break;

        default:
            logerr("GenerationCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
            break;
        }
        cmd_left--;
    }
}