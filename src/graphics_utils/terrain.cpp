#include "graphics_utils/terrain.h"
#include <cmath>
#include "common_utils/utility.h"
#include "tinyEngine/camera.h"
#include "graphics_utils/texture_manager.h"

    Heightmap::Heightmap(glm::vec3 pos, glm::vec2 size, float cell_size):
    Field_2d(pos, size, cell_size)
    {

    }
    Heightmap::Heightmap(glm::vec3 _pos, int _w, int _h):
    Field_2d(_pos, _w, _h)
    {

    }


    float Heightmap::get_height(glm::vec3 position)
    {
      return get_bilinear(position);
    }
    void Heightmap::random_generate(float base, float min, float max)
    {
      fill_perlin(base,min,max);
    }
    glm::vec2 Heightmap::get_grad(glm::vec3 position)
    {
        return get_grad_bilinear(position);
    }

    HeightmapTex::HeightmapTex(Heightmap &heightmap, int w, int h)
    {
        float pres = 1;
        glm::vec3 center = glm::vec3(0,0,0);
        float *data = safe_new<float>(w*h, "heightmap_data");
        for (int i=0;i<w;i++)
        {
            for (int j=-0;j<h;j++)
            {
                glm::vec3 pos = center + glm::vec3((i - w/2)*pres,0,(j- w/2)*pres);
                data[i*h + j] = (1e-3)*heightmap.get_height(pos);
            }
        }
        glGenTextures(1, &hmtex);
        glBindTexture(GL_TEXTURE_2D, hmtex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, w, h, 0, GL_RED, GL_FLOAT, data);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = {base_value, base_value, base_value, base_value};
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

        safe_delete<float>(data, "heightmap_data");
    }
    HeightmapTex::~HeightmapTex()
    {
        glDeleteTextures(1, &hmtex);
    }