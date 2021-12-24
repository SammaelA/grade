#include "graphics_utils/terrain.h"
#include <cmath>
#include "common_utils/utility.h"
#include "tinyEngine/camera.h"
#include "graphics_utils/texture_manager.h"
#include "common_utils/bbox.h"
#include "third_party/stb_image.h"
#include "tinyEngine/image.h"

    Heightmap::Heightmap(glm::vec3 pos, glm::vec2 size, float cell_size):
    Field_2d(pos, size, cell_size)
    {

    }
    Heightmap::Heightmap(glm::vec3 _pos, int _w, int _h):
    Field_2d(_pos, _w, _h)
    {

    }

    float Heightmap::get_height_simple(glm::vec3 position)
    {
        return get_bilinear(position);
    }

    float Heightmap::get_height(glm::vec3 position)
    {
        position.y = 0;
        if (!data)
            return base_val;
        glm::vec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;
        glm::ivec2 ps = glm::ivec2(rp.x > 0 ? rp.x : rp.x-1, rp.y > 0 ? rp.y : rp.y-1);
        float dx = rp.x - ps.x;
        float dy = rp.y - ps.y;
        int rdx = round(dx);
        float h1 = get(ps.x, ps.y);
        float h2 = get(ps.x + 1, ps.y + 1);
        float h3 = get(ps.x + rdx, ps.y + 1 - rdx);

        glm::vec3 p1 = pos + cell_size*glm::vec3(ps.x, 0, ps.y);
        p1.y = 0;
        glm::vec3 p2 = pos + cell_size*glm::vec3(ps.x + 1, 0, ps.y + 1);
        p2.y = 0;
        glm::vec3 p3 = pos + cell_size*glm::vec3(ps.x + rdx, 0, ps.y + 1 - rdx);
        p3.y = 0;

        glm::vec3 barycent = Barycentric(position, p1, p2, p3);

        return barycent.x*h1 + barycent.y*h2 + barycent.z*h3;
    }
    void Heightmap::random_generate(float base, float min, float max)
    {
      fill_perlin(base,min,max);
    }
    void Heightmap::load_from_image(float base, float min, float max, std::string texture_name)
    {
        base_val = base;
        min_val = base;
        max_val = base;
        int image_w = 0, image_h = 0, channels = 0;
        std::string filename = image::base_img_path + texture_name;
        auto *image_data = stbi_load(filename.c_str(), &image_w, &image_h, &channels, 3);
        if (!data || !image_data || !w || !h || !channels)
            return;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                glm::vec2 rp = glm::vec2(image_w*((float)(i + w)/(2*w + 1)), image_h*((float)(j + h)/(2*h + 1)));
                glm::ivec2 ps = rp;
                #define GET_F(a,b) (image_data[channels * (CLAMP(b,0,image_h-1) * image_w + CLAMP(a,0,image_w-1)) + 0]/255.0)
                float dx = rp.x - ps.x;
                float dy = rp.y - ps.y;
                float height = (dx*GET_F(ps.x, ps.y) + (1 - dx)*GET_F(ps.x + 1, ps.y))*(1 - dy) + 
                               (dx*GET_F(ps.x, ps.y + 1) + (1 - dx)*GET_F(ps.x + 1, ps.y + 1))*dy; 

                height = min + (max - min)*height;
                set(i,j,height);
            }
        }
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