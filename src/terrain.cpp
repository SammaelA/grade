#include "terrain.h"
#include <cmath>
#include "tinyEngine/utility.h"

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
        TerrainRenderer::TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step):
        terrain({"terrain_render.vs", "terrain_render.fs"}, {"in_Position", "in_Normal", "in_Tex"})
        {
            flat_terrain = new Model();
            int x = (2*size.x/step.x) + 1;
            int y = (2*size.y/step.y) + 1;
            glm::vec2 range = h.get_height_range();
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    int ind = flat_terrain->positions.size()/3;
                    glm::vec3 terr_pos = glm::vec3(pos.x - size.x + step.x*i,0,pos.z - size.y + step.y*j);
                terr_pos.y = h.get_height(terr_pos);
                //logerr("%f",terr_pos.y);
                flat_terrain->positions.push_back(terr_pos.x);
                flat_terrain->positions.push_back(terr_pos.y);
                flat_terrain->positions.push_back(terr_pos.z);

                flat_terrain->colors.push_back((0.0*j)/y);
                //flat_terrain->colors.push_back((terr_pos.y - range.x)/(range.y - range.x));
                flat_terrain->colors.push_back(0.1*length(h.get_grad(terr_pos)));
                flat_terrain->colors.push_back((0.0*i)/x);
                flat_terrain->colors.push_back(0);
                if (i != x - 1 && j != y - 1)
                {
                    flat_terrain->indices.push_back(ind);
                    flat_terrain->indices.push_back(ind + 1);
                    flat_terrain->indices.push_back(ind + y + 1);
                    flat_terrain->indices.push_back(ind);
                    flat_terrain->indices.push_back(ind + y + 1);
                    flat_terrain->indices.push_back(ind + y);
                }
                }
            }
            flat_terrain->update();
        }
    TerrainRenderer::~TerrainRenderer()
    {
        if (flat_terrain)
            delete flat_terrain;
    }
    void TerrainRenderer::render(glm::mat4 prc)
    {
        terrain.use();
        terrain.uniform("projectionCamera",prc);
        terrain.uniform("model",flat_terrain->model);

        flat_terrain->render();
    }