#include "terrain.h"
#include <cmath>

/* Function to linearly interpolate between a0 and a1
 * Weight w should be in the range [0.0, 1.0]
 */
float interpolate(float a0, float a1, float w) {
    /* // You may want clamping by inserting:
     * if (0.0 > w) return a0;
     * if (1.0 < w) return a1;
     */
    return (a1 - a0) * w + a0;
    /* // Use this cubic interpolation [[Smoothstep]] instead, for a smooth appearance:
     * return (a1 - a0) * (3.0 - w * 2.0) * w * w + a0;
     *
     * // Use [[Smootherstep]] for an even smoother result with a second derivative equal to zero on boundaries:
     * return (a1 - a0) * ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) + a0;
     */
}

typedef struct {
    float x, y;
} vector2;

/* Create random direction vector
 */
vector2 randomGradient(int ix, int iy) {
    // Random float. No precomputed gradients mean this works for any number of grid coordinates
    float random = 2920.f * sin(ix * 21942.f + iy * 171324.f + 8912.f) * cos(ix * 23157.f * iy * 217832.f + 9758.f);
    return (vector2) { .x = cos(random), .y = sin(random) };
}

// Computes the dot product of the distance and gradient vectors.
float dotGridGradient(int ix, int iy, float x, float y) {
    // Get gradient from integer coordinates
    vector2 gradient = randomGradient(ix, iy);

    // Compute the distance vector
    float dx = x - (float)ix;
    float dy = y - (float)iy;

    // Compute the dot-product
    return (dx*gradient.x + dy*gradient.y);
}

// Compute Perlin noise at coordinates x, y
float perlin(float x, float y) {
    // Determine grid cell coordinates
    int x0 = (int)x;
    int x1 = x0 + 1;
    int y0 = (int)y;
    int y1 = y0 + 1;

    // Determine interpolation weights
    // Could also use higher order polynomial/s-curve here
    float sx = x - (float)x0;
    float sy = y - (float)y0;

    // Interpolate between grid point gradients
    float n0, n1, ix0, ix1, value;

    n0 = dotGridGradient(x0, y0, x, y);
    n1 = dotGridGradient(x1, y0, x, y);
    ix0 = interpolate(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);
    ix1 = interpolate(n0, n1, sx);

    value = interpolate(ix0, ix1, sy);
    return value;
}
    Heightmap::Heightmap(glm::vec3 pos, glm::vec2 size, float cell_size):
    Heightmap(pos, ceil(size.x/cell_size), ceil(size.y/cell_size))
    {
        this->cell_size = cell_size;
    }
    Heightmap::Heightmap(glm::vec3 _pos, int _w, int _h)
    {
        w = _w;
        h = _h;
        pos = _pos;

        data = new float[(2*w + 1)*(2*h + 1)];
    }
    Heightmap::~Heightmap()
    {
        if (data)
            delete[] data;
    }
    float Heightmap::get(int x, int y)
    {
        if (x >= -w && x <= w && y >= -h && y <= h)
        {
            return data[(2*w + 1)*(y + h) + x + w];
        }
        else
            return base_height;
    }
    float Heightmap::get_height(glm::vec3 position)
    {
        if (!data)
            return base_height;
        glm::vec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;
        glm::ivec2 ps = rp;
        float dx = rp.x - ps.x;
        float dy = rp.y - ps.y;
        return (dx*get(ps.x, ps.y) + (1 - dx)*get(ps.x + 1, ps.y))*(1 - dy) + 
               (dx*get(ps.x, ps.y + 1) + (1 - dx)*get(ps.x + 1, ps.y + 1))*dy;
        
    }
    void Heightmap::random_generate(float base, float min, float max)
    {
        base_height = base;
        if (!data)
            return;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                data[(2*w + 1)*(j + h) + i + w] = min + 
                    0.5*(max - min)*(1 + perlin(10*(float)(i + w)/(2*w + 1), 10*(float)(j + h)/(2*h + 1)));
                //logerr("%d %f",(2*w + 1)*(j + h) + i + w,data[(2*w + 1)*(j + h) + i + w]);
            }
        }
    }
        TerrainRenderer::TerrainRenderer(Heightmap &h, glm::vec3 pos, glm::vec2 size, glm::vec2 step):
        terrain({"terrain_render.vs", "terrain_render.fs"}, {"in_Position", "in_Normal", "in_Tex"})
        {
            flat_terrain = new Model();
            int x = (2*size.x/step.x) + 1;
            int y = (2*size.y/step.y) + 1;
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

                flat_terrain->colors.push_back((0.5*j)/y);
                flat_terrain->colors.push_back(0.02*terr_pos.y);
                flat_terrain->colors.push_back((0.5*i)/x);
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