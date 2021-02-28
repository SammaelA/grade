#include "grove_generation_utils.h"
#include "distribution.h"

void PlanarShadowsMap::set_occluder(glm::vec3 position, float base_val, float r, float _pow)
{
    glm::ivec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;

    int n_cells = ceil(r/cell_size);

    for (int i=-n_cells;i<n_cells;i++)
    {
        for (int j=-n_cells; j<n_cells; j++)
        {
            set_safe(rp.x+i,rp.y+j,base_val/powf(sqrt(i*i + j*j),_pow));
        }
    }
}
void PlanarShadowsMap::clear()
{
    for (int i =0;i<(2*w + 1)*(2*h + 1);i++)
        data[i] = 0;
}
void HabitabilityMap::create(Heightmap &heightmap, GroveMask &mask)
{
    glm::vec2 minmax = heightmap.get_height_range();
    float delta = 0.5*(minmax.y - minmax.x);
    float base = minmax.x + delta;
    const float grad_mult = 0.05;

    fill_perlin(0,0,1,glm::ivec2(1/grad_mult,10*delta));

        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                glm::vec3 position = glm::vec3(pos.x + cell_size*i, 0, pos.z - cell_size*j);
                float noise = 0.3 + 0.7*get(i,j);
                float grad = CLAMP(grad_mult*length(heightmap.get_grad_bilinear(position)),0,1);
                float height = CLAMP(abs(heightmap.get_height(position) - delta)/delta,0.5,1);

                set(i,j,0.25*mask.get_bilinear(position)*(2*noise + grad + height));
            }
        }
}
void GroveMask::set_round(float r)
{
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                glm::vec3 position = glm::vec3(pos.x + cell_size*i, 0, pos.z - cell_size*j);
                if (length(position) < r) 
                    set(i,j,1);
                else
                    set(i,j,0);
            }
        }
}
double DensityMap::calc_sum()
{
    double sum = 0;
    for (int i =0;i<(2*w + 1)*(2*h + 1);i++)
        sum += data[i];
    return sum;
}
void DensityMap::create(HabitabilityMap &hm, PlanarShadowsMap &psm)
{
    const float shadow_q = 0.5;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                glm::vec3 position = glm::vec3(pos.x + cell_size*i, 0, pos.z - cell_size*j);
                float sh = psm.get_bilinear(position);
                float hab = hm.get_bilinear(position);
                float res = (1-shadow_q)*hab + shadow_q*(1/(1 + sh));
                set(i,j,res);
            }
        }
}
void DensityMap::choose_places_for_seeds(int count, std::vector<Seed> &seeds)
{
    if (count <= 0)
        return;
    double sum = calc_sum();
    int mult = 1;
    while (seeds.size() < count)
    {
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float rnd = urand(0,sum);
                if (mult*count*get(i,j) > rnd && seeds.size() < count)
                {
                    Seed s;
                    s.pos = glm::vec2(pos.x + cell_size*i, pos.z - cell_size*j);
                    s.roots_count = 1;
                    seeds.push_back(s);
                }
            }
        }
        mult *= 2;
    }
}