#include "generation/grove_generation_utils.h"
#include "generation_task.h"
#include "common_utils/distribution.h"
#include "core/tree.h"
#include "common_utils/body.h"

void PlanarShadowsMap::set_occluder(float3 position, float base_val, float r, float _pow)
{
    int2 rp = to_int2(float2(position.x - pos.x, position.z - pos.z)/cell_size);

    int n_cells = ceil(r/cell_size);

    for (int i=-n_cells;i<n_cells;i++)
    {
        for (int j=-n_cells; j<n_cells; j++)
        {
            add_safe(rp.x+i,rp.y+j,base_val/powf(sqrt(i*i + j*j) + 1,_pow));
        }
    }
}
void PlanarShadowsMap::clear()
{
    for (int i =0;i<(2*w + 1)*(2*h + 1);i++)
        data[i] = 0;
}
void PlanarShadowsMap::add(PlanarShadowsMap &src)
{
    int size = (2*w + 1)*(2*h + 1);
    int src_size = (2*src.w + 1)*(2*src.h + 1);
    if (size == src_size)
    {
        for (int i =0;i<size;i++)
            data[i] += src.data[i];
    }
    else
    {
        logerr("PSM add(..) from shadow map with different size is not implemented");
    }
}
void PlanarShadowsMap::set(PlanarShadowsMap &src)
{
    int size = (2*w + 1)*(2*h + 1);
    int src_size = (2*src.w + 1)*(2*src.h + 1);
    if (size == src_size)
    {
        for (int i =0;i<size;i++)
            data[i] = src.data[i];
    }
    else
    {
        logerr("PSM set(..) from shadow map with different size is not implemented");
    }
}
void HabitabilityMap::create(Heightmap &heightmap, GroveMask &mask)
{
    float2 minmax = heightmap.get_height_range();
    float delta = 0.5*(minmax.y - minmax.x);
    float base = minmax.x + delta;
    const float grad_mult = 0.05;

    fill_perlin(0,0,1,int2(1/grad_mult,10*delta));

        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float3 position = float3(pos.x + cell_size*i, 0, pos.z - cell_size*j);
                float noise = 0.3 + 0.7*get(i,j);
                float grad = CLAMP(grad_mult*length(heightmap.get_grad_bilinear(position)),0,1);
                float height = CLAMP(abs(heightmap.get_height(position) - delta)/delta,0.5,1);
                
                float res = 0.25*mask.get_bilinear(position)*(2*noise + grad + height);
                                       if (abs(res) > 1000)
                {
                    logerr("set wrong %f (%d %d) %f %f %f mask %f",res,i,j, noise, grad, height, mask.get_bilinear(position));
                }
                set(i,j, res);
            }
        }
}
void GroveMask::set_round(float r)
{
    base_val = 0;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float3 position = float3(cell_size*i, 0,cell_size*j);
                if (length(position) < r) 
                    set(i,j,1);
                else
                    set(i,j,0);
            }
        }
}
void GroveMask::set_round_min(float r, float val)
{
    base_val = 0;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float3 position = float3(cell_size*i, 0,cell_size*j);
                if (length(position) < r) 
                    set(i,j,MIN(val, get(i,j)));
                //else
                //    set(i,j,0);
            }
        }
}
void GroveMask::set_square(float x, float z)
{
        base_val = 0;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float3 position = float3(pos.x + cell_size*i, 0, pos.z - cell_size*j);
                if (abs(position.x) < x && abs(position.z) < z) 
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
                float3 position = float3(pos.x + cell_size*i, 0, pos.z + cell_size*j);
                float sh = psm.get_bilinear(position);
                float hab = hm.get_bilinear(position);
                float res = hab*(1/(1 + MAX(sh,0)));
                set(i,j,res);
                
                if (abs(res) > 1000)
                {
                    logerr("set wrong %f (%d %d) %f %f",res,i,j, sh, hab);
                }
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
        float rnd = urand(0,sum);
        for (int i = -w; i <= w; i++)
        {
            if (rnd < 0)
                break;
            for (int j = -h; j <= h; j++)
            {
                if (rnd < 0)
                    break;
                float r = get(i,j);
                if (r >= rnd && seeds.size() < count)
                {
                    Seed s;
                    s.pos = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                    s.roots_count = 1;
                    seeds.push_back(s);
                    rnd = -1;
                }
                else
                {
                    rnd -= r;
                }
            }
        }
    }
}