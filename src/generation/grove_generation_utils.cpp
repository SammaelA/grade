#include "generation/grove_generation_utils.h"
#include "common_utils/distribution.h"
#include "core/tree.h"
#include "core/body.h"
#include "tree_generators/generated_tree.h"

void PlanarShadowsMap::set_occluder(glm::vec3 position, float base_val, float r, float _pow)
{
    glm::ivec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;

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
                glm::vec3 position = glm::vec3(cell_size*i, 0,cell_size*j);
                if (length(position) < r) 
                    set(i,j,1);
                else
                    set(i,j,0);
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
                glm::vec3 position = glm::vec3(pos.x + cell_size*i, 0, pos.z - cell_size*j);
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
                glm::vec3 position = glm::vec3(pos.x + cell_size*i, 0, pos.z + cell_size*j);
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
                    s.pos = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
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
Seeder::Seeder(GroveGenerationData &ggd, float cell_size, const Heightmap *h):
Seeder(ggd.pos, ggd.size,cell_size, h)
{

}
Seeder::Seeder(glm::vec3 pos, glm::vec3 size, float cell_size, const Heightmap *h):
Countable(6),
mask(pos,glm::vec2(size.x, size.z),cell_size),
hm(pos,glm::vec2(size.x, size.z),cell_size),
psm(pos,glm::vec2(size.x, size.z),cell_size),
const_psm(pos,glm::vec2(size.x, size.z),cell_size),
dsm(pos,glm::vec2(size.x, size.z),cell_size)
{
    heightmap = (Heightmap *)h;
    //mask.set_round(0.85*ggd.size.x);
    mask.set_square(0.85*size.x,0.85*size.z);
    hm.create(*heightmap,mask);
    psm.set(const_psm);
    dsm.create(hm, psm);
}
void Seeder::add_body(Body *b, float opacity, bool solid)
{
    const_psm.add_body(b,opacity,solid);
}
void PlanarShadowsMap::add_body(Body *b, float opacity, bool solid)
{
    if (!b)
        return;
    const int iters = 100;
    BBox bb = b->get_Bbox();

        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float y_start = bb.position.y;
                float step = bb.sizes.y/iters;
                glm::vec3 ps = glm::vec3(pos.x + cell_size*i, y_start, pos.z + cell_size*j);
                for (int k = 0;k<iters;k++)
                {
                    if (b->in_body(ps))
                        Field_2d::add(i,j,opacity);
                }

            }
        }
}
int Seeder::joints_count(Branch *b)
{
    if (!b)
        return 0;
    int cnt = b->joints.size();
    for (Joint &j : b->joints)
    {
        for (Branch *br : j.childBranches)
            cnt += joints_count(br);
    }
    return cnt;
}
void Seeder::recalculate_planar_shadows(Branch *b, PlanarShadowsMap &psm, int level)
{
    if (!b || b->level > level)
        return;
    if (b->level < level)
    {
        for (Joint &j : b->joints)
        {
            for (Branch *br : j.childBranches)
                recalculate_planar_shadows(br,psm,level);
        }
    }
    else //b->level == level
    {
        if (b->joints.size() < 2)
            return;

        int cnt = joints_count(b);
        glm::vec3 pos = 0.5f*(b->joints.back().pos + b->joints.front().pos);
        float height = MAX(0,pos.y - heightmap->get_height(pos));
        float size = 0.5f*(length(b->joints.back().pos - b->joints.front().pos) + height);
        psm.set_occluder(pos,cnt,size,1);

    }
}
void Seeder::recalcuate_shadows(Tree *trees, int count)
{
    psm.set(const_psm);
    for (int i = 0; i < count; i++)
    {
        if (trees[i].valid)
            recalculate_planar_shadows(trees[i].root, psm, 1);
    }
   
    dsm.create(hm, psm);
}
void Seeder::choose_places_for_seeds(int count, std::vector<Seed> &seeds)
{
    dsm.choose_places_for_seeds(count, seeds);
}
void Seeder::add_tree_shadow(Tree &t)
{
    recalculate_planar_shadows(t.root, psm, 1);
    dsm.create(hm, psm);
}

void Seeder::recalcuate_shadows(mygen::Tree *trees, int count)
{
    psm.set(const_psm);
    for (int i = 0; i < count; i++)
    {
        recalculate_planar_shadows(trees[i].root, psm, 1);
    }
   
    dsm.create(hm, psm);
}
void Seeder::add_tree_shadow(mygen::Tree &t)
{
    recalculate_planar_shadows(t.root, psm, 1);
    dsm.create(hm, psm);
}
int Seeder::joints_count(mygen::Branch *b)
{
    if (!b)
        return 0;
    int cnt = b->joints.size();
    for (mygen::Joint &j : b->joints)
    {
        for (mygen::Branch *br : j.childBranches)
            cnt += joints_count(br);
    }
    return cnt;
}
void Seeder::recalculate_planar_shadows(mygen::Branch *b, PlanarShadowsMap &psm, int level)
{
    if (!b || b->level > level)
        return;
    if (b->level < level)
    {
        for (mygen::Joint &j : b->joints)
        {
            for (mygen::Branch *br : j.childBranches)
                recalculate_planar_shadows(br,psm,level);
        }
    }
    else //b->level == level
    {
        if (b->joints.size() < 2)
            return;

        int cnt = joints_count(b);
        glm::vec3 pos = 0.5f*(b->joints.back().pos + b->joints.front().pos);
        float height = MAX(0,pos.y - heightmap->get_height(pos));
        float size = 0.5f*(length(b->joints.back().pos - b->joints.front().pos) + height);
        psm.set_occluder(pos,cnt,size,1);

    }
}