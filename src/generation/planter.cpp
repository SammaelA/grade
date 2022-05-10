#include "planter.h"

Planter::Planter(LightVoxelsCube *_voxels, Heightmap *_heightmap, GroveMask *_mask, GroveMask *_bm,
                 glm::vec3 _center, glm::vec2 _size,
                 float base_density, int max_saplings, float _cell_size):
density(_center, _size, _cell_size),
occlusion(_center, _size, _cell_size)
{ 
    voxels = _voxels;
    heightmap = _heightmap;
    center = _center;
    size = _size;
    cell_size = _cell_size;
    mask = _mask;
    biome_mask = _bm;
    AABB box = voxels->get_bbox();
    glm::vec4 bord = heightmap->get_borders();//(min_pos,max_pos)
    if (center.x - size.x < box.min_pos.x ||
        center.z - size.y < box.min_pos.z ||
        center.x + size.x > box.max_pos.x ||
        center.z + size.y > box.max_pos.z)
    {
        logerr("planter error: plant area is outside the provided voxels array");
        return;
    }
    if (center.x - size.x < bord.x ||
        center.z - size.y < bord.y ||
        center.x + size.x > bord.z ||
        center.z + size.y > bord.w)
    {
        debug("planter warning: plant area is outside the provided heightmap\n");
    }

    density.fill_const(base_density);
    //подача маски, проекция маски на density
    occlusion.fill_const(0);
    saplings_left = MIN(max_saplings, base_density*size.x*size.y);
}

std::vector<glm::vec3> Planter::get_saplings()
{
    int cnt = saplings_left >= 2 ? saplings_left/2 : saplings_left;
    std::vector<glm::vec3> saplings;
    //спроецировать воксельный массив на occlusion
    std::function<float(glm::vec2 &)> func = [&](glm::vec2 &p) ->float
    {
        return voxels->get_occlusion_projection(glm::vec3(p.x,0,p.y));
    };
    occlusion.fill_func(func);

    for (int i=0;i<cnt;i++)
    {
        const int max_points = 50;
        const int max_tries = 50;
        int tries = 0;
        float sum = 0;
        std::vector<glm::vec3> points;
        while (points.size() < max_points && tries < max_tries)
        {
            float x = center.x + urand(-1,1)*size.x;
            float y = center.z + urand(-1,1)*size.y;
            glm::vec3 ps = glm::vec3(x,0,y);
            float occ = (1 /(1 + occlusion.get_bilinear(ps)))*density.get_bilinear(ps)*mask->get_bilinear(ps);
            if (biome_mask)
                occ *= biome_mask->get_bilinear(ps);
            if (occ > 1e-4)
            {
                points.push_back(glm::vec3(x,y,occ));
                sum += occ;
            }
            //logerr("occ is got %f %f %f",occlusion.get_bilinear(ps), density.get_bilinear(ps), mask->get_bilinear(ps));
            tries++;
        }
        if (points.empty())
        {
            debug("Warning: plant failed to find possible places for seed. Probably planting area is fully occluded\n");
        }
        else
        {
            float rnd = urand(0,sum);
            for (auto &p : points)
            {
                //logerr("p %f %f rnd %f %f",p.x,p.y,rnd,p.z);
                if (rnd < p.z)
                {
                    glm::vec3 pos = glm::vec3(p.x,0,p.y);
                    pos.y = heightmap->get_bilinear(pos);
                    saplings.push_back(pos);
                    saplings_planted.push_back(pos);
                    occlusion.set(pos,1e9);
                    debugl(1, "placed sapling %f %f %f\n",pos.x,pos.y,pos.z);
                    break;
                }
                else
                {
                    rnd -= p.z;
                }
            }
        }
    }
    saplings_left -= cnt;
    return saplings;
};