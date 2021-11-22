#include "grass_generator.h"
#include "scene_generator.h"
#include <algorithm>
#include "graphics_utils/texture_atlas.h"
#include "tinyEngine/postfx.h"
#include "graphics_utils/texture_manager.h"

void GrassGenerator::set_grass_types(const std::map<std::string, GrassType> &grass_types, Block &grass_settings)
{
    for (int i=0;i<grass_settings.size();i++)
    {
        std::string name = grass_settings.get_name(i);
        float q = grass_settings.get_double(i, -1);
        if (q > 0)
        {
            auto it = grass_types.find(name);
            if (it == grass_types.end())
            {
                logerr("cannot find grass type \"%s\" in types list",name.c_str());
            }
            else
            {
                used_grass_types.push_back(it->second);
                grass_quantity.push_back(q);
            }
        }
    }
}
void GrassGenerator::add_patch(std::vector<Cell> &cells, int cells_x, int cells_y, int cell_id, int type_id)
{
    grass_patches.emplace_back();
    Cell &main_cell = cells[cell_id];
    GrassType &type = used_grass_types[type_id];

    //choose parameters
    glm::vec2 pos = glm::vec2(urand(main_cell.bbox.min_pos.x,main_cell.bbox.max_pos.x),
                              urand(main_cell.bbox.min_pos.y,main_cell.bbox.max_pos.y));
    Normal n = Normal(0,1);
    float size = type.patch_size + type.patch_size_std_dev*n.get();
    float density = type.patch_density + type.patch_density_std_dev*n.get();
    grass_patches.back().sphere = Sphere2D(pos, size);
    grass_patches.back().grass_type = type_id;
    grass_patches.back().instances_max_cnt = MAX(1,density*SQR(size)/SQR(type.plant_size));

    //add patch to cells that intersects it

    for (auto &c : cells)
    {
        //TODO: make it faster
        if (c.bbox.intersects(grass_patches.back().sphere))
            c.grass_patches.push_back(grass_patches.size()-1);
    }

}
void GrassGenerator::prepare_grass_patches(std::vector<Cell> &cells, int cells_x, int cells_y)
{
    if (cells_x <= 0 || cells_y <= 0 || used_grass_types.empty())
        return;
    glm::vec2 cell_size = cells[0].bbox.max_pos - cells[0].bbox.min_pos;
    glm::vec2 full_size = glm::vec2(cells_x, cells_y)*cell_size;
    for (int i = 0;i<used_grass_types.size();i++)
    {
        auto &type = used_grass_types[i];
        int patches_count = grass_quantity[i]*(full_size.x*full_size.y)/(PI*SQR(type.patch_size));
        logerr("%f %d %f %f %f",grass_quantity[i],patches_count, full_size.x, full_size.y, type.patch_size);
        if (patches_count <= 0)
            continue;
        
        int full_iters = patches_count/(cells_x*cells_y);
        int patches_left = patches_count - full_iters*(cells_x*cells_y);

        for (int iter =0;iter<full_iters;iter++)
        {
            for (int c=0;c<cells.size();c++)
            {
                add_patch(cells, cells_x, cells_y, c, i);
            }
        }
        std::vector<int> rand_cells = std::vector<int>(cells.size(),0);
        for (int c=0;c<cells.size();c++)
        {
            rand_cells[c] = c;
        }
        std::random_shuffle(rand_cells.begin(),rand_cells.end());
        for (int j=0;j<patches_left;j++)
        {
            add_patch(cells, cells_x, cells_y, rand_cells[j], i);
        }
    }
    logerr("added %d grass patches",grass_patches.size());
}
void GrassGenerator::generate_grass_in_cell(Cell &cell, LightVoxelsCube *occlusion)
{
    bool test_occlusion = (occlusion != nullptr);
    test_occlusion = false;//TODO: test occlusion
    //int dots_mult = 3;
    Normal normal = Normal(0,1);
    for (int patch_id : cell.grass_patches)
    {
        GrassPatch &patch = grass_patches[patch_id];
        GrassType &type = used_grass_types[patch.grass_type];
        AABB2D pb = patch.sphere.intersect_bbox(cell.bbox);
        if (pb.empty())
            continue;
        float pb_sq = (pb.max_pos.x - pb.min_pos.x)*(pb.max_pos.y - pb.min_pos.y);
        int instances_cnt = round(patch.instances_max_cnt*(pb_sq/(PI*SQR(patch.sphere.r))));
        for (int i=0;i<instances_cnt;i++)
        {
            glm::vec2 pos = glm::vec2(urand(pb.max_pos.x, pb.min_pos.x),urand(pb.max_pos.y, pb.min_pos.y));
            float r = type.plant_size + type.plant_size_std_dev*normal.get();

            bool alive = patch.sphere.contains(pos);
            //TODO: test for occluson and intersection with other instances
            
            if (alive)
            {
                patch.grass_instances.push_back(Sphere2D(pos,r));
            }
        }
    }
}
void GrassGenerator::pack_all_grass(GrassPacked &grass_packed, Heightmap &h)
{  
    PostFx copy = PostFx("copy.fs");
    grass_packed.grass_instances = {};
    glm::ivec2 grass_tex_size = glm::ivec2(512, 512);
    
    {
        TextureAtlas atl = TextureAtlas(grass_tex_size.x,grass_tex_size.y*used_grass_types.size(),1);
        grass_packed.grass_textures = atl;
    }
    grass_packed.grass_textures.set_grid(grass_tex_size.x, grass_tex_size.y, false);
    for (auto &t : used_grass_types)
    {
        int tex_id = grass_packed.grass_textures.add_tex();
        grass_packed.grass_textures.target_slice(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex",t.texture);
        copy.render();

        grass_packed.grass_instances.push_back(std::pair<int, std::vector<GrassInstanceData>>(tex_id,{}));
    }
    grass_packed.grass_textures.gen_mipmaps();
    textureManager.save_bmp(grass_packed.grass_textures.tex(0),"grass_atlas");

    int in_cnt = 0;
    for (auto &patch : grass_patches)
    {
        std::vector<GrassInstanceData> &instances = grass_packed.grass_instances[patch.grass_type].second;
        instances.reserve(patch.grass_instances.size());
        for (auto &in : patch.grass_instances)
        {
            glm::vec3 pos3 = glm::vec3(in.pos.x, 0, in.pos.y);
            pos3.y = h.get_height(pos3);
            instances.push_back(GrassInstanceData(pos3,in.r,urand(0, 2*PI)));
            logerr("added grass instance %d (%f %f %f) r= %f",in_cnt,pos3.x,pos3.y,pos3.z,in.r);
            in_cnt++;
        }
    }
}