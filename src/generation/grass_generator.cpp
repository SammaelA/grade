#include "grass_generator.h"
#include "scene_generator.h"
#include <algorithm>
#include "graphics_utils/texture_atlas.h"
#include "tinyEngine/postfx.h"
#include "graphics_utils/texture_manager.h"
#include "graphics_utils/volumetric_occlusion.h"

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
        int patches_count = round(grass_quantity[i]*(full_size.x*full_size.y)/(PI*SQR(type.patch_size)));

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
    debugl(1, "added %d grass patches\n",grass_patches.size());
}

struct GridPoint
{
    static constexpr int point_size = 4;
    int patch_ids[point_size] = {-1,-1,-1,-1};
    std::list<Sphere2D>::iterator instances[point_size];
};
struct GrassGrid
{
public:
    static constexpr float precision = 1;

    GrassGrid(Cell &cell, float plant_av_sz)
    {
        cell_sz = (cell.bbox.max_pos - cell.bbox.min_pos);
        cell_pos = cell.bbox.min_pos;
        sz = glm::ivec2(precision/plant_av_sz*cell_sz) + glm::ivec2(1,1);
        p_sz = cell_sz/glm::vec2(sz);
        data = new GridPoint[sz.x*sz.y];
        //logerr("cell %d [%f %f] - [%f %f]",cell.id, cell.bbox.min_pos.x, cell.bbox.min_pos.y, cell.bbox.max_pos.x, cell.bbox.max_pos.y);
        //logerr("created grid %d %d [%f %f] - [%f %f]",sz.x, sz.y, center(0,0).x, center(0,0).y, center(sz.x-1,sz.x-1).x,
        //center(sz.x-1,sz.x-1).y);
    }
    ~GrassGrid()
    {
        if (data)
            delete[] data;
    }

    glm::ivec4 slice(AABB2D &box)
    {
        glm::vec2 grdf = (box.min_pos - cell_pos)/p_sz;
        int x0 = CLAMP(floor(grdf.x),0,sz.x);
        int y0 = CLAMP(floor(grdf.y),0,sz.y);

        grdf = (box.max_pos - cell_pos)/p_sz;
        int x1 = CLAMP(ceil(grdf.x),0,sz.x);
        int y1 = CLAMP(ceil(grdf.y),0,sz.y);
        //return glm::ivec3(0,0,sz.x-1,sz.y-1);
        return glm::ivec4(x0, y0, x1, y1);
    }
    GridPoint &get_point(int x, int y)
    {
        return data[y*sz.x + x];
    }
    glm::vec2 center(int x, int y)
    {
        return cell_pos + p_sz*glm::vec2(x + 0.5, y + 0.5);
    }

    GridPoint *data = nullptr;
    glm::ivec2 sz;
    glm::vec2 cell_sz;
    glm::vec2 cell_pos;
    glm::vec2 p_sz;
};

void GrassGenerator::generate_grass_in_cell(Cell &cell, Field_2d *occlusion)
{
    bool test_occlusion = (occlusion != nullptr);
    //test_occlusion = false;//TODO: test occlusion
    Normal normal = Normal(0,1);

    GrassGrid grid = GrassGrid(cell, used_grass_types[0].plant_size);

    for (int patch_id : cell.grass_patches)
    {
        GrassPatch &patch = grass_patches[patch_id];
        GrassType &type = used_grass_types[patch.grass_type];
        AABB2D pb = patch.sphere.intersect_bbox(cell.bbox);
        if (pb.empty())
            continue;
        float pb_sq = (pb.max_pos.x - pb.min_pos.x)*(pb.max_pos.y - pb.min_pos.y);
        int instances_cnt = round(patch.instances_max_cnt*(pb_sq/(PI*SQR(patch.sphere.r))));
        if (instances_cnt <= 0)
            continue;
        int new_instances = 0;
        float ipp_f_base = CLAMP(type.patch_density*(grid.p_sz.x*grid.p_sz.y), 1, GridPoint::point_size);
        float grid_diag = sqrt(SQR(grid.p_sz.x) + SQR(grid.p_sz.y));
        glm::ivec4 grid_slice = grid.slice(pb);
        //logerr("subgrid (%d %d) - (%d %d)",grid_slice.x, grid_slice.y, grid_slice.z, grid_slice.w);
        for (int x = grid_slice.x; x < grid_slice.z; x++)
        {
            for (int y = grid_slice.y; y < grid_slice.w; y++)
            {
                GridPoint &p = grid.get_point(x, y);
                glm::vec2 center = grid.center(x, y);
                if (patch.sphere.contains(center))
                {
                    float LS = type.light_sensivity;
                    float grow_chance = LS >= 0 ? 1 : 0;
                    if (test_occlusion && LS != 0)
                    {
                        float light = 1/(1 + occlusion->get_bilinear(glm::vec3(center.x, 0, center.y)));
                        if (light > 0.005)
                            grow_chance = LS > 0 ? pow(light, LS) : 1 - pow(light, -LS);
                        else
                            grow_chance = 0;
                    }
                    if (grow_chance < 0.01)
                        continue;
                    float ipp_f = ipp_f_base*grow_chance;
                    int ipp_i = floor(ipp_f);
                    float ipp_r = ipp_f - ipp_i;
                    int ipp = CLAMP(ipp_i + (urand() < ipp_r),0,GridPoint::point_size);
                    //logerr("cell %d patch %d (%d %d )trying to add %d instances %f %f %f",cell.id, patch_id, x, y, ipp,
                    //        grow_chance, center.x, center.y);

                    for (int point = 0; point < ipp; point++)
                    {
                        glm::vec2 pos = center + grid.p_sz*glm::vec2(urand(-0.5, 0.5), urand(-0.5, 0.5));
                        float r = MAX(type.plant_size + type.plant_size_std_dev*normal.get(),0.1f*type.plant_size);
                        bool placed = false;
                        int n = 0;
                        while (!placed && n < GridPoint::point_size)
                        {
                            if (p.patch_ids[n] < 0)
                            {
                                //add new instance
                                patch.grass_instances.push_back(Sphere2D(pos,r));
                                p.patch_ids[n] = patch_id;
                                p.instances[n] = --patch.grass_instances.end();
                                placed = true;
                                new_instances++;
                            }
                            else
                            {
                                GrassPatch &sec_p = grass_patches[p.patch_ids[n]];
                                GrassType &sec_type = used_grass_types[sec_p.grass_type];
                                float dist_q = (1 - length(pos - p.instances[n]->pos))/(grid_diag);// from 0 to 1
                                float push_a = type.push/(type.push + sec_type.push + 1e-6);
                                float push_b = 1 - push_a;

                                float rnd = urand();
                                if (rnd < dist_q*push_a)
                                {
                                    //replace existed instance with new
                                    sec_p.grass_instances.erase(p.instances[n]);
                                    patch.grass_instances.push_back(Sphere2D(pos,r));
                                    p.patch_ids[n] = patch_id;
                                    p.instances[n] = --patch.grass_instances.end();
                                    placed = true;
                                    new_instances++;
                                }
                                else if (rnd > 1 - dist_q*push_b || n == GridPoint::point_size - 1)
                                {
                                    //do not grow this instance
                                    placed = true;
                                }
                            }
                            n++;
                        }
                    }
                }
            }
        }
        /*
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
        */
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
    textureManager.save_png(grass_packed.grass_textures.tex(0),"grass_atlas");


    PostFx copy_alpha = PostFx("alpha_split_alpha.fs");
    TextureAtlas atl_alpha = TextureAtlas(grass_tex_size.x,grass_tex_size.y*used_grass_types.size(),1);
    atl_alpha.set_grid(grass_tex_size.x, grass_tex_size.y, false);
    for (auto &t : used_grass_types)
    {
        int tex_id = atl_alpha.add_tex();
        atl_alpha.target_slice(tex_id, 0);
        copy_alpha.use();
        copy_alpha.get_shader().texture("tex",t.texture);
        copy_alpha.render();
    }
    atl_alpha.gen_mipmaps();
    textureManager.save_png(atl_alpha.tex(0),"grass_atlas_alpha");


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
            //logerr("added grass instance %d (%f %f %f) r= %f",in_cnt,pos3.x,pos3.y,pos3.z,in.r);
            in_cnt++;
        }
    }
}