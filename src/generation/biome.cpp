#include "biome.h"
#include "save_utils/blk.h"
#include "common_utils/utility.h"
#include "metainfo_manager.h"
#include "common_utils/bbox.h"
#include "graphics_utils/texture_manager.h"
#include "grove_generation_utils.h"
#include "common_utils/distribution.h"

int gr = 0;
void load_biome_types(Block *bl, std::vector<std::pair<int, float>> &types)
{
    if (!bl)
    {
       logerr("error in reading biome description. Cannot find types block");
       return; 
    }

    for (int i=0;i<bl->size();i++)
    {
        float dens = bl->get_double(i,0);
        if (dens > 0)
        {
            std::string name = bl->get_name(i);
            int id = gr ? metainfoManager.get_grass_type_id_by_name(name) : metainfoManager.get_tree_type_id_by_name(name);
            if (id >= 0)
                types.push_back(std::pair<int, float>(id,dens));
        }
    }

    if (types.empty())
    {
       logerr("error in reading biome description. Types block have 0 valid types");   
    }
}

void Biome::load_from_blk(Block &b)
{
    for (int i=0;i<b.size();i++)
    {
        Block *bl = b.get_block(i);
        if (!bl)
        {
            logerr("error in reading biome description. Expected block but got %d",(int)b.get_type(i));
            continue;
        }
        VegClass &vc = (b.get_name(i) == "trees") ? trees : ((b.get_name(i) == "plants") ? plants : grass);
        if (b.get_name(i) == "grass")
            gr = 1;
        else
            gr = 0;
        std::string cov_type = bl->get_string("coverage","");
        if (cov_type == "full")
        {
            if (vc.main_density < 0)
            {
               vc.main_density = bl->get_double("density", -1);
               load_biome_types(bl->get_block("types"), vc.main_types); 
            }
            else
            {
                logerr("warning in reading biome description. Expected only 1 \"full\" coverage block");
            }
        }
        else if (cov_type == "patches")
        {
            vc.patch_descs.emplace_back();
            PatchDesc &p = vc.patch_descs.back();
            
            p.size = bl->get_double("size", p.size);
            p.size_std_dev = bl->get_double("size_std_dev", p.size_std_dev);
            p.density = bl->get_double("density", p.density);
            p.density_std_dev = bl->get_double("density_std_dev", p.density_std_dev);
            p.coverage_part = bl->get_double("coverage_part", p.coverage_part);
            p.push = bl->get_double("push", p.push);

            load_biome_types(bl->get_block("types"), p.types);
        }
        else
        {
            logerr("error in reading biome description. Expected coverage_type to be \"full\" or \"patches\". Got %s",cov_type.c_str());
        }
    }
}

BiomeMap::BiomeMap()
{

}
BiomeMap::~BiomeMap()
{
    if (data)
        delete[] data;
}
void BiomeMap::create(AABB2D _bbox, float _pixel_size)
{
    if (data)
        delete[] data;
    pixel_size = _pixel_size;
    bbox = _bbox;
    w = ceil((bbox.max_pos.x - bbox.min_pos.x)/pixel_size);
    h = ceil((bbox.max_pos.y - bbox.min_pos.y)/pixel_size);
    data = new biome_type_t[w*h];
}

    int BiomeMap::get(int x, int y)
    {
        return data[y*w + x];
    }
    int BiomeMap::get(glm::vec2 pos)
    {
        int x = (pos.x - bbox.min_pos.x)/pixel_size;
        int y = (pos.y - bbox.min_pos.y)/pixel_size;
        if (x < 0 || x >= w || y < 0 || y >= h)
            return -1;
        else
            return get(x, y);
    }
    int BiomeMap::get(glm::vec3 pos)
    {
        return get(glm::vec2(pos.x, pos.z));
    }
    
    void BiomeMap::set_rect(AABB2D box, int id)
    {
        int x0 = CLAMP(ceil(box.min_pos.x - bbox.min_pos.x)/pixel_size,0,w);
        int y0 = CLAMP(ceil(box.min_pos.y - bbox.min_pos.y)/pixel_size,0,h);
        int x1 = CLAMP(ceil(box.max_pos.x - bbox.min_pos.x)/pixel_size,0,w);
        int y1 = CLAMP(ceil(box.max_pos.y - bbox.min_pos.y)/pixel_size,0,h);

        for (int i=y0;i<y1;i++)
        {
            for (int j=x0;j<x1;j++)
            {
                data[i*w + j] = id;
            }
        }
    }
    void BiomeMap::set_round(glm::vec2 pos, float inner_r, float outer_r, int id)
    {
        AABB2D box = AABB2D(pos - outer_r*glm::vec2(1,1), pos + outer_r*glm::vec2(1,1));
        
        float ir2p = SQR(inner_r/pixel_size);
        float or2p = SQR(outer_r/pixel_size);

        int x0 = CLAMP(ceil(box.min_pos.x - bbox.min_pos.x)/pixel_size,0,w);
        int y0 = CLAMP(ceil(box.min_pos.y - bbox.min_pos.y)/pixel_size,0,h);
        int x1 = CLAMP(ceil(box.max_pos.x - bbox.min_pos.x)/pixel_size,0,w);
        int y1 = CLAMP(ceil(box.max_pos.y - bbox.min_pos.y)/pixel_size,0,h);

        float xc = 0.5*x0 + 0.5*x1;
        float yc = 0.5*y0 + 0.5*y1;

        if (or2p/ir2p > 1.01)
        {
            //smooth edges
            for (int i=y0;i<y1;i++)
            {
                for (int j=x0;j<x1;j++)
                {
                    float d2 = SQR(i - yc) + SQR(j - xc);
                    if (d2 > or2p)
                    {
                        continue;
                    }
                    else if (d2 < ir2p)
                    {
                        data[i*w + j] = id;
                    }
                    else if ((d2 - ir2p)/(or2p - ir2p) < urand())
                    {
                        data[i*w + j] = id;
                    }
                }
            }
        }
        else
        {
            //sharp edges
            for (int i=y0;i<=y1;i++)
            {
                for (int j=x0;j<=x1;j++)
                {
                    float d2 = SQR(i - xc) + SQR(j - yc);
                    if (d2 < ir2p)
                    {
                        data[i*w + j] = id;
                    }
                }
            }
        }
    }

    void BiomeMap::save_as_image(std::string name)
    {
        if (!data)
            return;
        unsigned char *image_data = new unsigned char[3*w*h];
        for (int i=0;i<h;i++)
        {
            for (int j=0;j<w;j++)
            {
                int type = ((int)data[i*w + j]) % 8 + 1;
                image_data[3*(i*w + j)] = 255*(type % 2);
                image_data[3*(i*w + j)+1] = 255*(type / 2 % 2);
                image_data[3*(i*w + j)+2] = 255*(type / 4 % 2);
                //image_data[4*(i*w + j)+3] = 1;
            }
        }
        textureManager.save_bmp_raw(image_data, w, h, 3, name);

        delete[] image_data;
    }

    void BiomeMap::get_stat(std::vector<std::pair<int,int>> &stat, AABB2D box)
    {
        int x0 = CLAMP(ceil(box.min_pos.x - bbox.min_pos.x)/pixel_size,0,w);
        int y0 = CLAMP(ceil(box.min_pos.y - bbox.min_pos.y)/pixel_size,0,h);
        int x1 = CLAMP(ceil(box.max_pos.x - bbox.min_pos.x)/pixel_size,0,w);
        int y1 = CLAMP(ceil(box.max_pos.y - bbox.min_pos.y)/pixel_size,0,h);

        stat = {};
        for (int i=y0;i<y1;i++)
        {
            for (int j=x0;j<x1;j++)
            {
                int b = data[i*w + j];
                for (auto &p : stat)
                {
                    if (p.first == b)
                    {
                        p.second++;
                        b = -1;
                        break;
                    }
                }
                if (b >= 0)
                    stat.push_back(std::pair<int,int>(b,1));
            }
        }
    }

    void BiomeMap::set_mask(GroveMask &mask, int biome_id)
    {
        std::function<float(glm::vec2 &)> func = [&](glm::vec2 &p) -> float
        {
            return (get(p) == biome_id) ? 1.0 : 0.0;
        };
        mask.fill_func(func);
    }