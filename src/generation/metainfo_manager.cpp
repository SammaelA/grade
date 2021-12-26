#include "metainfo_manager.h"
#include "save_utils/blk.h"
#include "tree_generators/GE_generator.h"
#include "tree_generators/mygen_parameters.h"
#include "tree_generators/simple_generator.h"
#include "tree_generators/simpliest_generator.h"

MetainfoManager metainfoManager;
TreeTypeData dummy1;
GrassType dummy2;
Biome dummy3;

void MetainfoManager::load_tree_types()
{
    int id = 0;
    BlkManager man;
    Block ge_gen_types, my_gen_types, simple_gen_types, simpliest_gen_types;
    man.load_block_from_file("ge_gen_presets.blk", ge_gen_types);
    man.load_block_from_file("my_gen_presets.blk", my_gen_types);
    man.load_block_from_file("simple_gen_presets.blk", simple_gen_types);
    man.load_block_from_file("simpliest_gen_presets.blk", simpliest_gen_types);

    for (int i = 0; i < ge_gen_types.size(); i++)
    {
        Block *bl = ge_gen_types.get_block(i);
        if (bl)
        {
            std::string name = ge_gen_types.get_name(i);
            std::string wood_tex_name = bl->get_string("wood_tex_name", "wood");
            std::string leaf_tex_name = bl->get_string("leaf_tex_name", "leaf");
            GETreeParameters *params = new GETreeParameters();
            params->load_from_blk(*bl);
            TreeTypeData type = TreeTypeData(id, params, wood_tex_name, leaf_tex_name);
            type.generator_name = "ge_gen";

            type.type_id = tree_types.size();
            tree_type_id_by_name.emplace(name, tree_types.size());
            tree_types.push_back(type);
            id++;
        }
    }

    for (int i = 0; i < my_gen_types.size(); i++)
    {
        Block *bl = my_gen_types.get_block(i);
        if (bl)
        {
            std::string name = my_gen_types.get_name(i);
            std::string wood_tex_name = bl->get_string("wood_tex_name", "wood");
            std::string leaf_tex_name = bl->get_string("leaf_tex_name", "leaf");
            TreeStructureParameters *params = new TreeStructureParameters();
            params->load_from_blk(*bl);
            TreeTypeData type = TreeTypeData(id, params, wood_tex_name, leaf_tex_name);
            type.generator_name = "my_gen";

            type.type_id = tree_types.size();
            tree_type_id_by_name.emplace(name, tree_types.size());
            tree_types.push_back(type);
            id++;
        }
    }
    
    for (int i = 0; i < simple_gen_types.size(); i++)
    {
        Block *bl = simple_gen_types.get_block(i);
        if (bl)
        {
            std::string name = simple_gen_types.get_name(i);
            std::string wood_tex_name = bl->get_string("wood_tex_name", "wood");
            std::string leaf_tex_name = bl->get_string("leaf_tex_name", "leaf");
            SimpleTreeStructureParameters *params = new SimpleTreeStructureParameters();
            params->load_from_blk(*bl);
            TreeTypeData type = TreeTypeData(id, params, wood_tex_name, leaf_tex_name);
            type.generator_name = "simple";

            type.type_id = tree_types.size();
            tree_type_id_by_name.emplace(name, tree_types.size());
            tree_types.push_back(type);
            id++;
        }
    }
        
    for (int i = 0; i < simpliest_gen_types.size(); i++)
    {
        Block *bl = simpliest_gen_types.get_block(i);
        if (bl)
        {
            std::string name = simpliest_gen_types.get_name(i);
            std::string wood_tex_name = bl->get_string("wood_tex_name", "wood");
            std::string leaf_tex_name = bl->get_string("leaf_tex_name", "leaf");
            SimpliestTreeStructureParameters *params = new SimpliestTreeStructureParameters();
            params->load_from_blk(*bl);
            TreeTypeData type = TreeTypeData(id, params, wood_tex_name, leaf_tex_name);
            type.generator_name = "simpliest";
            logerr("loaded simpliest type %s", name.c_str());

            type.type_id = tree_types.size();
            tree_type_id_by_name.emplace(name, tree_types.size());
            tree_types.push_back(type);
            id++;
        }
    }
}

void MetainfoManager::load_grass_types()
{
    BlkManager man;
    Block grass_types_blk;

    man.load_block_from_file("grass_presets.blk", grass_types_blk);
    for (int i = 0; i < grass_types_blk.size(); i++)
    {
        Block *bl = grass_types_blk.get_block(i);
        std::string name = grass_types_blk.get_name(i);
        if (bl)
        {
            GrassType gt;
            gt.load_from_blk(*bl);

            gt.id = grass_types.size();
            grass_type_id_by_name.emplace(name, grass_types.size());
            grass_types.push_back(gt);
        }
    }
}

void MetainfoManager::load_biomes()
{
    BlkManager man;
    Block biomes_blk;

    man.load_block_from_file("biomes.blk", biomes_blk);
    for (int i = 0; i < biomes_blk.size(); i++)
    {
        Block *bl = biomes_blk.get_block(i);
        std::string name = biomes_blk.get_name(i);
        if (bl)
        {
            Biome gt;
            gt.load_from_blk(*bl);

            gt.id = biomes.size();
            biome_id_by_name.emplace(name, biomes.size());
            biomes.push_back(gt);
        }
    }
}

bool MetainfoManager::reload_all(std::string dir)
{
    loaded = false;
    tree_type_id_by_name = {};
    tree_types = {};
    grass_type_id_by_name = {};
    grass_types = {};
    biome_id_by_name = {};
    biomes = {};

    load_tree_types();
    load_grass_types();
    load_biomes();

    loaded = tree_types.size() > 0 && grass_types.size() > 0 && biomes.size() > 0;
    return loaded;
}

TreeTypeData &MetainfoManager::get_tree_type(std::string name)
{
    auto it = tree_type_id_by_name.find(name);
    if (it == tree_type_id_by_name.end())
    {
        logerr("invalid tree type name = %s", name.c_str());
        return dummy1;
    }
    else
        return get_tree_type(it->second);
}

TreeTypeData &MetainfoManager::get_tree_type(int id)
{
    if (id < 0 || id >= tree_types.size())
    {
        logerr("invalid tree type id = %d", id);
        return dummy1;
    }
    else
        return tree_types[id];
}

int MetainfoManager::get_tree_type_id_by_name(std::string name)
{
    auto it = tree_type_id_by_name.find(name);
    if (it == tree_type_id_by_name.end())
    {
        logerr("invalid tree type name = %s", name.c_str());
        return -1;
    }
    else
        return it->second;
}

GrassType &MetainfoManager::get_grass_type(std::string name)
{
    auto it = grass_type_id_by_name.find(name);
    if (it == grass_type_id_by_name.end())
    {
        logerr("invalid grass type name = %s", name.c_str());
        return dummy2;
    }
    else
        return get_grass_type(it->second);
}

GrassType &MetainfoManager::get_grass_type(int id)
{
    if (id < 0 || id >= grass_types.size())
    {
        logerr("invalid tree type id = %d", id);
        return dummy2;
    }
    else
        return grass_types[id];
}

int MetainfoManager::get_grass_type_id_by_name(std::string name)
{
    auto it = grass_type_id_by_name.find(name);
    if (it == grass_type_id_by_name.end())
    {
        logerr("invalid grass type name = %s", name.c_str());
        return -1;
    }
    else
        return it->second;
}

Biome &MetainfoManager::get_biome(std::string name)
{
    auto it = biome_id_by_name.find(name);
    if (it == biome_id_by_name.end())
    {
        logerr("invalid biome name = %s", name.c_str());
        return dummy3;
    }
    else
        return get_biome(it->second);
}

Biome &MetainfoManager::get_biome(int id)
{
    if (id < 0 || id >= biomes.size())
    {
        logerr("invalid tree type id = %d", id);
        return dummy3;
    }
    else
        return biomes[id];
}

int MetainfoManager::get_biome_id_by_name(std::string name)
{
    auto it = biome_id_by_name.find(name);
    if (it == biome_id_by_name.end())
    {
        logerr("invalid grass type name = %s", name.c_str());
        return -1;
    }
    else
        return it->second;
}
