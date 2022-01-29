#pragma once
#include "core/grass.h"
#include "core/tree.h"
#include "biome.h"

class MetainfoManager
{
public:
    bool reload_all(std::string dir = ".");
    bool save_all(std::string dir = ".");
    
    void add_tree_type(TreeTypeData &ttd, std::string name);
    TreeTypeData &get_tree_type(std::string name);
    TreeTypeData &get_tree_type(int id);
    int get_tree_type_id_by_name(std::string name);
    std::vector<TreeTypeData> get_all_tree_types() {return tree_types;}

    GrassType &get_grass_type(std::string name);
    GrassType &get_grass_type(int id);
    int get_grass_type_id_by_name(std::string name);
    std::vector<GrassType> get_all_grass_types() {return grass_types;}

    Biome &get_biome(std::string name);
    Biome &get_biome(int id);
    int get_biome_id_by_name(std::string name);
    std::vector<Biome> get_all_biomes() {return biomes;}
private:
    void load_tree_types();
    void load_grass_types();
    void load_biomes();
    
    void save_tree_types();

    std::map<std::string, int> tree_type_id_by_name;
    std::vector<TreeTypeData> tree_types;
    
    std::map<std::string, int> grass_type_id_by_name;
    std::vector<GrassType> grass_types;

    std::map<std::string, int> biome_id_by_name;
    std::vector<Biome> biomes;

    bool loaded = false;
};

extern MetainfoManager metainfoManager;