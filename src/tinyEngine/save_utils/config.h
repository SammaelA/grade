#pragma once
#include "../../parameter.h"
#include "../../generators/mygen_parameters.h"
#include <map>

struct GroveGenerationData;
struct TreeTypeData;
class Config
{
public:
    Config();
    ~Config();
    bool load_config();
    bool load_ggds();
    TreeStructureParameters get(std::string name);
    GroveGenerationData get_ggd(std::string name);
private:
    std::map<std::string,TreeStructureParameters> presets;
    std::map<std::string,GroveGenerationData> ggds;
    std::vector<TreeStructureParameters *> params;
};