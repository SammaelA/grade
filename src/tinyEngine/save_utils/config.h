#pragma once
#include "../../parameter.h"
#include <map>
class Config
{
public:
    bool load_config();
    TreeStructureParameters get(std::string name);
private:
    std::map<std::string,TreeStructureParameters> presets;
};