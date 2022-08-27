#include "grass.h"
#include "common_utils/blk.h"
#include "tinyEngine/engine.h"
#include <string>

void GrassType::load_from_blk(Block &block)
{
    std::string texture_name = block.get_string("texture","grass");
    texture = engine::textureManager->get(texture_name);

    patch_size = block.get_double("patch_size",patch_size);
    patch_size_std_dev = block.get_double("patch_size_std_dev",patch_size_std_dev);
    patch_density = block.get_double("patch_density",patch_density);
    patch_density_std_dev = block.get_double("patch_density_std_dev",patch_density_std_dev);
    plant_size = block.get_double("plant_size",plant_size);
    plant_size_std_dev = block.get_double("plant_size_std_dev",plant_size_std_dev);
    light_sensivity = block.get_double("light_sensivity",light_sensivity);
    push = block.get_double("push",push);
    model_name = block.get_string("model_name",model_name);
}

GrassType::GrassType(): texture(engine::textureManager->empty())
{
    
}