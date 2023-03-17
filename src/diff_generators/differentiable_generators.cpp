#include "differentiable_generators.h"
#include "diff_geometry_generation.h"
#include "dishes_generator.h"
#include "building_generator.h"
#include <cppad/cppad.hpp>
#include <map>

namespace dgen
{
  GeneratorDescription get_generator_by_name(std::string gen_name)
  {
    static bool loaded = false;
    static std::map<std::string, GeneratorDescription> generators;

    if (!loaded)
    {
      loaded = true;
      {
        GeneratorDescription gen;
        std::string name = "dishes";
        gen.name = name;
        gen.generator = create_cup;
        gen.params_regularizer = parameters_cup_reg;
        gen.model_regularizer = default_model_reg;
        gen.generator_description_blk_path = "dishes_gen_parameters_description.blk";
        gen.presets_blk_path = "dishes_gen_presets.blk";
        generators.emplace(name, gen);
        debug("generator description %s loaded\n", name.c_str());
      }
      {
        GeneratorDescription gen;
        std::string name = "buildings";
        gen.name = name;
        gen.generator = create_building;
        gen.params_regularizer = default_parameters_reg;
        gen.model_regularizer = default_model_reg;
        gen.generator_description_blk_path = "buildings_gen_parameters_description.blk";
        gen.presets_blk_path = "buildings_gen_presets.blk";
        generators.emplace(name, gen);
        debug("generator description %s loaded\n", name.c_str());
      }
    }

    auto it = generators.find(gen_name);
    if (it == generators.end())
    {
      logerr("Cannot find generator %s", gen_name.c_str());
      logerr("Possible variants are: ");
      for (auto p : generators)
        logerr("%s",p.first.c_str());
      return GeneratorDescription();
    }
    else
      return it->second;
  }
};