#include "python_tree_gen.h"
#include "../common_utils/python_interaction.h"
#include "generation/grove_generation_utils.h"
#include "weber_penn_parameters.h"
#include <iosfwd>
#include <fstream>
#include <iostream>
#include <sstream>

void PythonTreeGen::plant_tree(glm::vec3 pos, TreeTypeData *type)
{
  tree_saplings.push_back(std::pair<glm::vec3, TreeTypeData *>(pos, type));
}

void PythonTreeGen::finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels)
{
  if (tree_saplings.empty())
    return;

  std::string scripts_path = "./scripts/tree_gen";
  std::string config_file_name = scripts_path + "/params.py";
  std::string params_string;
  WeberPennParameters dummy;
  WeberPennParameters *params = dynamic_cast<WeberPennParameters *>(tree_saplings[0].second->get_params());
  if (!params)
    params = &dummy;
  if (params->settings_already_in_file)
  {
    std::string settings_path = scripts_path + "/" + params->name + ".py";
    try
    {
      std::fstream f(settings_path);
      std::stringstream iss;
      if (f.fail())
      {
        logerr("unable to load file %s", settings_path.c_str());
        throw std::exception();
      }
      iss << f.rdbuf();
      params_string = iss.str();
      f.close();
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
      logerr("cannot open file %s with settings for tree %s. Default parameters are used.", settings_path.c_str(), params->name.c_str());
    }
  }
  if (params_string.empty())
  {
    params_string = params->convert_to_python_list();
  }
  std::ofstream out(config_file_name);
  out << params_string;
  out.close();

  PythonHelper ph;
  ph.init(scripts_path);

  for (int i = 0; i < tree_saplings.size(); i++)
  {
    std::string blk_path = trees_directory + "/" + "tree.bsg";

    ph.run_script("sa_test", "sa_test " + blk_path);
    ph.finish_script();

    Block b;

    glm::vec3 pos = tree_saplings[i].first;
    float scale = 10;
    b.add_vec3("position", pos);
    b.add_double("scale", scale);

    load_block_from_file(blk_path, b);
    load_tree(tree_saplings[i].second, trees_external[i], b);
    trees_taken++;
  }
}