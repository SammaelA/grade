#include "python_tree_gen.h"
#include "../common_utils/python_interaction.h"
#include "generation/grove_generation_utils.h"
#include "weber_penn_parameters.h"
#include <iosfwd>
#include <fstream>
#include <iostream>
#include <sstream>

void PythonTreeGen::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{

    std::string scripts_path = "./scripts/tree_gen"; 
    std::string config_file_name = scripts_path + "/params.py";
    std::string params_string;
    WeberPennParameters dummy;
    WeberPennParameters *params = dynamic_cast<WeberPennParameters *>(ggd.types[0].params);
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
                logerr("unable to load file %s",settings_path.c_str());
                throw std::exception();
            }
            iss << f.rdbuf();
            params_string = iss.str();
            f.close();
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            logerr("cannot open file %s with settings for tree %s. Default parameters are used.",settings_path.c_str(),params->name.c_str());
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
    BlkManager man;

    for (int i=0;i<ggd.trees_count;i++)
    {
        std::string blk_path = trees_directory + "/" + "tree.bsg";
        
        
        ph.run_script("sa_test","sa_test "+blk_path);
        ph.finish_script();

        Block b;

        //Seeder seeder = Seeder(ggd,10,&h);
        //std::vector<Seed> seeds;
        //seeder.choose_places_for_seeds(1,seeds);
        //if (seeds.size() >= 1)
        {
            glm::vec3 pos = glm::vec3(urand(-ggd.size.x,ggd.size.x),0,urand(-ggd.size.z,ggd.size.z));
            pos.y = h.get_height(pos);
            float scale = 10;
            b.add_vec3("position",pos);
            b.add_double("scale",scale);

            man.load_block_from_file(blk_path, b);
            load_tree(ggd,trees_external[i],h,b);
            trees_taken++;
        }
    }
}