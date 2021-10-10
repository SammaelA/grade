#include "python_tree_gen.h"
#include "tinyEngine/utility/python_interaction.h"
#include "grove_generation_utils.h"

void PythonTreeGen::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    PythonHelper ph;
    ph.init("./scripts/tree_gen");
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