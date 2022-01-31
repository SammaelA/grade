#include "all_generators.h"

AbstractTreeGenerator *get_generator(std::string &generator_name)
{
    AbstractTreeGenerator *gen;
    if (generator_name == "proctree")
        gen = new Proctree::ProctreeGenerator();
    else if (generator_name == "simple_gen")
        gen = new SimpleTreeGenerator();
    else if (generator_name == "simpliest_gen")
        gen = new SimpliestTreeGenerator();
    else if (generator_name == "load_from_file")
        gen = new TreeLoaderBlk();
    else if (generator_name == "python_tree_gen")
        gen = new PythonTreeGen();
    else if (generator_name == "ge_gen")
        gen = new GETreeGenerator();
    else if (generator_name == "ge_gen_simplified")
        gen = new GETreeGeneratorSimplified();
    else
        gen = new mygen::TreeGenerator();

    return gen;
}

ParametersSet *get_default_parameters(std::string &generator_name)
{
    ParametersSet *gen;
    if (generator_name == "proctree")
        gen = new Proctree::Properties();
    else if (generator_name == "simple_gen")
        gen = new SimpleTreeStructureParameters();
    else if (generator_name == "simpliest_gen")
        gen = new SimpliestTreeStructureParameters();
    else if (generator_name == "load_from_file")
        gen = nullptr;//not implemented
    else if (generator_name == "python_tree_gen")
        gen = new WeberPennParameters();
    else if (generator_name == "ge_gen")
        gen = new GETreeParameters();
    else
        gen = new TreeStructureParameters();

    return gen; 
}