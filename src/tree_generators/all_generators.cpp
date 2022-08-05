#include "all_generators.h"

AbstractTreeGenerator *get_generator(std::string &generator_name)
{
    AbstractTreeGenerator *gen;
    if (generator_name == "simple_gen")
        gen = new SimpleTreeGenerator();
    else if (generator_name == "simpliest_gen")
        gen = new SimpliestTreeGenerator();
    else if (generator_name == "load_from_file")
        gen = new TreeLoaderBlk();
    else if (generator_name == "ge_gen")
        gen = new GETreeGenerator();
    else if (generator_name == "ge_gen_simplified")
        gen = new GETreeGeneratorSimplified();
    else if (generator_name == "weber_penn_gen")
        gen = new WeberPennGenerator();
    else
        gen = nullptr;

    return gen;
}

ParameterSet *get_default_parameters(std::string &generator_name)
{
    ParameterSet *gen;
    if (generator_name == "simple_gen")
        gen = new SimpleTreeStructureParameters();
    else if (generator_name == "simpliest_gen")
        gen = new SimpliestTreeStructureParameters();
    else if (generator_name == "load_from_file")
        gen = nullptr;//not implemented
    else if (generator_name == "ge_gen" || generator_name == "ge_gen_simplified")
        gen = new GETreeParameters();
    else if (generator_name == "weber_penn_gen")
        gen = new WeberPennParametersNative();
    else
        gen = nullptr;

    return gen; 
}