#pragma once

#include "tree_generators/GE_generator.h"
#include "tree_generators/python_tree_gen.h"
#include "tree_generators/weber_penn_parameters.h"
#include "tree_generators/simple_generator.h"
#include "tree_generators/simpliest_generator.h"
#include "tree_generators/generated_tree.h"
#include "tree_generators/GE_generator_simplified.h"
#include "tree_generators/weber_penn_generator.h"

AbstractTreeGenerator *get_generator(std::string &generator_name);
ParametersSet *get_default_parameters(std::string &generator_name);