#pragma once

#include "tree_generators/GE_generator.h"
#include "tree_generators/simple_generator.h"
#include "tree_generators/simpliest_generator.h"
#include "tree_generators/load_tree_structure.h"
#include "tree_generators/GE_generator_simplified.h"
#include "tree_generators/weber_penn_generator.h"

AbstractTreeGenerator *get_generator(const std::string &generator_name);
ParameterSet *get_default_parameters(const std::string &generator_name);