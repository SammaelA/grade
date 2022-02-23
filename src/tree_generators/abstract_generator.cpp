#include "abstract_generator.h"

int AbstractTreeGenerator::joints_limit = 1000000;
std::atomic<int> AbstractTreeGenerator::branch_next_id(0);
std::atomic<int> AbstractTreeGenerator::tree_next_id(0); 