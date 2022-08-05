#pragma once

#include "load_tree_structure.h"

class PythonTreeGen : public TreeLoaderBlk
{
public:
  virtual void plant_tree(glm::vec3 pos, TreeTypeData *type) override;
  virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
private:
  std::vector<std::pair<glm::vec3, TreeTypeData *>> tree_saplings;
};