#pragma once

#include "graphics_utils/volumetric_occlusion.h"
#include "abstract_generator.h"
#include "common_utils/parameter.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "GE_generator_parameters.h"
#include "GE_generator.h"
#include <vector>
#include <list>
#include <atomic>


class GETreeGeneratorSimplified : public GETreeGenerator
{
public:
    virtual bool iterate(LightVoxelsCube &voxels) override;
private:

    void prepare_nodes_and_space_colonization(Tree &t, Branch &b, GETreeParameters &params, 
                                              std::vector<GrowPoint> &growth_points,
                                              int max_growth_per_node);
    virtual bool find_best_pos(LightVoxelsCube &voxels, float r, glm::vec3 pos,
                               glm::vec3 dir, float angle,
                               glm::vec3 &best_pos, float &best_occ) override;
};