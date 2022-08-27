#pragma once

#include "core/tree.h"
#include "graphics_utils/volumetric_occlusion.h"
#include <vector>
#include <map>
#include "common_utils/blk.h"
#include "graphics_utils/billboard_cloud.h"

glm::vec3 canonical_bbox();
bool get_dedicated_bbox(Branch *branch, BBox &bbox);
void voxelize_original_branch(Branch *b, LightVoxelsCube *light, int level_to, float scale);
void set_occlusion(Branch *b, LightVoxelsCube *light);
    struct Answer
    {
        bool exact;
        float from;
        float to;
        Answer(bool ex, float fr, float t)
        {
            exact = ex;
            from = fr;
            to = t;
        }
        Answer() : Answer(false, 0, 1){};
        Answer(const Answer&) = default;
        Answer(Answer&&) = default;
        Answer& operator=(const Answer&) = default;
        Answer& operator=(Answer&&) = default;
        Answer operator*(const float mult)
        {
            return Answer(exact, MIN(from*mult, to*mult), MAX(from*mult, to*mult));
        }
        Answer operator+(const Answer &add)
        {
            return Answer(exact && add.exact, from + add.from, to + add.to);
        }
        Answer operator-(const Answer &sub)
        {
            return Answer(exact && sub.exact, from - sub.to, to -sub.from);
        }
    };