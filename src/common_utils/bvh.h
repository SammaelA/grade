#pragma once
#include "bbox.h"
#include <list>
#include <vector>
#include <functional>

struct BVH
{
    explicit BVH(bool _simple_list)
    {
        simple_list = _simple_list;
    }

    bool simple_list = true;
    std::list<std::pair<AABB, uint64_t>> obj_bboxes;//<bbox, tag>
    //different bboxes can have the same tag is they will always be added and deleted together
    void add_bboxes(std::vector<AABB> &bboxes, uint64_t tag);
    void remove_bboxes(int64_t tag);
    void iterate_over_intersected_bboxes(AABB bbox, std::function<void(const std::pair<AABB, uint64_t> &)> func);
    bool contains(glm::vec3 point);
};