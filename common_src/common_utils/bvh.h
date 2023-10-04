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

    struct BVH_node
    {
        int bbox_idx = -1;
        int right_idx = -1;
        int left_idx = -1;
    };
    bool simple_list = true;
    std::vector<std::pair<AABB, uint64_t>> obj_bboxes;//<bbox, tag>

    std::vector<BVH_node> nodes;
    int root_node_idx = -1;
    int added_boxes_cnt = 0;
    int removed_boxes_cnt = 0;
    //different bboxes can have the same tag is they will always be added and deleted together
    void rebuild();
    void clear();
    int add_node_rec(std::vector<int> &boxes);
    void add_bboxes(std::vector<AABB> &bboxes, uint64_t tag);
    void remove_bboxes(int64_t tag);
    void remove_bboxes_iterate(int64_t tag, std::function<void(const std::pair<AABB, uint64_t> &)> func);
    void iterate_over_intersected_bboxes(AABB bbox, std::function<void(const std::pair<AABB, uint64_t> &)> func, bool debug = false) const;
    bool contains(glm::vec3 point) const;
};