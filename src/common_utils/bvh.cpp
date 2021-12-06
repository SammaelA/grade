#include "bvh.h"

void BVH::add_bboxes(std::vector<AABB> &bboxes, uint64_t tag)
{
    if (simple_list)
    {
        for (auto &box : bboxes)
            obj_bboxes.push_back(std::pair<AABB, uint64_t>(box, tag));
    }
}

void BVH::remove_bboxes(int64_t tag)
{
    if (simple_list)
    {
        auto it = obj_bboxes.begin();
        while (it != obj_bboxes.end())
        {
            if (it->second == tag)
                it = obj_bboxes.erase(it);
            else
                it++;
        }
    }
}

void BVH::iterate_over_intersected_bboxes(AABB bbox, std::function<void(const std::pair<AABB, uint64_t> &)> func)
{
    if (simple_list)
    {
        for (auto &p : obj_bboxes)
        {
            if (bbox.intersects(p.first))
                func(p);
        }
    }
}

bool BVH::contains(glm::vec3 point)
{
    if (simple_list)
    {
        for (auto &p : obj_bboxes)
        {
            if (p.first.contains(point))
                return true;
        }
    }

    return false;
}