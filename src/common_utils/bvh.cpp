#include "bvh.h"
#include <stack>
#include <algorithm>

void BVH::add_bboxes(std::vector<AABB> &bboxes, uint64_t tag)
{
    added_boxes_cnt += bboxes.size();
    for (auto &box : bboxes)
        obj_bboxes.push_back(std::pair<AABB, uint64_t>(box, tag));
}

void BVH::remove_bboxes(int64_t tag)
{
  for (auto &p : obj_bboxes)
  {
    if (p.second == tag)
    {
      p.second = 0;
      removed_boxes_cnt++;
    }
  }
}

void BVH::remove_bboxes_iterate(int64_t tag, std::function<void(const std::pair<AABB, uint64_t> &)> func)
{
  for (auto &p : obj_bboxes)
  {
    if (p.second == tag)
    {
      func(p);
      p.second = 0;
      removed_boxes_cnt++;
    }
  }
}

void BVH::iterate_over_intersected_bboxes(AABB bbox, std::function<void(const std::pair<AABB, uint64_t> &)> func, bool debug) const
{
    if (simple_list)
    {
        for (auto &p : obj_bboxes)
        {
            if (p.second && bbox.intersects(p.first))
                func(p);
        }
    }
    else
    {
        if (root_node_idx >= 0 && obj_bboxes[nodes[root_node_idx].bbox_idx].first.intersects(bbox))
        {
            std::stack<int> stack;
            stack.push(root_node_idx);
            while (!stack.empty())
            {
                auto &node = nodes[stack.top()];
                stack.pop();
                
                if ((node.right_idx < 0 && node.left_idx < 0) || debug)
                    func(obj_bboxes[node.bbox_idx]);

                if (node.right_idx >= 0 && obj_bboxes[nodes[node.right_idx].bbox_idx].first.intersects(bbox))
                    stack.push(node.right_idx);
                if (node.left_idx >= 0 && obj_bboxes[nodes[node.left_idx].bbox_idx].first.intersects(bbox))
                    stack.push(node.left_idx);
            }
        }
    }
}

bool BVH::contains(glm::vec3 point) const
{
    if (simple_list)
    {
        for (auto &p : obj_bboxes)
        {
            if (p.second && p.first.contains(point))
                return true;
        }
    }
    else
    {
        if (root_node_idx >= 0 && obj_bboxes[nodes[root_node_idx].bbox_idx].first.contains(point))
        {
            std::stack<int> stack;
            stack.push(root_node_idx);
            while (!stack.empty())
            {
                auto &node = nodes[stack.top()];
                stack.pop();
                
                if (node.right_idx < 0 && node.left_idx < 0)
                    return true;

                if (node.right_idx >= 0 && obj_bboxes[nodes[node.right_idx].bbox_idx].first.contains(point))
                    stack.push(node.right_idx);
                if (node.left_idx >= 0 && obj_bboxes[nodes[node.left_idx].bbox_idx].first.contains(point))
                    stack.push(node.left_idx);
            }
        }
    }
    return false;
}

int BVH::add_node_rec(std::vector<int> &boxes)
{
    if (boxes.size() == 0)
        return -1;
    nodes.emplace_back();
    int node_id = nodes.size() - 1;
    //auto &node = nodes.back();
    if (boxes.size() == 1)
    {
        nodes[node_id].bbox_idx = boxes.back();
    }
    else
    {
        glm::vec3 minp = glm::vec3(1e9,1e9,1e9);
        glm::vec3 maxp = glm::vec3(-1e9,-1e9,-1e9);
        for (int box_id : boxes)
        {
            minp = min(minp, obj_bboxes[box_id].first.min_pos);
            maxp = max(maxp, obj_bboxes[box_id].first.max_pos);
        }
        obj_bboxes.push_back(std::pair<AABB, uint64_t>(AABB(minp, maxp), 0));
        nodes[node_id].bbox_idx = obj_bboxes.size() - 1;

        int best_axis = (maxp.x - minp.x > maxp.z - minp.z) ? 0 : 2;
        std::sort(boxes.begin(), boxes.end(), [&](const int & a, const int & b) -> bool{    
            return obj_bboxes[a].first.min_pos[best_axis] < obj_bboxes[b].first.min_pos[best_axis];});
        float split = obj_bboxes[boxes[boxes.size()/2]].first.min_pos[best_axis];
        std::vector<int> right_idx = std::vector<int>(boxes.begin(), boxes.begin() + boxes.size()/2);
        std::vector<int> left_idx = std::vector<int>(boxes.begin() + boxes.size()/2, boxes.end());
        /*
        for (int box_id : boxes)
        {
            if (obj_bboxes[box_id].first.min_pos[best_axis] < split)
                left_idx.push_back(box_id);
            if (obj_bboxes[box_id].first.max_pos[best_axis] > split)
                right_idx.push_back(box_id);
        }
                                debug("boxes {");
    for (int box_id : boxes)
    {
        debug("%d ", box_id);
    }
    debug("}\n");
                        debug("right {");
    for (int box_id : right_idx)
    {
        debug("%d ", box_id);
    }
    debug("}\n");
                                debug("left {");
    for (int box_id : left_idx)
    {
        debug("%d ", box_id);
    }
    debug("}\n");
    */
        if (!right_idx.empty())
        {

            if (right_idx.size() == boxes.size())
            {
                logerr("right stuck with size %d split %f axis %d", right_idx.size(), split, best_axis);
                for (int box_id : boxes)
                {
                    AABB &box = obj_bboxes[box_id].first;
                    logerr("box %d %f %f - %f %f", box_id, box.min_pos.x, box.min_pos.z, box.max_pos.x, box.max_pos.z);
                }
            }
            nodes[node_id].right_idx = add_node_rec(right_idx);
        }
        if (!left_idx.empty())
        {

            if (left_idx.size() == boxes.size())
            {
                logerr("left stuck with size %d split %f axis %d", left_idx.size(), split, best_axis);
                for (int box_id : boxes)
                {
                    AABB &box = obj_bboxes[box_id].first;
                    logerr("box %d %f %f - %f %f", box_id, box.min_pos.x, box.min_pos.z, box.max_pos.x, box.max_pos.z);
                }
            }
            nodes[node_id].left_idx = add_node_rec(left_idx);
        }
    }
    //logerr("created node %d %d %d total nodes/bboxes %d %d",nodes[node_id].bbox_idx,nodes[node_id].left_idx,nodes[node_id].right_idx,
    //                    nodes.size(),obj_bboxes.size());
    return node_id;
}
void BVH::clear()
{
    for (auto &p : obj_bboxes)
    {
        p.second = 0;
    }
    rebuild();
}
void BVH::rebuild()
{
  auto old_boxes = obj_bboxes;
  obj_bboxes.clear();
  for (auto &p : old_boxes)
  {
    if (p.second > 0)
      obj_bboxes.push_back(p);
  }
  removed_boxes_cnt = 0;
  nodes = {};
  std::vector<int> boxes = std::vector<int>(obj_bboxes.size(), 0);
  for (int i = 0; i < obj_bboxes.size(); i++)
    boxes[i] = i;
  root_node_idx = add_node_rec(boxes);
  added_boxes_cnt = 0;
}