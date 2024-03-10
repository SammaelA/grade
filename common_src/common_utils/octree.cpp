#include "octree.h"
#include "common_utils/utility.h"

Octree::Node *Octree::new_node(AABB &_box) 
{
    return new Node(_box);
}

void Octree::delete_node(Node *n)
{
    if (n) delete n;
}

void Octree::insert_vector(std::vector<float3> &positions)
{
    for (auto &p : positions)
        root.insert(p);
}

Octree::Octree() : root(AABB(float3(0,0,0),float3(0,0,0))) {};
Octree::Octree(AABB _box) : root(_box) {};
void Octree::create(AABB _box)
{
    root.clear();
    root.box = _box;
}

bool Octree::Node::insert(float3 &pos)
{
    //logerr("insert (%f %f %f) in [%f %f %f]-[%f %f %f]", pos.x, pos.y, pos.z,
    //       box.min_pos.x,box.min_pos.y,box.min_pos.z,
    //       box.max_pos.x, box.max_pos.y, box.max_pos.z);
    if (!box.contains(pos))
        return false;
    
    if (child_nodes[0])
    {
        for (int i=0;i<8;i++) 
        {
            if (child_nodes[i]->insert(pos))
                return true;
        }
        /*
        logerr("error: malformed octree node %f %f %f",pos.x,pos.y,pos.z);
        for (int i=0;i<8;i++) 
        {
            logerr("ch node %f %f %f -- %f %f %f",child_nodes[i]->box.min_pos.x,child_nodes[i]->box.min_pos.y,
            child_nodes[i]->box.min_pos.z,child_nodes[i]->box.max_pos.x,child_nodes[i]->box.max_pos.y,
            child_nodes[i]->box.max_pos.z);
        }*/
        return false;
    }
    else if (points.size() < Octree::MAX_NODES)
    {
        points.push_back(pos);
        return true;
    }
    else
    {
        subdivide_half();
        
        for (int i=0;i<8;i++) 
        {
            if (child_nodes[i]->insert(pos))
                return true;
        }
    }
    return false;
}
void Octree::Node::subdivide_half()
{
    for (int x = 0;x<=1;x++)
    {
       for (int y = 0;y<=1;y++)
        {
            for (int z = 0;z<=1;z++)
            {
                float3 sizes = box.max_pos - box.min_pos;
                float3 new_min = box.min_pos + float3(x*0.5f*sizes.x, y*0.5f*sizes.y, z*0.5f*sizes.z);
                AABB _box(new_min, new_min + 0.5f*sizes);
                child_nodes[GetN(x,y,z)] = Octree::new_node(_box);
            }
        } 
    }
    for (auto &p : points)
    {
        for (int i=0;i<8;i++) 
        {
            if (child_nodes[i]->insert(p))
                break;
        }

    }
    points = {};
}
void Octree::Node::apply_to_neighbours_AABB(AABB &_box, std::function<void(float3 &)> func)
{
    if (!box.intersects(_box))
        return;
    if (child_nodes[0])
    {
        for (int i=0;i<8;i++) 
        {
            child_nodes[i]->apply_to_neighbours_AABB(_box, func);
        }
    }
    else
    {
        for (auto &p : points)
            func(p);
    }
}
void Octree::Node::apply_to_neighbours_sphere(AABB &_box, float r, float3 &center, 
                                              std::function<void(float3 &)> func)
{
    if (!box.intersects(_box))
        return;
    //logerr("apply (%f %f %f) %f to [%f %f %f]-[%f %f %f]", center.x, center.y, center.z, r,
    //      box.min_pos.x,box.min_pos.y,box.min_pos.z,
    //      box.max_pos.x, box.max_pos.y, box.max_pos.z);
    if (child_nodes[0])
    {
        for (int i=0;i<8;i++) 
        {
            child_nodes[i]->apply_to_neighbours_sphere(_box, r, center, func);
        }
    }
    else
    {
        float r_sq = r*r;
        for (auto &p : points)
        {
            if (dot(p - center, p - center) <= r_sq)
                func(p);
        }
    }
}
void Octree::Node::remove_in_sphere(AABB &_box, float r, float3 &center)
{
    if (!box.intersects(_box))
        return;
    if (child_nodes[0])
    {
        for (int i=0;i<8;i++) 
        {
            child_nodes[i]->remove_in_sphere(_box, r, center);
        }
    }
    else
    {
        float r_sq = r*r;
        auto it = points.begin();
        while (it != points.end())
        {
            if (dot(*it - center, *it - center) <= r_sq)
                it = points.erase(it);
            else 
                it++;
        }
    }
}

Octree::Node::~Node()
{
    clear();
}

void Octree::Node::clear()
{
    points.clear();
    if (child_nodes[0])
    {
        for (int i=0;i<8;i++)
        {
            child_nodes[i]->clear();
            Octree::delete_node(child_nodes[i]);
            child_nodes[i] = nullptr;
        }
    }
}