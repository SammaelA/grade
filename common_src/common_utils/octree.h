#pragma once 
#include "bbox.h"
#include <list>
#include <functional>
#include <vector>

/*
Octree built from a point cloud. It can be used to find 
closest neighbours or iterate over all points in some region.
It's implementation is very basic and not really efficient,
but still better than brute force search.
*/
class PointCloudOctree
{
public:
    static const int MAX_NODES = 8; //maximum points in leaf
    PointCloudOctree();
    PointCloudOctree(AABB _box);
    void create(AABB _box);
    ~PointCloudOctree(){clear();}

    void insert(float3 &pos) {root.insert(pos);}
    void clear() {root.clear();}
    void insert_vector(std::vector<float3> &positions);
    void apply_to_neighbours_AABB(AABB &_box, std::function<void(float3 &)> func)
    {
        root.apply_to_neighbours_AABB(_box, func);
    }
    void apply_to_neighbours_sphere(AABB &_box, float r, float3 &center, std::function<void(float3 &)> func)
    {
        root.apply_to_neighbours_sphere(_box, r, center, func);
    }
    void remove_in_sphere(AABB &_box, float r, float3 &center)
    {
        root.remove_in_sphere(_box, r, center);
    }
private:
    #define GetN(x,y,z) (z*4 + y*2 + x)
    struct Node
    {
        AABB box;
        std::list<float3> points;
        Node *child_nodes[8];
        Node(AABB _box): box(_box) {for (int i=0;i<8;i++) child_nodes[i] = nullptr;}
        ~Node();

        void clear();
        void subdivide_median();
        void subdivide_half();
        bool insert(float3 &pos);
        void apply_to_neighbours_AABB(AABB &_box, std::function<void(float3 &)> func);
        void apply_to_neighbours_sphere(AABB &_box, float r, float3 &center, std::function<void(float3 &)> func);
        void remove_in_sphere(AABB &_box, float r, float3 &center);
    };

    static Node *new_node(AABB &_box);
    static void delete_node(Node *n);

    Node root;
};