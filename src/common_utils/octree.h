#pragma once 
#include "bbox.h"
#include <list>
#include <functional>
#include <vector>

class Octree
{
public:
    static const int MAX_NODES = 8;
    Octree();
    Octree(AABB _box);
    void create(AABB _box);
    ~Octree(){clear();}

    void insert(glm::vec3 &pos) {root.insert(pos);}
    void clear() {root.clear();}
    void insert_vector(std::vector<glm::vec3> &positions);
    void apply_to_neighbours_AABB(AABB &_box, std::function<void(glm::vec3 &)> func)
    {
        root.apply_to_neighbours_AABB(_box, func);
    }
    void apply_to_neighbours_sphere(AABB &_box, float r, glm::vec3 &center, std::function<void(glm::vec3 &)> func)
    {
        root.apply_to_neighbours_sphere(_box, r, center, func);
    }
    void remove_in_sphere(AABB &_box, float r, glm::vec3 &center)
    {
        root.remove_in_sphere(_box, r, center);
    }
private:
    #define GetN(x,y,z) (z*4 + y*2 + x)
    struct Node
    {
        AABB box;
        std::list<glm::vec3> points;
        Node *child_nodes[8];
        Node(AABB _box): box(_box) {for (int i=0;i<8;i++) child_nodes[i] = nullptr;}
        ~Node();

        void clear();
        void subdivide_median();
        void subdivide_half();
        bool insert(glm::vec3 &pos);
        void apply_to_neighbours_AABB(AABB &_box, std::function<void(glm::vec3 &)> func);
        void apply_to_neighbours_sphere(AABB &_box, float r, glm::vec3 &center, std::function<void(glm::vec3 &)> func);
        void remove_in_sphere(AABB &_box, float r, glm::vec3 &center);
    };

    static Node *new_node(AABB &_box);
    static void delete_node(Node *n);

    Node root;
};