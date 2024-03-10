#pragma once
#include "impostor_similarity.h"
#include "tree_generators/abstract_generator.h"
#include <chrono>

class DebugGraph
{
public:
    struct Node
    {
        float2 pos;
        float3 color;
        float radius;
        Node(float2 p, float3 c, float r): pos(p), color(c), radius(r) {};
    };
    struct Edge
    {
        int node_from;
        int node_to;
        float thickness;
        float3 color;
        Edge(int from, int to, float th, float3 c = float3(0.5,0.5,0.5)): 
        node_from(from), node_to(to), thickness(th), color(c) {};
    };

    void add_node(Node n);
    void add_edge(Edge e);
    void clear();
    void save_as_image(std::string name, int pix_x, int pix_y);
private:

    std::vector<Node> nodes;
    std::vector<Edge> edges;
};