#include "GA_utils.h"
using namespace glm;
void DebugGraph::add_node(Node n)
{
    nodes.push_back(n);
}

void DebugGraph::add_edge(Edge e)
{
    edges.push_back(e);
}

void DebugGraph::clear()
{
    nodes.clear();
    edges.clear();
}

void DebugGraph::save_as_image(std::string name, int pix_x, int pix_y)
{
    vec2 min_coord = vec2(1e9,1e9);
    vec2 max_coord = vec2(-1e9,-1e9);

    for (auto &n : nodes)
    {
        min_coord.x = MIN(min_coord.x, n.pos.x);
        min_coord.y = MIN(min_coord.y, n.pos.y);
        max_coord.x = MAX(max_coord.x, n.pos.x);
        max_coord.y = MAX(max_coord.y, n.pos.y);
    } 
    if (length(max_coord - min_coord) < 1e-6)
        return;
    vec2 center = 0.5f*(max_coord + min_coord);
    vec2 sz_2 = 0.5f*(max_coord - min_coord);
    min_coord = center - 1.2f*sz_2;
    max_coord = center + 1.2f*sz_2;
    vec2 sz = (max_coord - min_coord);
    int pr_x = pix_y/sz.y * sz.x;
    int pr_y = pix_x/sz.x *sz.y;
    pix_x = MAX(pix_x, pr_x);
    pix_y = MAX(pix_y, pr_y);

    unsigned char *image = new unsigned char[4*pix_x*pix_y];
    unsigned char def[4] = {0,0,0,255};
    std::fill_n((int32_t *)image, pix_x*pix_y, *((int32_t*)def));
    vec2 sizes = vec2(pix_x, pix_y);

    for (Edge &e : edges)
    {
        Node &n1 = nodes[e.node_from];
        Node &n2 = nodes[e.node_to];
        ivec2 pix_n1 = ivec2(sizes*(n1.pos - min_coord)/(max_coord - min_coord));
        ivec2 pix_n2 = ivec2(sizes*(n2.pos - min_coord)/(max_coord - min_coord));
        vec3 line;//A,B,C in Ax + By + C = 0
        if (pix_n1.x == pix_n2.x)
        {
            if (pix_n1.y == pix_n2.y)
                continue;
            else
            {
                line = vec3(1, 0, pix_n1.x);
            }
        }
        else
        {
            line = vec3(pix_n2.y - pix_n1.y, -(pix_n2.x - pix_n1.x),-pix_n1.x*(pix_n2.y - pix_n1.y) + pix_n1.y*(pix_n2.x - pix_n1.x));
        }
        float ld = 1/sqrt(SQR(line.x) + SQR(line.y));
        glm::vec2 th = sizes*e.thickness/(max_coord - min_coord);
        for (int y = MIN(pix_n1.y, pix_n2.y); y <= MAX(pix_n1.y, pix_n2.y);y++)
        {
            for (int x = MIN(pix_n1.x, pix_n2.x); x <= MAX(pix_n1.x, pix_n2.x);x++)
            {
                float d = abs(line.x*x + line.y*y + line.z)*ld;
                if (d < 0.5*(th.x + th.y))
                {
                    //logerr("write edge %d %d",x, y);
                    image[4*(y*pix_x + x)  ] = 255*e.color.x;
                    image[4*(y*pix_x + x)+1] = 255*e.color.y;
                    image[4*(y*pix_x + x)+2] = 255*e.color.z;
                    image[4*(y*pix_x + x)+3] = 255;  
                }
            }            
        }
    }
    for (Node &n : nodes)
    {
        ivec2 pix_pos = ivec2(sizes*(n.pos - min_coord)/(max_coord - min_coord));
        ivec2 pix_r = max(vec2(1,1), n.radius* sizes/(max_coord - min_coord));
        vec2 pix_ir2 = vec2(1.0f/SQR(pix_r.x), 1.0f/SQR(pix_r.y));

        for (int y = MAX(0, pix_pos.y - pix_r.y); y < MIN(pix_y, pix_pos.y + pix_r.y); y++)
        {
            for (int x = MAX(0, pix_pos.x - pix_r.x); x < MIN(pix_x, pix_pos.x + pix_r.x); x++)
            {
                if (SQR(x - pix_pos.x)*pix_ir2.x + SQR(y - pix_pos.y)*pix_ir2.y <= 1)
                {
                    //logerr("write %d %d",x, y);
                    image[4*(y*pix_x + x)  ] = 255*n.color.x;
                    image[4*(y*pix_x + x)+1] = 255*n.color.y;
                    image[4*(y*pix_x + x)+2] = 255*n.color.z;
                    image[4*(y*pix_x + x)+3] = 255;
                }
            }
        }
    }

    textureManager.save_png_raw(image, pix_x, pix_y, 4, name);
    delete[] image;
}