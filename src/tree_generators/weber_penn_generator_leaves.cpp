#include "weber_penn_generator.h"
#include "turtle.h"

using glm::vec3;
using glm::vec4;
using glm::mat4;

void leaves(int n, std::vector<glm::vec3> &out_verts, std::vector<std::vector<int>> &out_indicies);
void blossom(int n, std::vector<glm::vec3> &out_verts, std::vector<std::vector<int>> &out_indicies);

void WeberPennGenerator::Leaf::get_shape(int leaf_type, float g_scale, float scale, float scale_x, BaseLeafMesh &blm)
{
    if (leaf_type < 0) //blossom
    {
        if (leaf_type < -3)
            leaf_type = -1; 
        blossom(abs(leaf_type + 1), blm.verts, blm.indicies);
    }
    else
    {
        if (leaf_type < 1 || leaf_type > 10)
            leaf_type = 8;
        leaves(leaf_type - 1, blm.verts, blm.indicies);
    }
    for (auto &vert : blm.verts)
    {
        vert *= scale * g_scale;
        vert.x *= scale_x;
    }
}
void WeberPennGenerator::Leaf::get_mesh(float bend, BaseLeafMesh &base_shape, int index, std::vector<glm::vec3> &out_verts, 
                                        std::vector<std::vector<int>> &out_indicies)
{
    auto trf = to_track_quat_ZY(direction);
    out_verts = base_shape.verts;
    out_indicies = base_shape.indicies;
    for (auto &v : out_verts)
    {
        v = trf*v;
    }
    LiteMath::quat tr1 = LiteMath::quat(1,0,0,0);
    LiteMath::quat tr2 = LiteMath::quat(1,0,0,0);
    if (bend > 0)
    {
        calc_bend_trf(bend, tr1, tr2);
        for (auto &v : out_verts)
            v = tr1*v;
    }
    for (auto &v : out_verts)
    {
        v += position;
    }
    index *= out_verts.size();

    //Not sure that it is needed
    for (auto &inds : out_indicies)
    {
        for (auto &ind : inds)
            ind += index;
    }
}
void WeberPennGenerator::Leaf::calc_bend_trf(float bend, LiteMath::quat &bend_trf_1, LiteMath::quat &bend_trf_2)
{
    glm::vec3 normal = glm::cross(direction, right);
    float theta_pos = atan2(position.y, position.x);
    float theta_bend = theta_pos - atan2(normal.y, normal.x);
    bend_trf_1 = LiteMath::angleAxis(theta_bend * bend, glm::vec3(0,0,1));
    //I think this is what the paper says but the second transform just looks stupid
    //so we just ignore it above

    direction = bend_trf_1*direction;
    right = bend_trf_1*right;
    normal = glm::cross(direction, right);

    //not sure it is equal to phi_bend = normal.declination()
    float phi_bend = LiteMath::to_radians(declination(normal));
    if (phi_bend > PI / 2)
        phi_bend = phi_bend - PI;
    bend_trf_2 = LiteMath::angleAxis(phi_bend * bend, right);
}

void leaves(int n, std::vector<glm::vec3> &out_verts, std::vector<std::vector<int>> &out_indicies)
{
    //we cannot handle other shapers that square, so no need to create more vertices
    n = 8;
    if (n == 0)
    {
        //1 = ovate
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.15, 0, 0.15),
                glm::vec3(0.25, 0, 0.3),
                glm::vec3(0.2, 0, 0.6),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.2, 0, 0.6),
                glm::vec3(-0.25, 0, 0.3),
                glm::vec3(-0.15, 0, 0.15),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {{0, 1, 9, 10}, {1, 2, 3, 4}, {4, 5, 6}, {6, 7, 8, 9}, {4, 6, 9, 1}};    
    }
    else if (n == 1)
    {
        //2 = linear
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.1, 0, 0.15),
                glm::vec3(0.1, 0, 0.95),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.1, 0, 0.95),
                glm::vec3(-0.1, 0, 0.15),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {{0, 1, 7, 8}, {1, 2, 3}, {3, 4, 5}, {5, 6, 7}, {1, 3, 5, 7}};    
    }
    else if (n == 2)
    {
        //3 = cordate
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.01, 0, 0.2),
                glm::vec3(0.2, 0, 0.1),
                glm::vec3(0.35, 0, 0.35),
                glm::vec3(0.25, 0, 0.6),
                glm::vec3(0.1, 0, 0.8),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.1, 0, 0.8),
                glm::vec3(-0.25, 0, 0.6),
                glm::vec3(-0.35, 0, 0.35),
                glm::vec3(-0.2, 0, 0.1),
                glm::vec3(-0.01, 0, 0.2),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {
                {0, 1, 11, 12},
                {1, 2, 3, 4},
                {11, 10, 9, 8},
                {11, 1, 4, 8},
                {8, 7, 6, 5, 4}
            };    
    }
    else if (n == 3)
    {
        //4 = maple
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.25, 0, 0.07),
                glm::vec3(0.2, 0, 0.18),
                glm::vec3(0.5, 0, 0.37),
                glm::vec3(0.43, 0, 0.4),
                glm::vec3(0.45, 0, 0.58),
                glm::vec3(0.3, 0, 0.57),
                glm::vec3(0.27, 0, 0.67),
                glm::vec3(0.11, 0, 0.52),
                glm::vec3(0.2, 0, 0.82),
                glm::vec3(0.08, 0, 0.77),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.08, 0, 0.77),
                glm::vec3(-0.2, 0, 0.82),
                glm::vec3(-0.11, 0, 0.52),
                glm::vec3(-0.27, 0, 0.67),
                glm::vec3(-0.3, 0, 0.57),
                glm::vec3(-0.45, 0, 0.58),
                glm::vec3(-0.43, 0, 0.4),
                glm::vec3(-0.5, 0, 0.37),
                glm::vec3(-0.2, 0, 0.18),
                glm::vec3(-0.25, 0, 0.07),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {
                {0, 1, 23, 24},
                {1, 2, 3, 4, 5},
                {23, 22, 21, 20, 19},
                {1, 5, 6, 7, 8},
                {23, 19, 18, 17, 16},
                {1, 8, 9, 10, 11},
                {23, 16, 15, 14, 13},
                {1, 11, 12, 13, 23}
            };    
    }
    else if (n == 4)
    {
        //5 = palmate
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.25, 0, 0.1),
                glm::vec3(0.5, 0, 0.3),
                glm::vec3(0.2, 0, 0.45),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.2, 0, 0.45),
                glm::vec3(-0.5, 0, 0.3),
                glm::vec3(-0.25, 0, 0.1),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {{0, 1, 9, 10}, {1, 2, 3, 4}, {1, 4, 5, 6, 9}, {9, 8, 7, 6}};    
    }
    else if (n == 5)
    {
        //6 = spiky oak
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.16, 0, 0.17),
                glm::vec3(0.11, 0, 0.2),
                glm::vec3(0.23, 0, 0.33),
                glm::vec3(0.15, 0, 0.34),
                glm::vec3(0.32, 0, 0.55),
                glm::vec3(0.16, 0, 0.5),
                glm::vec3(0.27, 0, 0.75),
                glm::vec3(0.11, 0, 0.7),
                glm::vec3(0.18, 0, 0.9),
                glm::vec3(0.07, 0, 0.86),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.07, 0, 0.86),
                glm::vec3(-0.18, 0, 0.9),
                glm::vec3(-0.11, 0, 0.7),
                glm::vec3(-0.27, 0, 0.75),
                glm::vec3(-0.16, 0, 0.5),
                glm::vec3(-0.32, 0, 0.55),
                glm::vec3(-0.15, 0, 0.34),
                glm::vec3(-0.23, 0, 0.33),
                glm::vec3(-0.11, 0, 0.2),
                glm::vec3(-0.16, 0, 0.17),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {
                {0, 1, 23, 24},
                {1, 2, 3},
                {3, 4, 5},
                {5, 6, 7},
                {7, 8, 9},
                {9, 10, 11},
                {1, 3, 5, 7, 9, 11, 12, 13, 15, 17, 19, 21, 23},
                {23, 22, 21},
                {21, 20, 19},
                {19, 18, 17},
                {17, 16, 15},
                {15, 14, 13}
            };    
    }
    else if (n == 6)
    {
        //7 = round oak
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.11, 0, 0.16),
                glm::vec3(0.11, 0, 0.2),
                glm::vec3(0.22, 0, 0.26),
                glm::vec3(0.23, 0, 0.32),
                glm::vec3(0.15, 0, 0.34),
                glm::vec3(0.25, 0, 0.45),
                glm::vec3(0.23, 0, 0.53),
                glm::vec3(0.16, 0, 0.5),
                glm::vec3(0.23, 0, 0.64),
                glm::vec3(0.2, 0, 0.72),
                glm::vec3(0.11, 0, 0.7),
                glm::vec3(0.16, 0, 0.83),
                glm::vec3(0.12, 0, 0.87),
                glm::vec3(0.06, 0, 0.85),
                glm::vec3(0.07, 0, 0.95),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.07, 0, 0.95),
                glm::vec3(-0.06, 0, 0.85),
                glm::vec3(-0.12, 0, 0.87),
                glm::vec3(-0.16, 0, 0.83),
                glm::vec3(-0.11, 0, 0.7),
                glm::vec3(-0.2, 0, 0.72),
                glm::vec3(-0.23, 0, 0.64),
                glm::vec3(-0.16, 0, 0.5),
                glm::vec3(-0.23, 0, 0.53),
                glm::vec3(-0.25, 0, 0.45),
                glm::vec3(-0.15, 0, 0.34),
                glm::vec3(-0.23, 0, 0.32),
                glm::vec3(-0.22, 0, 0.26),
                glm::vec3(-0.11, 0, 0.2),
                glm::vec3(-0.11, 0, 0.16),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {
                {0, 1, 33, 34},
                {1, 2, 3},
                {3, 4, 5, 6},
                {6, 7, 8, 9},
                {9, 10, 11, 12},
                {12, 13, 14, 15},
                {15, 16, 17},
                {1, 3, 6, 9, 12, 15, 17, 19, 22, 25, 28, 31, 33},
                {33, 32, 31},
                {31, 30, 29, 28},
                {28, 27, 26, 25},
                {25, 24, 23, 22},
                {22, 21, 20, 19},
                {19, 18, 17}
            };    
    }
    else if (n == 7)
    {
        //8 = elliptic (default)
        out_verts = {
                glm::vec3(0.005, 0, 0),
                glm::vec3(0.005, 0, 0.1),
                glm::vec3(0.15, 0, 0.2),
                glm::vec3(0.25, 0, 0.45),
                glm::vec3(0.2, 0, 0.75),
                glm::vec3(0, 0, 1),
                glm::vec3(-0.2, 0, 0.75),
                glm::vec3(-0.25, 0, 0.45),
                glm::vec3(-0.15, 0, 0.2),
                glm::vec3(-0.005, 0, 0.1),
                glm::vec3(-0.005, 0, 0)
            };
        out_indicies = {{0, 1, 9, 10}, {1, 2, 3, 4}, {4, 5, 6}, {6, 7, 8, 9}, {4, 6, 9, 1}};    
    }
    else if (n == 8)
    {
        //9 = rectangle
        out_verts = {
                glm::vec3(-0.5, 0, 0),
                glm::vec3(-0.5, 0, 1),
                glm::vec3(0.5, 0, 1),
                glm::vec3(0.5, 0, 0)
            };
        out_indicies = {{0, 1, 2, 3}};    
    }
    else if (n == 9)
    {
        if (n != 9)
        {
            logerr("wrong leaf type %d", n);
        }
        //10 = triangle
        out_verts = {glm::vec3(-0.5, 0, 0), glm::vec3(0, 0, 1), glm::vec3(0.5, 0, 0)};
        out_indicies = {{0, 1, 2}};    
    }
}
void blossom(int n, std::vector<glm::vec3> &out_verts, std::vector<std::vector<int>> &out_indicies)
{
    //we cannot handle other shapers that square, so no need to create more vertices
    n = 8;
    leaves(8, out_verts, out_indicies);
    return;
    if (n == 0)
    {
        //1 = cherry
        out_verts = {
                glm::vec3(0, 0, 0),
                glm::vec3(0.33, 0.45, 0.45),
                glm::vec3(0.25, 0.6, 0.6),
                glm::vec3(0, 0.7, 0.7),
                glm::vec3(-0.25, 0.6, 0.6),
                glm::vec3(-0.33, 0.45, 0.45),
                glm::vec3(0.49, 0.42, 0.6),
                glm::vec3(0.67, 0.22, 0.7),
                glm::vec3(0.65, -0.05, 0.6),
                glm::vec3(0.53, -0.17, 0.45),
                glm::vec3(0.55, -0.33, 0.6),
                glm::vec3(0.41, -0.57, 0.7),
                glm::vec3(0.15, -0.63, 0.6),
                glm::vec3(0, -0.55, 0.45),
                glm::vec3(-0.15, -0.63, 0.6),
                glm::vec3(-0.41, -0.57, 0.7),
                glm::vec3(-0.55, -0.33, 0.6),
                glm::vec3(-0.53, -0.17, 0.45),
                glm::vec3(-0.65, -0.05, 0.6),
                glm::vec3(-0.67, 0.22, 0.7),
                glm::vec3(-0.49, 0.42, 0.6)
            };
        out_indicies = {
                {0, 1, 2, 3},
                {0, 3, 4, 5},
                {0, 1, 6, 7},
                {0, 7, 8, 9},
                {0, 9, 10, 11},
                {0, 11, 12, 13},
                {0, 13, 14, 15},
                {0, 15, 16, 17},
                {0, 17, 18, 19},
                {0, 19, 20, 5}
            };    
    }
    else if (n == 1)
    {
        //2 = orange
        out_verts = {
                glm::vec3(0, 0, 0),
                glm::vec3(-0.055, 0.165, 0.11),
                glm::vec3(-0.125, 0.56, 0.365),
                glm::vec3(0, 0.7, 0.45),
                glm::vec3(0.125, 0.56, 0.365),
                glm::vec3(0.055, 0.165, 0.11),
                glm::vec3(0.14, 0.10, 0.11),
                glm::vec3(0.495, 0.29, 0.365),
                glm::vec3(0.665, 0.215, 0.45),
                glm::vec3(0.57, 0.055, 0.36),
                glm::vec3(0.175, 0, 0.11),
                glm::vec3(0.14, -0.1, 0.11),
                glm::vec3(0.43, -0.38, 0.365),
                glm::vec3(0.41, -0.565, 0.45),
                glm::vec3(0.23, -0.53, 0.365),
                glm::vec3(0.05, -0.165, 0.11),
                glm::vec3(-0.14, -0.1, 0.11),
                glm::vec3(-0.43, -0.38, 0.365),
                glm::vec3(-0.41, -0.565, 0.45),
                glm::vec3(-0.23, -0.53, 0.365),
                glm::vec3(-0.05, -0.165, 0.11),
                glm::vec3(-0.14, 0.10, 0.11),
                glm::vec3(-0.495, 0.29, 0.365),
                glm::vec3(-0.665, 0.215, 0.45),
                glm::vec3(-0.57, 0.055, 0.36),
                glm::vec3(-0.175, 0, 0.11),
                glm::vec3(0.1, -0.1, 0.4),
                glm::vec3(-0.1, -0.1, 0.4),
                glm::vec3(-0.1, 0.1, 0.4),
                glm::vec3(0.1, 0.1, 0.4)
            };
        out_indicies = {
                {0, 1, 2, 3},
                {0, 3, 4, 5},
                {0, 6, 7, 8},
                {0, 8, 9, 10},
                {0, 11, 12, 13},
                {0, 13, 14, 15},
                {0, 16, 17, 18},
                {0, 18, 19, 20},
                {0, 21, 22, 23},
                {0, 23, 24, 25},
                {0, 26, 27},
                {0, 27, 28},
                {0, 28, 29},
                {0, 29, 26}
            };    
    }
    else if (n == 2)
    {
        //3 = magnolia
        out_verts = {
                glm::vec3(0, 0, 0),
                glm::vec3(0.19, -0.19, 0.06),
                glm::vec3(0.19, -0.04, 0.06),
                glm::vec3(0.34, -0.11, 0.35),
                glm::vec3(0.3, -0.3, 0.6),
                glm::vec3(0.11, -0.34, 0.35),
                glm::vec3(0.04, -0.19, 0.06),
                glm::vec3(0.19, 0.19, 0.06),
                glm::vec3(0.19, 0.04, 0.06),
                glm::vec3(0.34, 0.11, 0.35),
                glm::vec3(0.3, 0.3, 0.6),
                glm::vec3(0.11, 0.34, 0.35),
                glm::vec3(0.04, 0.19, 0.06),
                glm::vec3(-0.19, -0.19, 0.06),
                glm::vec3(-0.19, -0.04, 0.06),
                glm::vec3(-0.34, -0.11, 0.35),
                glm::vec3(-0.3, -0.3, 0.6),
                glm::vec3(-0.11, -0.34, 0.35),
                glm::vec3(-0.04, -0.19, 0.06),
                glm::vec3(-0.19, 0.19, 0.06),
                glm::vec3(-0.19, 0.04, 0.06),
                glm::vec3(-0.34, 0.11, 0.35),
                glm::vec3(-0.3, 0.3, 0.6),
                glm::vec3(-0.11, 0.34, 0.35),
                glm::vec3(-0.04, 0.19, 0.06),
                glm::vec3(0, -0.39, 0.065),
                glm::vec3(0.15, -0.23, 0.065),
                glm::vec3(0.23, -0.46, 0.39),
                glm::vec3(0, -0.62, 0.65),
                glm::vec3(-0.23, -0.46, 0.39),
                glm::vec3(-0.15, -0.23, 0.065),
                glm::vec3(0, 0.39, 0.065),
                glm::vec3(0.15, 0.23, 0.065),
                glm::vec3(0.23, 0.46, 0.39),
                glm::vec3(0, 0.62, 0.65),
                glm::vec3(-0.23, 0.46, 0.39),
                glm::vec3(-0.15, 0.23, 0.065),
                glm::vec3(-0.39, 0, 0.065),
                glm::vec3(-0.23, 0.15, 0.065),
                glm::vec3(-0.46, 0.23, 0.39),
                glm::vec3(-0.62, 0, 0.65),
                glm::vec3(-0.46, -0.23, 0.39),
                glm::vec3(-0.23, -0.15, 0.065),
                glm::vec3(0.39, 0, 0.065),
                glm::vec3(0.23, 0.15, 0.065),
                glm::vec3(0.46, 0.23, 0.39),
                glm::vec3(0.62, 0, 0.65),
                glm::vec3(0.46, -0.23, 0.39),
                glm::vec3(0.23, -0.15, 0.065)
            };
        out_indicies = {
                {0, 1, 2},
                {1, 2, 3},
                {1, 3, 4},
                {1, 4, 5},
                {1, 5, 6},
                {1, 6, 0},
                {0, 7, 8},
                {7, 8, 9},
                {7, 9, 10},
                {7, 10, 11},
                {7, 11, 12},
                {7, 12, 0},
                {0, 13, 14},
                {13, 14, 15},
                {13, 15, 16},
                {13, 16, 17},
                {13, 17, 18},
                {13, 18, 0},
                {0, 19, 20},
                {19, 20, 21},
                {19, 21, 22},
                {19, 22, 23},
                {19, 23, 24},
                {19, 24, 0},
                {0, 25, 26},
                {25, 26, 27},
                {25, 27, 28},
                {25, 28, 29},
                {25, 29, 30},
                {25, 30, 0},
                {0, 31, 32},
                {31, 32, 33},
                {32, 33, 34},
                {31, 34, 35},
                {31, 35, 36},
                {31, 36, 0},
                {0, 37, 38},
                {37, 38, 39},
                {37, 39, 40},
                {37, 40, 41},
                {37, 41, 42},
                {37, 42, 0},
                {0, 43, 44},
                {43, 44, 45},
                {43, 45, 46},
                {43, 46, 47},
                {43, 47, 48},
                {43, 48, 0}
            };    
    }
    else 
    {
        logerr("wrong blossom type %d", n);
        leaves(8, out_verts, out_indicies);
    }
}
