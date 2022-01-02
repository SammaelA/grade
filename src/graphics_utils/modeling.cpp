#include "modeling.h"
#include "common_utils/utility.h"
#include "graphics_utils/texture_manager.h"

#define PI 3.14159265f
using namespace glm;

void Visualizer::leaf_to_model(Leaf &l, Model *m, float scale)
{
    if (l.edges.size() < 4)
    {
        return;
    }
    glm::vec3 a = l.edges[0];
    glm::vec3 b = l.edges[1];
    glm::vec3 c = l.edges[2];
    glm::vec3 n = glm::normalize(glm::cross(a - b, c - b));
    std::vector<float> tex_c{0, 0, 0, 1, 1, 1, 1, 0};
    int _b = m->positions.size() / 3;
    for (int i = 0; i < 4; i++)
    {
        glm::vec3 v = scale*(l.edges[i] - l.edges[0]) + l.edges[0];
        m->positions.push_back(v.x);
        m->positions.push_back(v.y);
        m->positions.push_back(v.z);
        m->normals.push_back(n.x);
        m->normals.push_back(n.y);
        m->normals.push_back(n.z);
        m->colors.push_back(tex_c[2 * i]);
        m->colors.push_back(tex_c[2 * i + 1]);
        m->colors.push_back(0);
        m->colors.push_back(1);
    }

    m->indices.push_back(_b);
    m->indices.push_back(_b + 1);
    m->indices.push_back(_b + 2);
    m->indices.push_back(_b + 2);
    m->indices.push_back(_b + 3);
    m->indices.push_back(_b);
}
void Visualizer::packed_leaf_to_model(PackedLeaf &l, Model *m, glm::vec2 tc_zw)
{
    if (l.edges.size() % 4 || l.edges.empty())
    {
        return;
    }
    for (int start_n = 0;start_n < l.edges.size();start_n += 4)
    {
        glm::vec3 a = l.edges[start_n + 0];
        glm::vec3 b = l.edges[start_n + 1];
        glm::vec3 c = l.edges[start_n + 2];
        glm::vec3 n = glm::normalize(glm::cross(a - b, c - b));
        std::vector<float> tex_c{0, 0, 0, 1, 1, 1, 1, 0};
        int _b = m->positions.size() / 3;
        for (int i = 0; i < 4; i++)
        {
            glm::vec3 v = l.edges[start_n + i];
            m->positions.push_back(v.x);
            m->positions.push_back(v.y);
            m->positions.push_back(v.z);
            m->normals.push_back(n.x);
            m->normals.push_back(n.y);
            m->normals.push_back(n.z);
            m->colors.push_back(tex_c[2 * i]);
            m->colors.push_back(tex_c[2 * i + 1]);
            m->colors.push_back(tc_zw.x);
            m->colors.push_back(tc_zw.y);
        }

        m->indices.push_back(_b);
        m->indices.push_back(_b + 1);
        m->indices.push_back(_b + 2);
        m->indices.push_back(_b + 2);
        m->indices.push_back(_b + 3);
        m->indices.push_back(_b);
    }
}
void Visualizer::get_ring(glm::vec3 &start, glm::vec3 &dir, float radius, SegmentVertexes &sv, int ring_size, 
                          float rel_ring_pos, std::vector<float> &mults, glm::vec3 p)
{
    sv.ringsize = ring_size;
    dir = glm::normalize(dir);
    if (abs(dot(dir,p)) > 0.9999)
    {
        p = glm::normalize(glm::vec3(1,0,0.001));
    }
    glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, p)), 1.0);

    glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)(2 * PI / ring_size), dir);

    if (ring_size == 0)
        logerr("0 ringsize !!!");
    for (int i = 0; i < ring_size; i++)
    {
        VertexData vd;
        vd.pos = start + radius*Branch::get_r_mult( i*2*PI/ring_size,mults) * glm::vec3(n.x, n.y, n.z);
        vd.normal = n;
        vd.tex_coord.x = ((float)i) / ring_size;
        vd.tex_coord.y = rel_ring_pos;
        sv.bigRing.push_back(vd);
        n = r * n;
    }
}
void Visualizer::get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale)
{
    sv.s = s;
    glm::vec3 start = s.begin;
    glm::vec3 end = s.end;
    glm::vec3 dir = end - start;
    get_ring(start, dir, scale*s.rel_r_begin, sv, ring_size, rel_ring_pos,s.mults);
}
void Visualizer::get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale)
{
    sv.s = s;
    glm::vec3 start = s.begin;
    glm::vec3 end = s.end;
    glm::vec3 dir = end - start;
    std::vector<VertexData> data = sv.bigRing;
    get_ring(end, dir, scale*s.rel_r_end, sv, ring_size, rel_ring_pos,s.mults);
    sv.smallRing = sv.bigRing;
    sv.bigRing = data;
}
void Visualizer::seg_vertexes_to_model(SegmentVertexes &sv, Model *m, glm::vec2 tc_zw)
{
    Model *h = m;
    int _b = h->positions.size() / 3;
    if (sv.smallRing.size() < sv.ringsize || sv.bigRing.size() < sv.ringsize)
        return;
    for (auto pos : sv.smallRing)
    {
        h->positions.push_back(pos.pos.x);
        h->positions.push_back(pos.pos.y);
        h->positions.push_back(pos.pos.z);
        h->normals.push_back(pos.normal.x);
        h->normals.push_back(pos.normal.y);
        h->normals.push_back(pos.normal.z);
        
        float col_x = 1 - abs(2.0f*pos.tex_coord.x - 1);
        h->colors.push_back(col_x);
        h->colors.push_back(pos.tex_coord.y);
        h->colors.push_back(tc_zw.x);
        h->colors.push_back(tc_zw.y);
    }
    int shift = 0;
    float best_match = glm::length(sv.bigRing.front().pos - sv.smallRing.front().pos);
    for (int i = 0;i<sv.bigRing.size();i++)
    {
        if (glm::length(sv.bigRing[i].pos - sv.smallRing.front().pos) < best_match)
        {
            best_match = glm::length(sv.bigRing[i].pos - sv.smallRing.front().pos);
            shift = i;
        }
    }
    shift = 0;
    for (int i = 0;i<sv.bigRing.size();i++)
    {
        int real_i = (i + shift) % sv.bigRing.size();
        VertexData pos = sv.bigRing[real_i];
        h->positions.push_back(pos.pos.x);
        h->positions.push_back(pos.pos.y);
        h->positions.push_back(pos.pos.z);
        
        h->normals.push_back(pos.normal.x);
        h->normals.push_back(pos.normal.y);
        h->normals.push_back(pos.normal.z);
        
        float col_x = 1 - abs(2.0f*i/sv.bigRing.size() - 1);

        h->colors.push_back(col_x);
        h->colors.push_back(pos.tex_coord.y);
        h->colors.push_back(tc_zw.x);
        h->colors.push_back(tc_zw.y);
    }
    int ringsize = sv.ringsize;
    for (int i = 0; i < ringsize; i++)
    {
        //Bottom Triangle
        h->indices.push_back(_b + i);
        h->indices.push_back(_b + (i + 1) % (ringsize));
        h->indices.push_back(_b + i + ringsize);
        //Upper Triangle
        h->indices.push_back(_b + i + ringsize);
        h->indices.push_back(_b + (i + 1) % (ringsize));
        h->indices.push_back(_b + (i + 1) % (ringsize) + ringsize);
    }
}
void Visualizer::joint_to_model(Joint &j, Model *m, bool leaves)
{
}
void Visualizer::recursive_branch_to_model(Branch &b, Model *m, bool leaves, float scale, int level_from, int level_to)
{
    if (b.level < 0 || b.level > level_to)
        return;

    if (b.level >= level_from)
    {
        if (!leaves)
        {
            std::vector<SegmentVertexes> vets;
            int i = 0;
            int ringsize = 3 * pow(2, MAX(3 - b.level,0));
            ringsize = 3;//this function is called only for billboards creation, so detailed branches are not needed
            float br = 0.5*(b.segments.front().rel_r_begin + b.segments.back().rel_r_end);
            for (auto &segment : b.segments)
            {
                SegmentVertexes vt;
                float dist = glm::length(segment.begin - b.segments.front().begin); 
                get_base_ring(segment, vt, ringsize, dist/(PI*br), scale);
                if (!vets.empty())
                    vets.back().smallRing = vt.bigRing;
                //segment_to_model(segment,m,leaves);
                vets.push_back(vt);
                i++;
            }
            if (!vets.empty())
            {
                float dist = glm::length(b.segments.back().end - b.segments.front().begin); 
                get_last_seg_vertexes(b.segments.back(), vets.back(), ringsize, dist/(PI*br), scale);
            }

            for (auto &vt : vets)
            {
                seg_vertexes_to_model(vt, m);
            }
        }
        else
        {
            for (auto &joint : b.joints)
            {
                if (joint.leaf)
                    leaf_to_model(*(joint.leaf), m, scale);
            }
        }
    }
    for (auto &joint : b.joints)
    {
        for (auto branch : joint.childBranches)
            recursive_branch_to_model(*branch, m, leaves, scale, level_from, level_to);
    }
}
int s_cnt(Branch *b)
{
    int cnt = b->segments.size();
    for (auto &joint : b->joints)
    {
        for (auto ch_b : joint.childBranches)
            cnt += s_cnt(ch_b);
    }
    return cnt;
}
void Visualizer::recursive_branch_to_model_fast(Branch &branch, Model *m, bool leaves, float scale, int level_from, 
                                                int level_to)
{
    /*
    constexpr int ringsize = 4;
    int cnt = s_cnt(&branch);
    m->indices.reserve(cnt*ringsize*6);
    m->positions.reserve((cnt+1)*ringsize*3);
    m->normals.reserve((cnt+1)*ringsize*3);
    m->colors.reserve((cnt+1)*ringsize*4);*/
    recursive_branch_to_model_fast_i(branch, m, leaves, scale, level_from, level_to);
}

void Visualizer::recursive_branch_to_model_fast_i(Branch &branch, Model *m, bool leaves, float scale, int level_from, 
                                                  int level_to)
{
    if (branch.level < 0 || branch.level > level_to)
        return;


    if (branch.level >= level_from)
    {
        if (!leaves)
        {   
            constexpr int ringsize = 4;
            //m->indices.reserve(m->indices.size() + branch.segments.size()*ringsize*6);
            //m->positions.reserve(m->positions.size() + (branch.segments.size()+1)*ringsize*3);
            //m->normals.reserve(m->normals.size() + (branch.segments.size()+1)*ringsize*3);
            //m->colors.reserve(m->colors.size() + (branch.segments.size()+1)*ringsize*4);
            float total_dist = 0;
            float base_r = 0.5f*(branch.segments.front().rel_r_begin + branch.segments.back().rel_r_end);
            for (auto &segment : branch.segments)
            {
                vec3 &pos = segment.begin;
                vec3 a = segment.end - segment.begin;
                float len = length(a);
                a = a / len;
                vec3 b = abs(a.y) > 1e-6 ? vec3(-a.y,a.x,0) : vec3(-a.z,0,a.x);
                b = normalize(b);
                vec3 c = cross(a,b);

                int _b = m->positions.size()/3;
                for (int i = 0; i < ringsize; i++)
                {
                    float angle = (2*PI*i)/ringsize;
                    vec3 n = b*cosf(angle) + c*sinf(angle);
                    vec3 p = pos + segment.rel_r_begin*n;
                    m->normals.push_back(n.x);
                    m->normals.push_back(n.y);
                    m->normals.push_back(n.z);

                    m->positions.push_back(p.x);
                    m->positions.push_back(p.y);
                    m->positions.push_back(p.z);

                    float tcy = 0.2*total_dist/base_r;
                    tcy = abs(tcy - (int)tcy - 0.5f);
                    m->colors.push_back(1 - abs(angle - PI)/PI);
                    m->colors.push_back(tcy);
                    m->colors.push_back(1);
                    m->colors.push_back(0);

                    //Bottom Triangle
                    m->indices.push_back(_b + i);
                    m->indices.push_back(_b + (i + 1) % (ringsize));
                    m->indices.push_back(_b + i + ringsize);
                    //Upper Triangle
                    m->indices.push_back(_b + i + ringsize);
                    m->indices.push_back(_b + (i + 1) % (ringsize));
                    m->indices.push_back(_b + (i + 1) % (ringsize) + ringsize);
                }
                total_dist += len;
            }
            {
            //last vertices
                auto &segment = branch.segments.back();
                            vec3 &pos = segment.end;
                vec3 a = segment.end - segment.begin;
                float len = length(a);
                a = a / len;
                vec3 b = abs(a.y) > 1e-6 ? vec3(-a.y,a.x,0) : vec3(-a.z,0,a.x);
                b = normalize(b);
                vec3 c = cross(a,b);

                for (int i = 0; i < ringsize; i++)
                {
                    float angle = 2*PI*i/ringsize;
                    vec3 n = b*cosf(angle) + c*sinf(angle);
                    vec3 p = pos + segment.rel_r_end*n;
                    m->normals.push_back(n.x);
                    m->normals.push_back(n.y);
                    m->normals.push_back(n.z);

                    m->positions.push_back(p.x);
                    m->positions.push_back(p.y);
                    m->positions.push_back(p.z);

                    float tcy = 0.2*total_dist/base_r;
                    tcy = abs(tcy - (int)tcy - 0.5f);
                    m->colors.push_back(1 - abs(angle - PI));
                    m->colors.push_back(tcy);
                    m->colors.push_back(1);
                    m->colors.push_back(0);
                }
            }
        }
        else
        {
            for (auto &joint : branch.joints)
            {
                if (joint.leaf)
                    leaf_to_model(*(joint.leaf), m, scale);
            }
        }
    }
    for (auto &joint : branch.joints)
    {
        for (auto ch_b : joint.childBranches)
            recursive_branch_to_model_fast_i(*ch_b, m, leaves, scale, level_from, level_to);
    }
}

void Visualizer::branch_to_model(Branch &b, Model *m, bool leaves)
{
    if (!leaves)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, MAX(3 - b.level,0));
        for (auto &segment : b.segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3, 1);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b.segments.back(), vets.back(), ringsize, (float)(i % 3) / 3, 1);
        }

        //seg_vertexes_to_model(vets.front(),m);
        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, m);
        }
        if (b.level == -1)
        {
            logerr("added %d/%d segments",vets.size(),b.segments.size());
            for (auto &vt : vets)
            {
                debug("(%f %f %f) %d (%f)",vt.smallRing.front().pos.x,vt.smallRing.front().pos.y,vt.smallRing.front().pos.z,
                ringsize,vt.s.end.x);
            }
            debugnl();
        }
    }
    else
    {
        for (auto &joint : b.joints)
            joint_to_model(joint, m, leaves);
    }
}
void Visualizer::segment_to_model(Segment &s, Model *m, bool leaves)
{

    glm::vec3 start = s.begin;
    glm::vec3 end = s.end;
    glm::vec3 dir = end - start;
    int ringsize = 8;
    Model *h = m;
    //Get Some Normal Vector
    glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
    glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

    //Add the Correct Number of Indices
    glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)3.141 / ringsize, dir);

    //Index Buffer
    int _b = h->positions.size() / 3;

    //GL TRIANGLES
    for (int i = 0; i < ringsize; i++)
    {
        //Bottom Triangle
        h->indices.push_back(_b + i * 2 + 0);
        h->indices.push_back(_b + (i * 2 + 2) % (2 * ringsize));
        h->indices.push_back(_b + i * 2 + 1);
        //Upper Triangle
        h->indices.push_back(_b + (i * 2 + 2) % (2 * ringsize));
        h->indices.push_back(_b + (i * 2 + 3) % (2 * ringsize));
        h->indices.push_back(_b + i * 2 + 1);
    }

    for (int i = 0; i < ringsize; i++)
    {

        h->positions.push_back(start.x + s.rel_r_begin * n.x);
        h->positions.push_back(start.y + s.rel_r_begin * n.y);
        h->positions.push_back(start.z + s.rel_r_begin * n.z);
        h->normals.push_back(n.x);
        h->normals.push_back(n.y);
        h->normals.push_back(n.z);
        n = r * n;

        h->positions.push_back(end.x + s.rel_r_end * n.x);
        h->positions.push_back(end.y + s.rel_r_end * n.y);
        h->positions.push_back(end.z + s.rel_r_end * n.z);
        h->normals.push_back(n.x);
        h->normals.push_back(n.y);
        h->normals.push_back(n.z);
        n = r * n;
    }
}
void Visualizer::add_branch_layer(Tree &t, int layer, Model *m)
{
    if (layer >= 0 && layer < t.branchHeaps.size())
    {
        for (auto &branch : t.branchHeaps[layer]->branches)
        {
            branch_to_model(branch, m, false);
        }
    }
}
void Visualizer::packed_branch_to_model(PackedBranch &b, Model *m, bool leaves, float precision, glm::vec2 tc_zw)
{
    if (!leaves)
    {
        if (b.joints.size() < 2)
            return;
        std::vector<SegmentVertexes> vets;
        std::vector<float> empty_mults;
        
        float diff = b.joints[0].r > 1 ? sqrt(b.joints[0].r) : - sqrt(1/b.joints[0].r);
        int ringsize = CLAMP(pow(2, precision + 1) + diff, 4, pow(2, precision + 2));
        int tex_step = 4;
        if (ringsize % 2 == 1)
            ringsize++;
        glm::vec3 dir = b.joints[1].pos - b.joints[0].pos;
        float br = 0.5*(b.joints.front().r + b.joints.back().r);
        for (int i = 0; i < b.joints.size(); i++)
        {
            if (i != 0)
                dir = b.joints[i].pos - b.joints[i - 1].pos;
            SegmentVertexes vt;
            std::vector<float> &mults = b.r_mults.size() > i ? b.r_mults[i] : empty_mults;
            float dist = glm::length(b.joints[i].pos - b.joints.front().pos); 
            get_ring(b.joints[i].pos, dir, b.joints[i].r, vt, ringsize, dist/(2*PI*br), mults, b.plane_coef);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            vets.push_back(vt);
        }
        SegmentVertexes vt;
        int i = b.joints.size() - 1;
        std::vector<float> &mults = b.r_mults.size() > i ? b.r_mults[i] : empty_mults;
        float dist = glm::length(b.joints[i].pos - b.joints.front().pos); 
        get_ring(b.joints[i].pos, dir, b.joints[i].r, vt, ringsize, dist/(2*PI*br), mults, b.plane_coef);
        vets.back().smallRing = vt.bigRing;

        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, m, tc_zw);
        }
    }
    else
    {
        for (auto &l : b.leaves)
        {
            packed_leaf_to_model(l,m, tc_zw);
        }
    }
}
Visualizer::Visualizer():
tree_tex(textureManager.empty()),
leaves_tex(textureManager.empty())
{

}
Visualizer::Visualizer(Texture _tree_tex, Texture _leaves_tex, Shader *_tree_shader):
tree_tex(_tree_tex),
leaves_tex(_leaves_tex)

{
    tree_shader = _tree_shader;
}
void Visualizer::box_to_model(Box *box, Model *m)
{
    std::vector<glm::vec3> pos;
    std::vector<glm::vec2> colors;
    for (int i = 0; i <= 1; i++)
    {
        for (int j = 0; j <= 1; j++)
        {
            for (int k = 0; k <= 1; k++)
            {
                pos.push_back(box->pos + (float)i * box->a + (float)j * box->b + (float)k * box->c);
                colors.push_back(glm::vec2(k, j));
            }
        }
    }
    int indicies[36] = {0, 2, 1, 3, 1, 2, 2, 3, 6, 7, 6, 3, 7, 5, 6, 4, 6, 5, 0, 1, 4, 5, 1, 4, 0, 4, 2, 6, 2, 4, 1, 5, 3, 7, 3, 5};
    int _b = m->positions.size() / 3;
    for (int i = 0; i < pos.size(); i++)
    {
        m->positions.push_back(pos[i].x);
        m->positions.push_back(pos[i].y);
        m->positions.push_back(pos[i].z);

        m->colors.push_back(colors[i].x);
        m->colors.push_back(colors[i].y);
        m->colors.push_back(0.0);
        m->colors.push_back(1.0);
    }
    for (int i = 0; i < 36; i++)
    {
        m->indices.push_back(_b + indicies[i]);
    }
}
void Visualizer::ellipsoid_to_model(Ellipsoid *b, Model *m, int sectors, int stacks, bool smooth)
{
    float x, y, z, xy;                              // vertex position
    float nx, ny, nz, lengthInv = 1.0f;             // normal
    float s, t;                                     // texCoord

    float sectorStep = 2 * PI / sectors;
    float stackStep = PI / stacks;
    float sectorAngle, stackAngle;
    int _b = m->positions.size() / 3;

    for(int i = 0; i <= stacks; ++i)
    {
        stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = cosf(stackAngle);             // r * cos(u)
        z = sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for(int j = 0; j <= sectors; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            glm::vec3 pos = glm::inverse(b->transform)*glm::vec4(x,y,z,1); 
            m->positions.push_back(pos.x);
            m->positions.push_back(pos.y);
            m->positions.push_back(pos.z);

            // normalized vertex normal
            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;
            m->normals.push_back(nx);
            m->normals.push_back(ny);
            m->normals.push_back(nz);

            // vertex tex coord between [0, 1]
            s = (float)j / sectors;
            t = (float)i / stacks;
            m->colors.push_back(s);
            m->colors.push_back(t);
            m->colors.push_back(0.0);
            m->colors.push_back(1.0);
        }
    }

    // indices
    //  k1--k1+1
    //  |  / |
    //  | /  |
    //  k2--k2+1
    unsigned int k1, k2;

    for(int i = 0; i < stacks; ++i)
    {
        k1 = i * (sectors + 1);     // beginning of current stack
        k2 = k1 + sectors + 1;      // beginning of next stack

        for(int j = 0; j < sectors; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding 1st and last stacks
            if(i != 0)
            {
                m->indices.push_back(_b + k1);
                m->indices.push_back(_b + k2);
                m->indices.push_back(_b + k1 + 1);
            }

            if(i != (stacks-1))
            {
                m->indices.push_back(_b + k1 + 1);
                m->indices.push_back(_b + k2);
                m->indices.push_back(_b + k2 + 1);
            }
        }
    }
}
void Visualizer::cylinder_to_model(Cylinder *b, Model *m, int sectors)
{
    int v0 = m->positions.size()/3;
    for (int i=0;i<sectors;i++)
    {
        glm::vec3 vert = b->pos - b->c + cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b;
        m->positions.push_back(vert.x);
        m->positions.push_back(vert.y);
        m->positions.push_back(vert.z);

        glm::vec3 n = glm::normalize(cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b);
        m->normals.push_back(n.x);
        m->normals.push_back(n.y);
        m->normals.push_back(n.z);

        float col_x = 1 - abs(2.0f*i/sectors - 1);
        m->colors.push_back(col_x);
        m->colors.push_back(0.0);
        m->colors.push_back(0.0);
        m->colors.push_back(1.0);
    }
    int v1 = m->positions.size()/3;
    for (int i=0;i<sectors;i++)
    {
        glm::vec3 vert = b->pos + b->c + cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b;
        m->positions.push_back(vert.x);
        m->positions.push_back(vert.y);
        m->positions.push_back(vert.z);

        glm::vec3 n = glm::normalize(cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b);
        m->normals.push_back(n.x);
        m->normals.push_back(n.y);
        m->normals.push_back(n.z);

        float col_x = 1 - abs(2.0f*i/sectors - 1);
        m->colors.push_back(col_x);
        m->colors.push_back(1.0);
        m->colors.push_back(0.0);
        m->colors.push_back(1.0);
    }
    int v_down = m->positions.size()/3;

    glm::vec3 vert = b->pos - b->c;
    m->positions.push_back(vert.x);
    m->positions.push_back(vert.y);
    m->positions.push_back(vert.z);

    glm::vec3 n = glm::normalize(-b->c);
    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);
    
    m->colors.push_back(0.0);
    m->colors.push_back(1.0);
    m->colors.push_back(0.0);
    m->colors.push_back(1.0);

    int v_up = m->positions.size()/3;

    vert = b->pos + b->c;
    m->positions.push_back(vert.x);
    m->positions.push_back(vert.y);
    m->positions.push_back(vert.z);

    n = glm::normalize(b->c);
    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);
    
    m->colors.push_back(0.0);
    m->colors.push_back(0.0);
    m->colors.push_back(0.0);
    m->colors.push_back(1.0);

    for (int i=0;i<sectors;i++)
    {
        int i1 = (i+1)%sectors;
        m->indices.push_back(v0 + i);
        m->indices.push_back(v1 + i1);
        m->indices.push_back(v0 + i1);

        m->indices.push_back(v1 + i);
        m->indices.push_back(v1 + i1);
        m->indices.push_back(v0 + i);

        m->indices.push_back(v_down);
        m->indices.push_back(v0 + i1);
        m->indices.push_back(v0 + i);

        m->indices.push_back(v_up);
        m->indices.push_back(v1 + i);
        m->indices.push_back(v1 + i1);
    }
}
void Visualizer::body_to_model(Body *b, Model *m, bool fixed_tc, glm::vec4 tc)
{
    int last_tc = m->colors.size();
    Box *box = dynamic_cast<Box *>(b);
    if (box)
        box_to_model(box,m);
    Ellipsoid *el = dynamic_cast<Ellipsoid *>(b);
    if (el)
        ellipsoid_to_model(el,m,20,20);
    Cylinder *cyl = dynamic_cast<Cylinder *>(b);
    if (cyl)
        cylinder_to_model(cyl,m,20);
    if (fixed_tc)
    {
        for (int i=last_tc;i<m->colors.size() - 3;i+=4)
        {
            m->colors[i] = tc.x;
            m->colors[i+1] = tc.y;
            m->colors[i+2] = tc.z;
            m->colors[i+3] = tc.w;
        }
    }
}

int HMLOD = 0;
glm::vec3 get_n(Heightmap &h, glm::vec3 terr_pos, glm::vec2 step)
{
    terr_pos.y = 0;
    glm::vec3 terr_pos1 = terr_pos + glm::vec3(step.x,0,0);
    glm::vec3 terr_pos2 = terr_pos + glm::vec3(0,0,step.y);
    terr_pos.y = h.get_height(terr_pos);
    terr_pos1.y = h.get_height(terr_pos1);
    terr_pos2.y = h.get_height(terr_pos2);
    glm::vec3 n = glm::normalize(glm::cross(terr_pos1 - terr_pos, terr_pos2 - terr_pos));
    if (dot(n, glm::vec3(0,1,0)) < 0)
        n = -n;
    return n;
}

void add_vertex(Heightmap &h, glm::vec3 &terr_pos, Model *m, glm::vec2 step)
{
    terr_pos.y = h.get_height(terr_pos);
    glm::vec3 n = get_n(h, terr_pos, step);

    m->positions.push_back(terr_pos.x);
    m->positions.push_back(terr_pos.y);
    m->positions.push_back(terr_pos.z);

    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);

    float l = 0.1 * length(h.get_grad(terr_pos));
    glm::vec2 tc = glm::vec2(2.0f*abs(fract(glm::vec2(terr_pos.x, terr_pos.z)/47.1f) - 0.5f));
    tc = glm::vec2(terr_pos.x + urand(-HMLOD,HMLOD), terr_pos.z + urand(-HMLOD,HMLOD))/27.1f;
    m->colors.push_back(tc.x);
    m->colors.push_back(tc.y);
    m->colors.push_back(0);
    m->colors.push_back(0);
}

void Visualizer::heightmap_to_model(Heightmap &h, Model *m, glm::vec2 detailed_size, glm::vec2 full_size, 
                                    float precision, int LODs)
{
    int x = MAX(pow(2, round(log2(2*detailed_size.x / precision))), 4);
    int y = MAX(pow(2, round(log2(2*detailed_size.y / precision))), 4);
    x = h.get_grid_size().x;
    y = h.get_grid_size().y;
    //x = 4;
    //y = 4;
    glm::vec2 step = 2.0f * detailed_size / glm::vec2(x, y);
    glm::vec3 center_pos = h.get_pos();
    glm::vec2 total_size = detailed_size;
    HMLOD = 0;
    //detailed part of hmap
    for (int i = 0; i <= x; i++)
    {
        for (int j = 0; j <= y; j++)
        {
            int ind = m->positions.size() / 3;
            glm::vec3 terr_pos = glm::vec3(center_pos.x - detailed_size.x + step.x * i, 0, 
                                           center_pos.z - detailed_size.y + step.y * j);
            add_vertex(h, terr_pos, m, step);

            if (i != x && j != y)
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 1);
                m->indices.push_back(ind + (y+1) + 1);
                m->indices.push_back(ind);
                m->indices.push_back(ind + (y+1) + 1);
                m->indices.push_back(ind + (y+1));
            }
        }
    }
    bool sides_spec = y % 2;
    x = x/2 + 2;//we can safly do it, as x and y are powers of 2
    y = y/2;
    step *= 2.0f;
    HMLOD = 1;
    vec3 pos;
    glm::vec3 start_pos_top = glm::vec3(center_pos.x - detailed_size.x - step.x, 0, center_pos.z + detailed_size.y + step.y);
    glm::vec3 start_pos_right = glm::vec3(center_pos.x + detailed_size.x + step.x, 0, center_pos.z - detailed_size.y);
    glm::vec3 start_pos_bottom = glm::vec3(center_pos.x - detailed_size.x - step.x, 0, center_pos.z - detailed_size.y - step.y);
    glm::vec3 start_pos_left = glm::vec3(center_pos.x - detailed_size.x - step.x, 0, center_pos.z - detailed_size.y);
    glm::vec3 end_pos_bottom = glm::vec3(center_pos.x + detailed_size.x + step.x, 0, center_pos.z - detailed_size.y);
    glm::vec3 end_pos_top = glm::vec3(center_pos.x + detailed_size.x + step.x, 0, center_pos.z + detailed_size.y);

    while (total_size.x < full_size.x || total_size.y < full_size.y)
    {
       
        vec3 start_pos = start_pos_top;
        for (int i = 0;i<x;i++)
        {
            bool less = (i == 0) || (i == x-1);
            int ind = m->positions.size() / 3;

            pos = start_pos + vec3(i*step.x,0,0);
            add_vertex(h, pos, m, step);
            pos = start_pos + vec3(i*step.x,0,-step.y);
            add_vertex(h, pos, m, step);
            if (!less)
            {
              pos = start_pos + vec3((i+0.5)*step.x,0,-step.y);
              add_vertex(h, pos, m, step);
            }
            if (i == x - 1)
                pos = end_pos_top;
            else
                pos = start_pos + vec3((i+1)*step.x,0,-step.y);
            add_vertex(h, pos, m, step);
            if (i == x - 1)
                pos = end_pos_top + vec3(0,0,step.y);
            else
                pos = start_pos + vec3((i+1)*step.x,0,0);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);

            if (less)
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
            else
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 4);
                
                m->indices.push_back(ind + 4);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
        }

        start_pos = start_pos_bottom;
        for (int i = 0;i<x;i++)
        {
            bool less = (i == 0) || (i == x-1);
            int ind = m->positions.size() / 3;

            pos = start_pos + vec3(i*step.x,0,0);
            add_vertex(h, pos, m, step);
            pos = start_pos + vec3(i*step.x,0,step.y);
            add_vertex(h, pos, m, step);
            if (!less)
            {
                pos = start_pos + vec3((i+0.5)*step.x,0,step.y);
                add_vertex(h, pos, m, step);
            }
            
            if (i == x - 1)
                pos = end_pos_bottom;
            else
                pos = start_pos + vec3((i+1)*step.x,0,step.y);
            add_vertex(h, pos, m, step);
            if (i == x - 1)
                pos = end_pos_bottom + vec3(0,0,-step.y);
            else
                pos = start_pos + vec3((i+1)*step.x,0,0);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            if (less)
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
            else
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 4);
                
                m->indices.push_back(ind + 4);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
        }

    
        
        start_pos = start_pos_left;
        for (int i=0;i<y - sides_spec;i++)
        {
            int ind = m->positions.size() / 3;

            pos = start_pos + vec3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + vec3(step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + vec3(step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + vec3(step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + vec3(0,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);
        }
        if (sides_spec)
        {
            int i = y -1;
            int ind = m->positions.size() / 3;

            pos = start_pos + vec3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + vec3(step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + vec3(step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + vec3(step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos_top + vec3(0,0,-step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos_top + vec3(step.x,0,-step.y);
            add_vertex(h, pos, m, step);

            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);

            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 3);
            m->indices.push_back(ind + 5);
        }

          start_pos = start_pos_right;
        for (int i=0;i<y - sides_spec;i++)
        {
            int ind = m->positions.size() / 3;

            pos = start_pos + vec3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + vec3(-step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + vec3(-step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + vec3(-step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + vec3(0,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);
        }
        if (sides_spec)
        {
            int i = y -1;
            int ind = m->positions.size() / 3;

            pos = start_pos + vec3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + vec3(-step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + vec3(-step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + vec3(-step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = end_pos_top;
            add_vertex(h, pos, m, step);
            
            pos = end_pos_top + vec3(-step.x,0,0);
            add_vertex(h, pos, m, step);

            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);

            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 3);
            m->indices.push_back(ind + 5);
        }

        //break;
        start_pos_top += glm::vec3(-2.0f*step.x, 0, 2.0f*step.y);
        start_pos_bottom += glm::vec3(-2.0f*step.x, 0, -2.0f*step.y);
        end_pos_top += glm::vec3(2.0f*step.x, 0, step.y);
        start_pos_right += glm::vec3(2.0f*step.x, 0, -step.y);
        start_pos_left += glm::vec3(-2.0f*step.x, 0, -step.y);
        end_pos_bottom += glm::vec3(2.0f*step.x, 0, -step.y);

        sides_spec = sides_spec || y % 2;

        x = x/2 + 2;
        y = y/2 + 1;
        step *= 2.0f;
        HMLOD*=2;
        total_size += step;
    }

    m->update();
}