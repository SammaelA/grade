#include "tree_modeling.h"
#include "common_utils/utility.h"
#include "tinyEngine/engine.h"

#define PI 3.14159265f
using namespace glm;
namespace visualizer
{
  void leaf_to_model(Leaf &l, Mesh *m, float scale)
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
void packed_leaf_to_model(PackedLeaf &l, Mesh *m, glm::vec2 tc_zw, bool need_tangent)
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
        glm::vec3 tangent = glm::normalize(a-b);
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
            if (need_tangent)
            {
                m->tangents.push_back(tangent.x);
                m->tangents.push_back(tangent.y);
                m->tangents.push_back(tangent.z);
            }
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
void get_ring(glm::vec3 &start, glm::vec3 &dir, float radius, SegmentVertexes &sv, int ring_size, 
                          float rel_ring_pos, std::vector<float> &mults, glm::vec3 p)
{
    sv.ringsize = ring_size;
    dir = glm::normalize(dir);
    if (abs(dot(dir,p)) > 0.9999)
    {
        p = glm::normalize(glm::vec3(1,0,0.001));
    }
    glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, p)), 1.0);

    glm::mat4 r = LiteMath::rotate(glm::mat4(1.0), (float)(2 * PI / ring_size), dir);

    if (ring_size == 0)
        logerr("0 ringsize !!!");
    for (int i = 0; i < ring_size; i++)
    {
        VertexData vd;
        vd.pos = start + radius*Branch::get_r_mult( i*2*PI/ring_size,mults) * glm::vec3(n.x, n.y, n.z);
        vd.normal = n;
        vd.tangent = glm::cross(dir, vd.normal);
        vd.tex_coord.x = ((float)i) / ring_size;
        vd.tex_coord.y = rel_ring_pos;
        sv.bigRing.push_back(vd);
        n = r * n;
    }
}
void get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale)
{
    sv.s = s;
    glm::vec3 start = s.begin;
    glm::vec3 end = s.end;
    glm::vec3 dir = end - start;
    get_ring(start, dir, scale*s.rel_r_begin, sv, ring_size, rel_ring_pos,s.mults);
}
void get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos, float scale)
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
void seg_vertexes_to_model(SegmentVertexes &sv, Mesh *m, glm::vec2 tc_zw, bool need_tangents)
{
    Mesh *h = m;
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
        if (need_tangents)
        {
            h->tangents.push_back(pos.tangent.x);
            h->tangents.push_back(pos.tangent.y);
            h->tangents.push_back(pos.tangent.z);       
        }
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
        if (need_tangents)
        {
            h->tangents.push_back(pos.tangent.x);
            h->tangents.push_back(pos.tangent.y);
            h->tangents.push_back(pos.tangent.z);            
        }
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
void joint_to_model(Joint &j, Mesh *m, bool leaves)
{
}
void recursive_branch_to_model(Branch &b, Mesh *m, bool leaves, float scale, int level_from, int level_to)
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
void recursive_branch_to_model_fast(Branch &branch, Mesh *m, bool leaves, float scale, int level_from, 
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

void recursive_branch_to_model_fast_i(Branch &branch, Mesh *m, bool leaves, float scale, int level_from, 
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

void branch_to_model(Branch &b, Mesh *m, bool leaves)
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
void segment_to_model(Segment &s, Mesh *m, bool leaves)
{

    glm::vec3 start = s.begin;
    glm::vec3 end = s.end;
    glm::vec3 dir = end - start;
    int ringsize = 8;
    Mesh *h = m;
    //Get Some Normal Vector
    glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
    glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

    //Add the Correct Number of Indices
    glm::mat4 r = LiteMath::rotate(glm::mat4(1.0), (float)3.141 / ringsize, dir);

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
void add_branch_layer(Tree &t, int layer, Mesh *m)
{
    if (layer >= 0 && layer < t.branchHeaps.size())
    {
        for (auto &branch : t.branchHeaps[layer]->branches)
        {
            branch_to_model(branch, m, false);
        }
    }
}
void packed_branch_to_model(PackedBranch &b, Mesh *m, bool leaves, float precision, glm::vec2 tc_zw, bool need_tangents)
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
            seg_vertexes_to_model(vt, m, tc_zw, need_tangents);
        }
    }
    else
    {
        for (auto &l : b.leaves)
        {
            packed_leaf_to_model(l,m, tc_zw, need_tangents);
        }
    }
}
}