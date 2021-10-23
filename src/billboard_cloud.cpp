#include "billboard_cloud.h"
#include <algorithm>
#include <vector>
#include <glm/gtc/matrix_transform.hpp>
#include "tree.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "distribution.h"

using namespace glm;
#define TEX_ATLAS_LAYERS 4
BillboardCloudRaw::BillboardCloudRaw() :
Countable(2),
rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"})                            
{   
}

BillboardCloudRaw::~BillboardCloudRaw()
{
    if (atlas)
        delete atlas;
}
BillboardCloudRaw::AtlasParams BillboardCloudRaw::set_atlas_params(int quality, int cnt)
{
    const int fallback_size = 1024, fallback_layers = 4;
    int groups = cnt;
    GLint max_tex_size, max_layers;

    glGetIntegerv(GL_MAX_TEXTURE_SIZE,&max_tex_size);
    glGetIntegerv(GL_MAX_ARRAY_TEXTURE_LAYERS, &max_layers);

    int cnt_x = MIN(max_tex_size/quality,ceil(sqrt(groups)));
    int cnt_y = MIN(max_tex_size/quality,ceil(sqrt(groups)));
    if (cnt_x*cnt_y == 0)
    {
        logerr("Unable to create texture array with billboard size = %d. It is more that max texture size",quality);
        AtlasParams par;
        par.x = fallback_size;
        par.y = fallback_size;
        par.layers = fallback_layers;
        par.valid = false;
        return par;
    }
    int tex_x = cnt_x * quality;
    int tex_y = cnt_y * quality;
    int layers = ceil((float)groups/(cnt_x*cnt_y));
    if (layers > max_layers)
    {
        logerr("Unable to create texture atlas for %d billboards. It requires too many layers",groups);
        AtlasParams par;
        par.x = fallback_size;
        par.y = fallback_size;
        par.layers = fallback_layers;
        par.valid = false;
        return par;
    }

    AtlasParams par;
    par.x = tex_x;
    par.y = tex_y;
    par.layers = layers;
    par.valid = true;
    par.grid_x = quality;
    par.grid_y = quality;
    debugl(10,"successfully created atlas %dx%dx%d for %d (max = %d) billboards\n",
          par.x,par.y,par.layers,groups,cnt_x*cnt_y*layers);
    
    return par;
}

bool BillboardCloudRaw::BPD_comp(BranchProjectionData &a, BranchProjectionData &b)
{
    return a.projection_err > b.projection_err;
}
float BillboardCloudRaw::projection_error_rec(Branch *b, vec3 &n, float d)
{
    if (!b || b->joints.size() == 0)
        return 0;
    float err = 0.0;
    for (auto &j : b->joints)
    {
        err += abs(dot(j.pos, n) + d);
        for (auto br : j.childBranches)
            err += projection_error_rec(br, n, d);
    }
    return err;
}

void BillboardCloudRaw::create_billboard(TreeTypeData &ttd, Branch *branch, BBox &min_bbox, Visualizer &tg, int num, Billboard &bill,
                                         TextureAtlas &atlas, BillboardGenerationParams params)
{
    if (num < 0)
    {
        logerr("too many billboards = %d", billboard_count);
        return;
    }
    mat4 transl = translate(mat4(1.0f), -1.0f * min_bbox.position);
    mat4 SC = scale(mat4(1.0f), min_bbox.sizes);
    mat4 SC_inv = inverse(SC);
    mat4 rot_inv(vec4(min_bbox.a, 0), vec4(min_bbox.b, 0), vec4(min_bbox.c, 0), vec4(0, 0, 0, 1));
    mat4 rot = inverse(rot_inv);
    mat4 ort = ortho(-1, 1, -1, 1, 1, -1);
    vec4 rb = transl * rot * vec4(branch->joints.front().pos,1);
    mat4 tex_sh = scale(mat4(1), vec3(2, 2, 2));
    mat4 tex_tr = translate(mat4(1), vec3(-1, -1, -1));
    mat4 atlas_tr = atlas.tex_transform(num);
    mat4 result = ort * tex_tr * tex_sh * atlas_tr * SC_inv;
    Model bm;
    std::function<void(Model *)> _c_wood = [&](Model *h) { 
        tg.recursive_branch_to_model(*branch, &bm, false, params.wood_scale, params.level_from, params.level_to); 
    };
    std::function<void(Model *)> _c_leaves = [&](Model *h) { 
        tg.recursive_branch_to_model(*branch, &bm, true, params.leaf_scale, params.level_from, params.level_to); 
    };

    int tex_count = params.normals_needed ? 1 : atlas.tex_count();
    for (int k = 0; k<tex_count;k++)
    {
        bool transparent_pass = (k == 0 && params.leaf_opacity < 1);

        atlas.target(num,k);
        rendererToTexture.use();

        bm.construct(_c_wood);
        rendererToTexture.texture("tex", ttd.wood);
        rendererToTexture.uniform("model", transl * rot * bm.model);
        rendererToTexture.uniform("projectionCamera", result);
        rendererToTexture.uniform("state", params.monochrome ? -1 : k);
        rendererToTexture.uniform("projection_zero", rb.z);
        if (params.monochrome)
            rendererToTexture.uniform("fixed_color",glm::vec4(1,0,0,1));
        bm.render(GL_TRIANGLES);
        
        if (params.leaf_opacity > 0)
        {
            if (transparent_pass)
            {
                glEnable(GL_BLEND); 
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
                glDisable(GL_DEPTH_TEST); 
            }
            bm.construct(_c_leaves);

            Texture &leaf = ttd.leaf;
            glBindTexture(leaf.type,leaf.texture);
            glTexParameteri(leaf.type, GL_TEXTURE_BASE_LEVEL, 0);
            glTexParameteri(leaf.type, GL_TEXTURE_MAX_LEVEL, 0);

            rendererToTexture.texture("tex", ttd.leaf);
            rendererToTexture.uniform("model", transl * rot * bm.model);
            rendererToTexture.uniform("projectionCamera", result);
            rendererToTexture.uniform("projection_zero", rb.z);
            if (params.monochrome)
                rendererToTexture.uniform("fixed_color",glm::vec4(0,1,0,params.leaf_opacity));
            bm.render(GL_TRIANGLES);

            glTexParameteri(leaf.type, GL_TEXTURE_BASE_LEVEL, 0);
            glTexParameteri(leaf.type, GL_TEXTURE_MAX_LEVEL, 1000);
            if (transparent_pass)
            {
                glDisable(GL_BLEND);   
                glEnable(GL_DEPTH_TEST); 
            }
        }

        billboards.push_back(bill);
    }
}

void BillboardCloudRaw::create_billboard(TreeTypeData &ttd, Branch *branch, BBox &min_bbox, Visualizer &tg, int num,
                                         Billboard &bill, float leaf_scale,  float wood_scale, 
                                         bool monochrome, int level_from, int level_to)
{
    BillboardGenerationParams bgp;
    bgp.leaf_scale = leaf_scale;
    bgp.wood_scale = wood_scale;
    bgp.monochrome = monochrome;
    bgp.level_from = level_from;
    bgp.level_to = level_to;

    create_billboard(ttd, branch, min_bbox, tg, num, bill, *atlas, bgp);
}
void center_of_mass(Branch *b,glm::vec3 &sum, double &mass)
{
    for (Segment &s : b->segments)
    {
        float v = glm::length(s.begin - s.end)*(s.rel_r_begin*s.rel_r_begin + s.rel_r_begin*s.rel_r_end + 
                                                s.rel_r_end*s.rel_r_end);
        sum += 0.5f*v*(s.begin + s.end);
        mass += v;
    }
    for (Joint &j : b->joints)
    {
        for (Branch *br : j.childBranches)
        {
            center_of_mass(br, sum, mass);
        }
    }
}
BBox BillboardCloudRaw::get_minimal_bbox_fixed_dirs(Branch *branch)
{
    vec3 a(0, 0, 0);
    vec3 b;
    vec3 c;
    a = (branch->segments.back().end - branch->segments.front().begin);
    a = normalize(a); //average branch direction
    double sum = 1e-6;
    vec3 v(0,0,0);
    center_of_mass(branch, v, sum);
    v = ((float)(1.0/sum))*v - branch->segments.front().begin;

    b = v - a*dot(a,v);
    if (dot(b,b) < 1e-6)
    {
        if (abs(1 - dot(a,vec3(1,0,0))) > 1e-6)
            b = vec3(1,0,0);
        else
            b = vec3(0,0,1); 
    }
    b = normalize(b);
    c = cross(a, b);
    return BillboardCloudRaw::get_bbox(branch, a, b, c);
}
BBox BillboardCloudRaw::get_minimal_bbox(Branch *branch)
{
    return get_minimal_bbox_fixed_dirs(branch);
    int iterations = 360;
    vec3 a(0, 0, 0);
    vec3 b;
    vec3 c;
    a = (branch->segments.back().end - branch->segments.front().begin);
    a = normalize(a); //average branch direction

    b.x = urand();
    b.y = urand();
    b.z = urand();
    b = normalize(b - a * dot(a, b));
    c = cross(a, b);
    mat4 br = rotate(mat4(1.0), (float)(2 * PI / iterations), a);
    BBox min_bbox;
    float min_minside = 1e10;
    for (int i = 0; i < iterations; i++)
    {
        BBox box = get_bbox(branch, a, b, c);
        float minside = MIN(MIN(box.sizes.x, box.sizes.y), box.sizes.z);
        if (minside < min_minside)
        {
            min_bbox = box;
            min_minside = minside;
        }
        b = br * vec4(b, 1);
        c = cross(a, b);
    }
    return min_bbox;
}
BBox BillboardCloudRaw::get_bbox(Branch *branch, glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
    vec4 bias = vec4(1, 1, 1, 0);
    mat4 rot_inv(vec4(a, 0), vec4(b, 0), vec4(c, 0), vec4(0, 0, 0, 1));
    mat4 rot = inverse(rot_inv);
    //transform from model to bbox coordinates
    vec4 mx(-1e10, -1e10, -1e10, 1);
    vec4 mn(1e10, 1e10, 1e10, 1);
    BBox box;
    update_bbox(branch, rot, mn, mx);
    mn -= bias;
    mx += bias;
    box.sizes = mx - mn;
    box.position = mn;
    box.a = a;
    box.b = b;
    box.c = c;
    return box;
}
void BillboardCloudRaw::update_bbox(Branch *branch, mat4 &rot, vec4 &mn, vec4 &mx)
{
    for (Joint &j : branch->joints)
    {
        vec4 pos = rot * vec4(j.pos, 1);
        mn = min(mn, pos);
        mx = max(mx, pos);
        if (j.leaf)
        {
            for (auto &vert : j.leaf->edges)
            {
                pos = rot * vec4(vert, 1);
                mn = min(mn, pos);
                mx = max(mx, pos);
            }
        }
        for (Branch *br : j.childBranches)
            update_bbox(br, rot, mn, mx);
    }
}

std::vector<glm::vec3> Billboard::get_tc(TextureAtlas &atlas)
{
    if (positions.size() == 4)
    {
        std::vector<glm::vec3> tcs;
        std::vector<float> tex_c{0,0, 1,0, 0,1, 1,1};
        for (int i = 0; i < 4; i++)
        {
            tcs.push_back(vec3(tex_c[2 * i], tex_c[2 * i + 1],0));
            atlas.process_tc(id, tcs.back());
        }
        return tcs;
    }
}
void Billboard::to_model(Model *m, TextureAtlas &atlas)
{
    if (positions.size() == 4)
    {
        int _b = m->positions.size() / 3;
        glm::vec3 a = positions[0];
        glm::vec3 b = positions[1];
        glm::vec3 c = positions[2];
        glm::vec3 n = glm::normalize(glm::cross(a - b, c - b));
        std::vector<float> tex_c{0, 0, 1, 0, 1, 1, 0, 1};
        for (int i = 0; i < 4; i++)
        {
            glm::vec3 v = positions[i];
            glm::vec3 tc = vec3(tex_c[2 * i], tex_c[2 * i + 1],0);
            atlas.process_tc(id, tc);
            m->positions.push_back(v.x);
            m->positions.push_back(v.y);
            m->positions.push_back(v.z);
            m->normals.push_back(n.x);
            m->normals.push_back(n.y);
            m->normals.push_back(n.z);
            m->colors.push_back(tc.x);
            m->colors.push_back(tc.y);
            m->colors.push_back(tc.z);
            m->colors.push_back(1);
        }

        m->indices.push_back(_b);
        m->indices.push_back(_b + 1);
        m->indices.push_back(_b + 3);
        m->indices.push_back(_b + 1);
        m->indices.push_back(_b + 2);
        m->indices.push_back(_b + 3);
    }
}
Billboard::Billboard(const BBox &box, int id, int branch_id, int type, glm::vec3 base_joint, bool _instancing)
{
    this->id = id;
    this->branch_id = branch_id;
    this->instancing = _instancing;
    mat4 rot_inv(vec4(box.a, 0), vec4(box.b, 0), vec4(box.c, 0), vec4(0, 0, 0, 1));
    mat4 rot = inverse(rot_inv);
    vec3 base_joint_rel = rot * vec4(base_joint, 1.0f);
    base_joint_rel -= box.position;
    vec4 pos = rot_inv * vec4(box.position, 1.0f);
    if (type == 1 || type == 0)
    {
        vec3 npos = vec3(pos.x, pos.y, pos.z);
        if (type == 1)
            npos += base_joint_rel.z * box.c;
        positions.push_back(npos);
        positions.push_back(npos + box.sizes.x * box.a);
        positions.push_back(npos + box.sizes.x * box.a + box.sizes.y * box.b);
        positions.push_back(npos + box.sizes.y * box.b);

        float d = -dot(box.c, npos);
        planeCoef = vec4(box.c.x, box.c.y, box.c.z, d);
    }
}
void BillboardCloudRaw::split_IDA_by_type(InstanceDataArrays &IDA, std::vector<InstanceDataArrays> &res)
{
    std::map<int, int> arr_num_by_id;
    for (int i=0;i<IDA.centers_self.size();i++)
    {
        auto it = arr_num_by_id.find(IDA.type_ids[i]);
        if (it == arr_num_by_id.end())
        {
            arr_num_by_id.emplace(IDA.type_ids[i], res.size());
            res.push_back(InstanceDataArrays());
            it = arr_num_by_id.find(IDA.type_ids[i]);
        }

        res[it->second].type_ids.push_back(IDA.type_ids[i]);
        res[it->second].tree_ids.push_back(IDA.tree_ids[i]);
        res[it->second].centers_par.push_back(IDA.centers_par[i]);
        res[it->second].centers_self.push_back(IDA.centers_self[i]);
        res[it->second].transforms.push_back(IDA.transforms[i]);
    }
    return;
    int sz = IDA.type_ids.size();
    if (IDA.centers_par.size() != sz || IDA.centers_self.size() != sz || IDA.transforms.size() != sz)
    {
        logerr("trying to split corrupted IDA");
        return;
    }
    std::vector<bool> need_change(sz,true);
    bool work_left = (sz != 0);
    while (work_left)
    {
        work_left = false;
        int cur_type = -1;
        res.push_back(InstanceDataArrays());
        for (int i=0;i<sz;i++)
        {
            if (need_change[i] && (cur_type == -1 || cur_type == IDA.type_ids[i]))
            {
                work_left = true;
                cur_type = IDA.type_ids[i];
                logerr("cur type %d",cur_type);
                need_change[i] = false;
                res.back().type_ids.push_back(cur_type);
                res.back().tree_ids.push_back(IDA.tree_ids[i]);
                res.back().centers_par.push_back(IDA.centers_par[i]);
                res.back().centers_self.push_back(IDA.centers_self[i]);
                res.back().transforms.push_back(IDA.transforms[i]);
            }
        }
    }
}
void BillboardCloudRaw::prepare(int branch_level, std::vector<ClusterData> &clusters, 
                                BillboardCloudData *data)
{
    std::vector<std::list<BillboardData>::iterator> out_billboards;
    for (auto &cl : clusters)
        prepare(quality, branch_level, cl, ttd, data, out_billboards);
    return;

    std::map<int, InstanceDataArrays> prev_transforms;
    std::map<int, InstanceDataArrays> all_transforms;
    std::vector<Branch> base_branches;
    BranchHeap heap;
    LeafHeap l_heap;
    for (int i = 0; i < clusters.size(); i++)
    {
        InstanceDataArrays IDA = clusters[i].IDA;
        Branch *b = clusters[i].base;
        if (!b)
            continue;
        prev_transforms.emplace(i, IDA);
    }

    //cluster can contain branches with different types. We should 
    //split it is groups with same type to make billboards correct

    for (auto &pr : prev_transforms)
    {
        std::vector<InstanceDataArrays> idas;
        split_IDA_by_type(pr.second,idas);
        for (InstanceDataArrays &group_ida : idas)
        {
            if (group_ida.type_ids.empty())
                continue;
            int k = all_transforms.size();
            all_transforms.emplace(k,group_ida);
            InstanceDataArrays IDA = clusters[pr.first].IDA;
            Branch *b = clusters[pr.first].base;
            base_branches.push_back(Branch());
            base_branches.back().deep_copy(b, heap, &l_heap);
            base_branches.back().mark_A = k;
            base_branches.back().type_id = group_ida.type_ids.front();
        }
    }
    prepare(branch_level, base_branches);
    
    std::map<int,int> proj;
    if (data)
    {
        data->level = branch_level;
        data->valid = true;
        data->atlas = *atlas;
        data->billboards.clear();
        for (auto it = all_transforms.begin(); it != all_transforms.end(); it++)
        {
            proj.emplace(it->first,data->billboards.size());
            data->billboards.push_back(BillboardData());
            data->billboards.back().IDA = it->second;
            data->billboards.back().base_position = base_branches[it->first].joints.front().pos;
        }
        std::vector<std::list<BillboardData>::iterator> iters;
        auto it = data->billboards.begin();
        while (it != data->billboards.end())
        {
            iters.push_back(it);
            it++;
        }
        for (Billboard &b : billboards)
        {
            b.instancing = true;
            iters[proj.at(b.branch_id)]->billboards.push_back(b);
        }
    }
}
void BillboardCloudRaw::prepare(int _quality, int branch_level, ClusterData &cluster, std::vector<TreeTypeData> &_ttd,
             BillboardCloudData *data, std::vector<std::list<BillboardData>::iterator> &out_billboards)
{
    quality = _quality;
    ttd = _ttd;
    
    if (!cluster.base)
    {
        logerr("cannot make billboard from cluster with empty or damaged branch");
        return;
    }
    std::map<int, InstanceDataArrays> prev_transforms;
    std::map<int, InstanceDataArrays> all_transforms;
    std::vector<Branch> base_branches;
    BranchHeap heap;
    LeafHeap l_heap;
    prev_transforms.emplace(0, cluster.IDA);

    //cluster can contain branches with different types. We should 
    //split it is groups with same type to make billboards correct

    for (auto &pr : prev_transforms)
    {
        std::vector<InstanceDataArrays> idas;
        split_IDA_by_type(pr.second,idas);
        for (InstanceDataArrays &group_ida : idas)
        {
            if (group_ida.type_ids.empty())
                continue;
            int k = all_transforms.size();
            all_transforms.emplace(k,group_ida);
            Branch *b = cluster.base;
            base_branches.push_back(Branch());
            base_branches.back().deep_copy(b, heap, &l_heap);
            base_branches.back().mark_A = k;
            base_branches.back().type_id = group_ida.type_ids.front();
        }
    }

    if (data)
    {
        if (!data->atlas.valid)
        {
            AtlasParams params;
            params.grid_x = 8;
            params.grid_y = 8;
            params.layers = 2;
            params.valid = true;
            params.x = params.grid_x*(int)quality;
            params.y = params.grid_y*(int)quality;
            TextureAtlas atl = TextureAtlas(params.x, params.y, params.layers);
            data->atlas = atl;
            data->atlas.set_grid((int)quality,(int)quality,true);
        }
        prepare(branch_level, base_branches, &(data->atlas));
    }
    else
    {
        prepare(branch_level, base_branches, nullptr);
    }
    std::map<int,int> proj;
    if (data)
    {
        data->level = branch_level;
        data->valid = true;
        if (atlas)
            data->atlas = *atlas;
        //data->billboards.clear();
        for (auto it = all_transforms.begin(); it != all_transforms.end(); it++)
        {
            proj.emplace(it->first,data->billboards.size());
            data->billboards.push_back(BillboardData());
            data->billboards.back().IDA = it->second;
            data->billboards.back().base_position = base_branches[it->first].joints.front().pos;
            auto bit = data->billboards.end();
            bit--;
            out_billboards.push_back(bit);
        }
        std::vector<std::list<BillboardData>::iterator> iters;
        auto it = data->billboards.begin();
        while (it != data->billboards.end())
        {
            iters.push_back(it);
            it++;
        }
        for (Billboard &b : billboards)
        {
            b.instancing = true;
            iters[proj.at(b.branch_id)]->billboards.push_back(b);
        }
    }
}
void BillboardCloudRaw::extend(int quality, int branch_level, ClusterData &cluster, std::vector<TreeTypeData> &_ttd,
             BillboardCloudData *data, std::vector<std::list<BillboardData>::iterator> &billboards)
{
    std::vector<InstanceDataArrays> idas;
    split_IDA_by_type(cluster.IDA,idas);
    if (idas.size() == billboards.size())
    {
        //all billboard types are already here. Just expand IDA
        for (auto &bill : billboards)
        {
            for (auto &ida : idas)
            {
                if (bill->IDA.type_ids.front() == ida.type_ids.front())
                {
                    bill->IDA = ida;
                    break;
                }
            }
        }
    }
    else
    {
        //expand billboards list with new types
        logerr("adding billboard for new type %d %d", idas.size(), billboards.size());
        InstanceDataArrays cluster_ida = cluster.IDA;
        for (auto &ida : idas)
        {
            bool found = false;
            for (auto &bill : billboards)
            {
                if (bill->IDA.type_ids.front() == ida.type_ids.front())
                {
                    bill->IDA = ida;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                //new type
                cluster.IDA = ida;
                prepare(quality, branch_level, cluster, _ttd, data, billboards);
            }
        }
        cluster.IDA = cluster_ida;
    }
}

void expand_branches(std::vector<Branch> &oldb, std::vector<Branch> &newb, int d)
{
    if (d == 0)
    {
        return;
    }
    newb.clear();
    for (Branch &branch : oldb)
    {
        for (Joint &j : branch.joints)
        {
            for (Branch *br : j.childBranches)
            {
                newb.push_back(*br);
                newb.back().mark_A = branch.mark_A;
                newb.back().type_id = branch.type_id;
            }
        }
    }
    if (d>1)
        expand_branches(newb, oldb, d-1);
}
void BillboardCloudRaw::prepare(int branch_level, std::vector<Branch> &old_branches, TextureAtlas *_atlas)
{
    if (old_branches.empty())
        return;
    int layer = old_branches.front().level;
    if (layer > branch_level)
        return;
    
    std::vector<Branch> branches;
    if (layer < branch_level)
    {
        expand_branches(old_branches,branches, branch_level - layer);
        if (branches.empty())
            branches = old_branches;
        layer = branch_level;
    }
    else if (layer == branch_level)
        branches = old_branches;
    std::vector<BillboardBox> billboard_boxes;
    for (Branch &branch : branches)
    {
        if (branch.joints.empty())
            continue;
        BBox min_bbox = get_minimal_bbox(&branch);
        vec3 base_joint = branch.joints.front().pos;
        billboard_boxes.push_back(BillboardBox(&branch, min_bbox, base_joint, -1));
    }

    int add_billboards_count = 0; //1 + layer*4;
    Visualizer tg(ttd[0].wood, ttd[0].leaf, nullptr);

    billboards.clear();
    bool external_atlas = false;
    if (!_atlas || !_atlas->is_valid())
    {
        AtlasParams params = set_atlas_params(quality, billboard_boxes.size());
        atlas = new TextureAtlas(params.x,params.y,params.layers);
        atlas->set_grid(params.grid_x,params.grid_y);
        atlas->set_clear_color(glm::vec4(0, 0, 0, 0));
        atlas->clear();
        logerr("new atlas %d",_atlas->is_valid() );
    }
    else
    {
        atlas = _atlas;
        external_atlas = true;
    }
    int atlas_capacity = atlas->capacity();
    std::vector<BranchProjectionData> projectionData;

    int i = 0;
    for (auto &p : billboard_boxes)
    {
        vec3 base_joint = p.b->joints.front().pos;
        Billboard b(p.min_bbox, 0, 0, 1, base_joint);
        vec3 plane_n = vec3(b.planeCoef.x, b.planeCoef.y, b.planeCoef.z);
        for (Joint &j : p.b->joints)
        {
            for (Branch *br : j.childBranches)
            {
                float err = projection_error_rec(br, plane_n, b.planeCoef.w);
                projectionData.push_back(BranchProjectionData(err, br, i));
                projectionData.back().br->type_id = p.b->type_id;
            }
        }
        i++;
    }
    std::sort(projectionData.begin(), projectionData.end(), BPD_comp);
    add_billboards_count = MIN(atlas_capacity - billboard_boxes.size(), projectionData.size());
    add_billboards_count = 0;
    int k = 0;
    for (auto &proj : projectionData)
    {
        if (k < add_billboards_count)
        {
            proj.br->level = -1;
            BBox min_bbox = get_minimal_bbox(proj.br);
            billboard_boxes.push_back(BillboardBox(proj.br, min_bbox, vec3(0, 0, 0), proj.parent_n));
        }
        else
        {
            break;
        }
        k++;
    }
    for (auto &p : billboard_boxes)
    {
        int num = atlas->add_tex();
        vec3 base_joint = p.b->joints.front().pos;
        Billboard b(p.min_bbox, num, p.b->mark_A, 1, base_joint);
        if (p.parent >= 0 && p.parent < billboards.size())
        {
            //it is a secondary billboard
            //it should be attached to projection of base joint
            p.b->level = layer;
            Billboard parent_billboard = billboards[p.parent];
            vec3 n = parent_billboard.planeCoef;
            vec3 proj = p.base_joint - (dot(n, p.base_joint) + parent_billboard.planeCoef.w) * n;
            p.base_joint = proj;
            b.branch_id = parent_billboard.branch_id;
        }
        create_billboard(ttd[p.b->type_id], p.b, p.min_bbox, tg, num, b, 1.5);
    }
    debugl(8,"created %d billboards\n", billboard_boxes.size());
    atlas->gen_mipmaps();
    if (external_atlas)
        atlas = nullptr;
}
BillboardCloudRenderer::BillboardCloudRenderer(BillboardCloudData *data):
Countable(3),
rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                            {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"})
{
    this->data = data;
    if (!data || !data->valid)
    {
        return;
    }
    for (BillboardData &bill : data->billboards)
    {
        Model *m = new Model();
        for (Billboard &b : bill.billboards)
        {
            b.instancing = true;
            b.to_model(m, data->atlas);
        }
        for (int i=0;i<bill.IDA.centers_par.size();i++)
        {
            if (data->level > GroveRenderer::base_level)
                bill.IDA.centers_par[i] = bill.IDA.centers_self[i];
        }

        instances.push_back(m);
    }
}
BillboardCloudRenderer::~BillboardCloudRenderer()
{
    debugl(11,"deleting BCR\n");
    if (cloud)
        delete (cloud);
    debugl(11,"cloud deleted\n");
    for (int i=0;i<instances.size();i++)
    {
        delete instances[i];
       debugl(11,"instance %d deleted\n",i);
    }
}
void BillboardCloudRenderer::render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light, 
                                    glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                                    glm::vec4 screen_size, bool to_shadow, GroveRendererDebugParams dbgpar)
{
    if (to_shadow)
        return;
    if (!data || !data->valid)
        return;
    
    if (renderMode == ONLY_SINGLE || renderMode == BOTH)
    {
        std::function<void(Model *)> _ce = [&](Model *h) {
            for (BillboardData &bill : data->billboards)
            {
                for (Billboard &b : bill.billboards)
                    b.to_model(h, data->atlas);
            }
        };
        cloud->construct(_ce);
        billboardRenderer.use();
        billboardRenderer.texture("tex", data->atlas.tex(0));
        billboardRenderer.uniform("model", cloud->model);
        billboardRenderer.uniform("projectionCamera", projection *view);
     
        cloud->render(GL_TRIANGLES);
    }
    if (renderMode == ONLY_INSTANCES || renderMode == BOTH)
    {
        billboardRendererInstancing.use();
        billboardRendererInstancing.uniform("camera_pos", camera_pos);
        billboardRendererInstancing.uniform("screen_size",screen_size);
        billboardRendererInstancing.texture("color_tex", data->atlas.tex(0));
        billboardRendererInstancing.texture("normal_tex", data->atlas.tex(1));
        billboardRendererInstancing.texture("noise",textureManager.get("noise"));
        billboardRendererInstancing.uniform("projection", projection);
        billboardRendererInstancing.uniform("view", view);
        billboardRendererInstancing.uniform("type_id", (uint)mdrd.type_id);
        billboardRendererInstancing.uniform("debug_model_id",dbgpar.need_focus_model ? dbgpar.model_focused : -1);
        
        glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                            mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
    }
}