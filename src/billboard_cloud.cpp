#include "billboard_cloud.h"
#include <algorithm>
#include <vector>
#include <glm/gtc/matrix_transform.hpp>
#include "tree.h"
#include "generated_tree.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "distribution.h"

using namespace glm;
#define TEX_ATLAS_LAYERS 4
BillboardCloudRaw::BillboardCloudRaw(int tex_w, int tex_h, std::vector<TreeTypeData> &_ttd): 
                                                       atlas(tex_w, tex_h,TEX_ATLAS_LAYERS),
                                                       wood(textureManager.empty()),
                                                       rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
                                                       billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
                                                       billboardRendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                                                                                   {"in_Position", "in_Normal", "in_Tex", "in_Model"})
{
    ttd = _ttd;
    cloud = new Model();
}
BillboardCloudRaw::~BillboardCloudRaw()
{
    if (cloud)
        delete (cloud);
    for (int i=0;i<instances.size();i++)
    {
        delete instances[i]->m;
        delete instances[i];
    }
}
void BillboardCloudRaw::setup_preparation()
{
}
void matprint()
{
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
void BillboardCloudRaw::create_billboard(std::vector<TreeTypeData> &ttd, std::map<int, InstanceDataArrays> &all_transforms,
                                         std::vector<Branch> &brs, BBox &min_bbox, Visualizer &tg, int num,
                                         Billboard &bill, float leaf_scale)
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

    mat4 tex_sh = scale(mat4(1), vec3(2, 2, 2));
    mat4 tex_tr = translate(mat4(1), vec3(-1, -1, -1));
    mat4 atlas_tr = atlas.tex_transform(num);
    mat4 result = ort * tex_tr * tex_sh * atlas_tr * SC_inv * transl * rot;
    
    Branch *branch = nullptr;
    std::function<void(Model *)> _c_wood = [&](Model *h) { if (branch) tg.recursive_branch_to_model(*branch, h, false); };
    std::function<void(Model *)> _c_leaves = [&](Model *h) { if (branch) tg.recursive_branch_to_model(*branch, h, true, leaf_scale); };

    atlas.target(num);
    rendererToTexture.use();
    for (Branch &br : brs)
    {
        Texture &wood = ttd[br.type_id].wood;
        Texture &leaf = ttd[br.type_id].leaf;
        InstanceDataArrays &ida = all_transforms.at(br.base_seg_n);
        Model bm;
        branch = &br;
        bm.construct(_c_wood);
        rendererToTexture.texture("tex", wood);
        rendererToTexture.uniform("projectionCamera", result);
        for (glm::mat4 &model : ida.transforms)
        {
            rendererToTexture.uniform("model", model);
            bm.render(GL_TRIANGLES);
        }
        bm.construct(_c_leaves);

        glBindTexture(leaf.type,leaf.texture);
        glTexParameteri(leaf.type, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(leaf.type, GL_TEXTURE_MAX_LEVEL, 0);

        rendererToTexture.texture("tex", leaf);
        rendererToTexture.uniform("projectionCamera", result);
        for (glm::mat4 &model : ida.transforms)
        {
            rendererToTexture.uniform("model", model);
            bm.render(GL_TRIANGLES);
        }

        glTexParameteri(leaf.type, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(leaf.type, GL_TEXTURE_MAX_LEVEL, 1000);
    }
    billboards.push_back(bill);
}
void BillboardCloudRaw::create_billboard(TreeTypeData &ttd, Branch *branch, BBox &min_bbox, Visualizer &tg, int num,
                                         Billboard &bill, float leaf_scale)
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

    mat4 tex_sh = scale(mat4(1), vec3(2, 2, 2));
    mat4 tex_tr = translate(mat4(1), vec3(-1, -1, -1));
    mat4 atlas_tr = atlas.tex_transform(num);
    mat4 result = ort * tex_tr * tex_sh * atlas_tr * SC_inv * transl * rot;
    Model bm;
    std::function<void(Model *)> _c_wood = [&](Model *h) { tg.recursive_branch_to_model(*branch, &bm, false); };
    std::function<void(Model *)> _c_leaves = [&](Model *h) { tg.recursive_branch_to_model(*branch, &bm, true, leaf_scale); };

    atlas.target(num);
    rendererToTexture.use();

    bm.construct(_c_wood);
    rendererToTexture.texture("tex", ttd.wood);
    rendererToTexture.uniform("model", bm.model);
    rendererToTexture.uniform("projectionCamera", result);
    bm.render(GL_TRIANGLES);

    bm.construct(_c_leaves);

    Texture &leaf = ttd.leaf;
    glBindTexture(leaf.type,leaf.texture);
    glTexParameteri(leaf.type, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(leaf.type, GL_TEXTURE_MAX_LEVEL, 0);

    rendererToTexture.texture("tex", ttd.leaf);
    rendererToTexture.uniform("model", bm.model);
    rendererToTexture.uniform("projectionCamera", result);
    bm.render(GL_TRIANGLES);

    billboards.push_back(bill);

    glTexParameteri(leaf.type, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(leaf.type, GL_TEXTURE_MAX_LEVEL, 1000);
}
BBox BillboardCloudRaw::get_minimal_bbox(Branch *branch)
{
    int iterations = 360;
    vec3 a(0, 0, 0);
    vec3 b;
    vec3 c;
    for (Segment &seg : branch->segments)
    {
        a += (seg.end - seg.begin);
    }
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
    if (branch->dead)
        return;
    for (Joint &j : branch->joints)
    {
        vec4 pos = rot * vec4(j.pos, 1);
        mn = min(mn, pos);
        mx = max(mx, pos);
        if (j.leaf && !(j.leaf->dead))
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
void BillboardCloudRaw::render(mat4 &projectionCamera)
{
    if (renderMode == ONLY_SINGLE || renderMode == BOTH)
    {
        std::function<void(Model *)> _ce = [&](Model *h) {
            for (Billboard bill : billboards)
            {
                bill.to_model(h, atlas);
            }
        };
        cloud->construct(_ce);
        billboardRenderer.use();
        billboardRenderer.texture("tex", atlas.tex());
        billboardRenderer.uniform("model", cloud->model);
        billboardRenderer.uniform("projectionCamera", projectionCamera);
        cloud->render(GL_TRIANGLES);
    }
    if (renderMode == ONLY_INSTANCES || renderMode == BOTH)
    {
        billboardRendererInstancing.use();
        billboardRendererInstancing.texture("tex", atlas.tex());
        billboardRendererInstancing.uniform("projectionCamera", projectionCamera);
        for (Instance *in : instances)
        {
            Model *m = (Model *)(in->m);
            m->update();
            in->render(GL_TRIANGLES);
        }
    }
}
void BillboardCloudRaw::set_textures(Texture _wood)
{
    wood = _wood;
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
void split_IDA_by_type(InstanceDataArrays &IDA, std::vector<InstanceDataArrays> &res)
{
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
                need_change[i] = false;
                res.back().type_ids.push_back(cur_type);
                res.back().centers_par.push_back(IDA.centers_par[i]);
                res.back().centers_self.push_back(IDA.centers_self[i]);
                res.back().transforms.push_back(IDA.transforms[i]);
            }
        }
    }
}
void BillboardCloudRaw::prepare(Tree &t, int branch_level, std::vector<ClusterData> &clusters, 
                                BillboardCloudData *data)
{
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
            base_branches.back().base_seg_n = k;
            base_branches.back().type_id = group_ida.type_ids.front();
        }
    }
    prepare(t, branch_level, base_branches);
    
    std::map<int,int> proj;
    if (data)
    {
        data->level = branch_level;
        data->valid = true;
        data->atlas = atlas;
        data->billboards.clear();
        for (auto it = all_transforms.begin(); it != all_transforms.end(); it++)
        {
            proj.emplace(it->first,data->billboards.size());
            data->billboards.push_back(BillboardData());
            data->billboards.back().IDA = it->second;
        }
    }
    for (Billboard &b : billboards)
    {
        b.instancing = true;
        Model *m = new Model();
        b.to_model(m, atlas);
        m->update();
        Instance *in = new Instance(m);
        in->addBuffer(all_transforms[b.branch_id].transforms);
        instances.push_back(in);
        if (data)
            data->billboards[proj.at(b.branch_id)].billboards.push_back(b);
    }
}
void BillboardCloudRaw::prepare(Tree &t, int branch_level, int layer)
{
    if (ready)
        return;

    std::vector<Branch> branches;
    BranchHeap heap;
    LeafHeap l_heap;
    for (Branch &b : t.branchHeaps[layer]->branches)
    {
        branches.push_back(Branch());
        branches.back().deep_copy(&b, heap, &l_heap);
    }
    prepare(t, branch_level, branches);
    ready = true;
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
                newb.back().base_seg_n = branch.base_seg_n;
                newb.back().type_id = branch.type_id;
            }
        }
    }
    if (d>1)
        expand_branches(newb, oldb, d-1);
}
void BillboardCloudRaw::prepare(Tree &t, int branch_level, std::vector<Branch> &old_branches)
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
        if (branch.dead || branch.joints.empty())
            continue;
        BBox min_bbox = get_minimal_bbox(&branch);
        vec3 base_joint = branch.joints.front().pos;
        billboard_boxes.push_back(BillboardBox(&branch, min_bbox, base_joint, -1));
    }

    int add_billboards_count = 0; //1 + layer*4;
    Visualizer tg(t.wood, t.leaf, nullptr);
    tg.set_params(t.params());
    billboards.clear();
    atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    glm::ivec4 sizes = atlas.get_sizes();
    int cnt = ceil(sqrt(billboard_boxes.size()/atlas.layers_count() + 1));
    int tex_size = (sizes.x) / cnt - 2;
    atlas.set_grid(tex_size, tex_size);
    atlas.clear();
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
                if (br->dead)
                    continue;
                float err = projection_error_rec(br, plane_n, b.planeCoef.w);
                projectionData.push_back(BranchProjectionData(err, br, i));
                projectionData.back().br->type_id = p.b->type_id;
            }
        }
        i++;
    }
    std::sort(projectionData.begin(), projectionData.end(), BPD_comp);
    add_billboards_count = MIN(atlas.layers_count() * cnt * cnt - billboard_boxes.size(), projectionData.size());
    int k = 0;
    for (auto &proj : projectionData)
    {
        if (k < add_billboards_count)
        {
            proj.br->level = -1;
            //proj.br->type_id = proj.
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
        int num = atlas.add_tex();
        vec3 base_joint = p.b->joints.front().pos;
        Billboard b(p.min_bbox, num, p.b->base_seg_n, 1, base_joint);
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
    atlas.gen_mipmaps();
}
BillboardCloudRenderer::BillboardCloudRenderer(BillboardCloudData *data):
rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                            {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"})
{
    this->data = data;
    if (!data || !data->valid)
    {
        //logerr("empty billboard data %uud",data);
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
       // m->update();
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
    debugl(11,"deleting BCR");
    if (cloud)
        delete (cloud);
    debugl(11,"cloud deleted");
    for (int i=0;i<instances.size();i++)
    {
        delete instances[i];
       debugl(11,"instance %d deleted",i);
    }
}
void BillboardCloudRenderer::render(MultiDrawRendDesc &mdrd, glm::mat4 &projectionCamera, glm::vec3 camera_pos,
                                    glm::vec4 screen_size)
{
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
        billboardRenderer.texture("tex", data->atlas.tex());
        billboardRenderer.uniform("model", cloud->model);
        billboardRenderer.uniform("projectionCamera", projectionCamera);
        cloud->render(GL_TRIANGLES);
    }
    if (renderMode == ONLY_INSTANCES || renderMode == BOTH)
    {
        billboardRendererInstancing.use();
        billboardRendererInstancing.uniform("camera_pos",camera_pos);
        billboardRendererInstancing.uniform("screen_size",screen_size);
        billboardRendererInstancing.texture("tex", data->atlas.tex());
        billboardRendererInstancing.texture("noise",textureManager.get("noise"));
        billboardRendererInstancing.uniform("projectionCamera", projectionCamera);
        billboardRendererInstancing.uniform("type_id", (uint)mdrd.type_id);
        
        glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                            mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
    }
}