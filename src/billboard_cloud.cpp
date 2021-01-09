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

using namespace glm;
BillboardCloudRaw::BillboardCloudRaw(int tex_w, int tex_h) : atlas(tex_w, tex_h),
                                                       wood(textureManager.empty()),
                                                       rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
                                                       billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
                                                       billboardRendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                                                                                   {"in_Position", "in_Normal", "in_Tex", "in_Model"})
{
    cloud = new Model();
}
BillboardCloudRaw::~BillboardCloudRaw()
{
    delete (cloud);
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
void BillboardCloudRaw::create_billboard(Tree &t, Branch *branch, BBox &min_bbox, Visualizer &tg, int num, Billboard &bill)
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
    std::function<void(Model *)> _c_leaves = [&](Model *h) { tg.recursive_branch_to_model(*branch, &bm, true); };

    atlas.target(num);
    rendererToTexture.use();

    bm.construct(_c_wood);
    rendererToTexture.texture("tex", t.wood);
    rendererToTexture.uniform("model", bm.model);
    rendererToTexture.uniform("projectionCamera", result);
    bm.render(GL_TRIANGLES);

    bm.construct(_c_leaves);
    rendererToTexture.texture("tex", t.leaf);
    rendererToTexture.uniform("model", bm.model);
    rendererToTexture.uniform("projectionCamera", result);
    bm.render(GL_TRIANGLES);

    billboards.push_back(bill);
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
    b.x = (float)rand() / RAND_MAX;
    b.y = (float)rand() / RAND_MAX;
    b.z = (float)rand() / RAND_MAX;
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
            glm::vec2 tc = vec2(tex_c[2 * i], tex_c[2 * i + 1]);
            atlas.process_tc(id, tc);
            m->positions.push_back(v.x);
            m->positions.push_back(v.y);
            m->positions.push_back(v.z);
            m->normals.push_back(n.x);
            m->normals.push_back(n.y);
            m->normals.push_back(n.z);
            m->colors.push_back(tc.x);
            m->colors.push_back(tc.y);
            m->colors.push_back(0);
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
    if (type == 1)
    {
        vec3 npos = vec3(pos.x, pos.y, pos.z);
        npos += base_joint_rel.z * box.c;
        positions.push_back(npos);
        positions.push_back(npos + box.sizes.x * box.a);
        positions.push_back(npos + box.sizes.x * box.a + box.sizes.y * box.b);
        positions.push_back(npos + box.sizes.y * box.b);

        float d = -dot(box.c, npos);
        planeCoef = vec4(box.c.x, box.c.y, box.c.z, d);
    }
}
void BillboardCloudRaw::prepare(Tree &t, std::vector<Clusterizer::Cluster> &clusters, std::list<int> &numbers, BillboardCloudData *data)
{
    std::map<int, std::vector<glm::mat4>> all_transforms;
    std::vector<Branch> base_branches;
    BranchHeap heap;
    LeafHeap l_heap;
    for (int i : numbers)
    {
        std::vector<glm::mat4> transforms;
        Branch *b = clusters[i].prepare_to_replace(transforms);
        if (!b)
            continue;
        all_transforms.emplace(i, transforms);
        base_branches.push_back(Branch());
        base_branches.back().deep_copy(b, heap, &l_heap);
        base_branches.back().base_seg_n = i;
    }
    prepare(t, base_branches);
    if (data)
    {
        data->valid = true;
        data->atlas = &atlas;
    }
    for (Billboard &b : billboards)
    {
        b.instancing = true;
        Model *m = new Model();
        b.to_model(m, atlas);
        m->update();
        Instance *in = new Instance(m);
        in->addBuffer(all_transforms[b.branch_id]);
        instances.push_back(in);
        if (data)
        {
            data->billboards.push_back(BillboardCloudData::BillboardData());
            data->billboards.back().billboard = b;
            data->billboards.back().transforms = all_transforms[b.branch_id];
        }
    }
}
void BillboardCloudRaw::prepare(Tree &t, int layer)
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
    prepare(t, branches);
    ready = true;
}
void BillboardCloudRaw::prepare(Tree &t, std::vector<Branch> &branches)
{
    if (branches.empty())
        return;
    int layer = branches.front().level;
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
    tg.set_params(t.params);
    billboards.clear();
    atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    glm::ivec4 sizes = atlas.get_sizes();
    int cnt = ceil(sqrt(billboard_boxes.size()) + add_billboards_count);
    int tex_size = sizes.x / cnt - 2;
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
            }
        }
        i++;
    }
    std::sort(projectionData.begin(), projectionData.end(), BPD_comp);
    add_billboards_count = MIN(cnt * cnt - billboard_boxes.size(), projectionData.size());
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
        create_billboard(t, p.b, p.min_bbox, tg, num, b);
    }
    logerr("created %d billboards\n", billboard_boxes.size());
    glGenerateTextureMipmap(atlas.tex().texture);
    //atlas.gen_mipmaps();
}
BillboardCloudRenderer::BillboardCloudRenderer(BillboardCloudData *data):
rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                            {"in_Position", "in_Normal", "in_Tex", "in_Model"})
{
    this->data = data;
    if (!data || !data->valid)
    {
        logerr("empty billboard data %uud",data);
        return;
    }
    for (BillboardCloudData::BillboardData &bill : data->billboards)
    {
        bill.billboard.instancing = true;
        Model *m = new Model();
        bill.billboard.to_model(m, *data->atlas);
        m->update();
        Instance *in = new Instance(m);
        in->addBuffer(bill.transforms);
        instances.push_back(in);
    }
}
void BillboardCloudRenderer::render(glm::mat4 &projectionCamera)
{
    if (!data || !data->valid)
        return;
    
    if (renderMode == ONLY_SINGLE || renderMode == BOTH)
    {
        std::function<void(Model *)> _ce = [&](Model *h) {
            for (BillboardCloudData::BillboardData &bill : data->billboards)
            {
                bill.billboard.to_model(h, *data->atlas);
            }
        };
        cloud->construct(_ce);
        billboardRenderer.use();
        billboardRenderer.texture("tex", data->atlas->tex());
        billboardRenderer.uniform("model", cloud->model);
        billboardRenderer.uniform("projectionCamera", projectionCamera);
        cloud->render(GL_TRIANGLES);
    }
    if (renderMode == ONLY_INSTANCES || renderMode == BOTH)
    {
        billboardRendererInstancing.use();
        billboardRendererInstancing.texture("tex", data->atlas->tex());
        billboardRendererInstancing.uniform("projectionCamera", projectionCamera);
        for (Instance *in : instances)
        {
            Model *m = (Model *)(in->m);
            m->update();
            in->render(GL_TRIANGLES);
        }
    }
}