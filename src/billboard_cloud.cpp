#include "billboard_cloud.h"
#include <algorithm>
#include <vector>
#include <glm/gtc/matrix_transform.hpp>
#include "tree.h"
#include "generated_tree.h"
#include "tinyEngine/utility/shader.h"
#include "tinyEngine/utility.h"
#include "visualizer.h"
using namespace glm;
BillboardCloud::BillboardCloud(int tex_w, int tex_h):
atlas(tex_w,tex_h),
rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{
    cloud = new Model();
}
BillboardCloud::~BillboardCloud()
{
    delete(cloud);
}
void BillboardCloud::setup_preparation()
{

}
void matprint()
{

}


bool BillboardCloud::BPD_comp(BranchProjectionData &a, BranchProjectionData &b)
{
    return a.projection_err > b.projection_err;   
}
float BillboardCloud::projection_error_rec(Branch *b, vec3 &n, float d)
{
    if (!b || b->joints.size() == 0)
        return 0;
    float err = 0.0;
    for (auto &j : b->joints)
    {
        err += abs(dot(j.pos,n) + d);
        for (auto br : j.childBranches)
            err += projection_error_rec(br,n,d);
    }
    return err;
}
void BillboardCloud::create_billboard(Tree &t, Branch *branch, BBox &min_bbox, Visualizer &tg, int num, Billboard &bill)
{
    if (num < 0)
    {
        fprintf(stderr, "too many billboards = %d", billboard_count);
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
    if (t.wood)
        rendererToTexture.texture("tex", *(t.wood));
    rendererToTexture.uniform("model", bm.model);
    rendererToTexture.uniform("projectionCamera", result);
    bm.render(GL_TRIANGLES);

    bm.construct(_c_leaves);
    if (t.leaf)
        rendererToTexture.texture("tex", *(t.leaf));
    rendererToTexture.uniform("model", bm.model);
    rendererToTexture.uniform("projectionCamera", result);
    bm.render(GL_TRIANGLES);

    billboards.push_back(bill);
}
void BillboardCloud::prepare(Tree &t, int layer)
{
    if (ready)
       return;

    std::vector<BillboardBox> billboard_boxes;
    for (Branch &branch : t.branchHeaps[layer]->branches)
    {
        if (branch.dead || branch.joints.empty())
            continue;
        BBox min_bbox = get_minimal_bbox(&branch);
        vec3 base_joint = branch.joints.front().pos;
        billboard_boxes.push_back(BillboardBox(&branch,min_bbox,base_joint,-1));
    }    
    
    int add_billboards_count = 1 + layer*4;
    Visualizer tg(t.wood, t.leaf, nullptr);
    tg.set_params(t.params);
    billboards.clear();
    atlas.set_clear_color(glm::vec4(0,0,0,0));
    glm::ivec4 sizes = atlas.get_sizes();
    int cnt = ceil(sqrt(billboard_boxes.size()) + add_billboards_count);
    int tex_size = sizes.x/cnt-2;
    atlas.set_grid(tex_size,tex_size);
    atlas.clear();
    std::vector<BranchProjectionData> projectionData;

    int i = 0;
    for (auto &p : billboard_boxes)
    {
        vec3 base_joint = p.b->joints.front().pos;
        Billboard b(p.min_bbox, 0, 1, base_joint);
        vec3 plane_n = vec3(b.planeCoef.x,b.planeCoef.y,b.planeCoef.z);
        for (Joint &j : p.b->joints)
        {
            for (Branch *br : j.childBranches)
            {
                if (br->dead)
                    continue;
                float err = projection_error_rec(br,plane_n,b.planeCoef.w);
                projectionData.push_back(BranchProjectionData(err,br,i));
            }
        }
        i++;
    }
    std::sort(projectionData.begin(),projectionData.end(),BPD_comp);
    add_billboards_count = MIN(cnt*cnt - billboard_boxes.size(), projectionData.size());
    int k = 0;
    for (auto &proj : projectionData)
    {
        if (k<add_billboards_count)
        {
            proj.br->level = -1;
            BBox min_bbox = get_minimal_bbox(proj.br);
            billboard_boxes.push_back(BillboardBox(proj.br, min_bbox,vec3(0,0,0),proj.parent_n));
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
        Billboard b(p.min_bbox, num, 1, base_joint);
        if (p.parent>=0 && p.parent < billboards.size())
        {
            //it is a secondary billboard
            //it should be attached to projection of base joint
            p.b->level = layer;
            Billboard parent_billboard = billboards[p.parent];
            vec3 n = parent_billboard.planeCoef;
            vec3 proj = p.base_joint - (dot(n,p.base_joint) + parent_billboard.planeCoef.w)*n;
            p.base_joint = proj;
        }
        create_billboard(t,p.b,p.min_bbox,tg,num,b);
    }
    glGenerateTextureMipmap(atlas.tex().texture);
    //atlas.gen_mipmaps();
    ready = true;
}
BillboardCloud::BBox BillboardCloud::get_minimal_bbox(Branch *branch)
{
    int iterations = 360;
    vec3 a(0,0,0);
    vec3 b;
    vec3 c;
    for (Segment &seg : branch->segments)
    {
        a+=(seg.end - seg.begin);
    }
    a = normalize(a);//average branch direction
    b.x = (float)rand()/RAND_MAX;
    b.y = (float)rand()/RAND_MAX;
    b.z = (float)rand()/RAND_MAX;
    b = normalize(b);
    c = cross(a,b);
    mat4 br = rotate(mat4(1.0), (float)(2*PI/iterations), a);
    BBox min_bbox;
    float min_minside = 1e10;
    for (int i=0;i<iterations;i++)
    {
        BBox box = get_bbox(branch,a,b,c);
        float minside = MIN(MIN(box.sizes.x,box.sizes.y),box.sizes.z);
        if (minside<min_minside)
        {
            min_bbox = box;
            min_minside = minside;
        }
        b = br*vec4(b,1);
        c = cross(a,b);
    }
    return min_bbox;
}
BillboardCloud::BBox BillboardCloud::get_bbox(Branch *branch, glm::vec3 a, glm::vec3 b, glm::vec3 c)
{   
    vec4 bias = vec4(1,1,1,0);
    mat4 rot_inv(vec4(a,0),vec4(b,0),vec4(c,0),vec4(0,0,0,1));
    mat4 rot = inverse(rot_inv);
    //transform from model to bbox coordinates
    vec4 mx(-1e10,-1e10,-1e10,1);
    vec4 mn(1e10,1e10,1e10,1);
    BBox box;
    update_bbox(branch,rot,mn,mx);
    mn-=bias;
    mx+=bias;
    box.sizes = mx - mn;
    box.position = mn;
    box.a = a;
    box.b = b;
    box.c = c;
    //fprintf(stderr,"pos size (%f %f %f) (%f %f %f)",mn.x,mn.y,mn.z,box.sizes.x,box.sizes.y,box.sizes.z);
    return box;
}
void BillboardCloud::update_bbox(Branch *branch, mat4 &rot, vec4 &mn, vec4 &mx)
{
    if (branch->dead)
        return;
    for (Joint &j : branch->joints)
    {
        vec4 pos = rot*vec4(j.pos,1);
        mn = min(mn,pos);
        mx = max(mx,pos);
        if (j.leaf && !(j.leaf->dead))
        {
            for (auto &vert : j.leaf->edges)
            {
                pos = rot*vec4(vert,1);
                mn = min(mn,pos);
                mx = max(mx,pos);
            }
        }
        for (Branch *br : j.childBranches)
            update_bbox(br,rot,mn,mx);
    }
}
void BillboardCloud::render(mat4 &projectionCamera)
{
    std::function<void(Model*)> _ce = [&](Model* h)
    {
        for (Billboard bill : billboards)
        {
            bill.to_model(h,atlas);
        }
    };
    cloud->construct(_ce);
    billboardRenderer.use();
    billboardRenderer.texture("tex",atlas.tex());
    billboardRenderer.uniform("model", cloud->model);
    billboardRenderer.uniform("projectionCamera", projectionCamera);
    cloud->render(GL_TRIANGLES);
}
void BillboardCloud::set_textures(Texture *wood)
{
    pwood = wood;
}
void BillboardCloud::Billboard::to_model(Model *m, TextureAtlas &atlas)
{
    if (positions.size() == 4)
    {
        int _b = m->positions.size()/3;
        glm::vec3 a = positions[0];
        glm::vec3 b = positions[1];
        glm::vec3 c = positions[2];
        glm::vec3 n = glm::normalize(glm::cross(a-b,c-b));
        std::vector<float> tex_c{0,0,1,0,1,1,0,1};
        for (int i=0;i<4;i++)
        {
            glm::vec3 v = positions[i];
            glm::vec2 tc = vec2(tex_c[2*i],tex_c[2*i+1]);
            atlas.process_tc(id,tc);
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
        m->indices.push_back(_b+1);
        m->indices.push_back(_b+3);
        m->indices.push_back(_b+1);
        m->indices.push_back(_b+2);
        m->indices.push_back(_b+3);
    }
}
BillboardCloud::Billboard::Billboard(const BBox &box, int id, int type, glm::vec3 base_joint)
{
    this->id = id;
    mat4 rot_inv(vec4(box.a,0),vec4(box.b,0),vec4(box.c,0),vec4(0,0,0,1));
    mat4 rot = inverse(rot_inv);
    vec3 base_joint_rel = rot*vec4(base_joint,1.0f);
    base_joint_rel -= box.position;
    if (base_joint_rel.x > box.sizes.x || base_joint_rel.y > box.sizes.y || base_joint_rel.z > box.sizes.z ||
        base_joint_rel.x < 0 || base_joint_rel.y < 0 || base_joint_rel.z < 0 || true)
    {
        //fprintf(stderr,"wrong %f %f %f sizes + %f %f %f %f\n",base_joint_rel.x,base_joint_rel.y,base_joint_rel.z,
        //box.sizes.x,box.sizes.y,box.sizes.z, box.sizes.y/box.sizes.z);
    }
    vec4 pos = rot_inv*vec4(box.position,1.0f);
    if (type == 1)
    {
        vec3 npos = vec3(pos.x,pos.y,pos.z);
        //if (base_joint_rel.x > 10)
        //    return;//broken billboard TODO: fix me
        npos += base_joint_rel.z*box.c;
        positions.push_back(npos);
        positions.push_back(npos+box.sizes.x*box.a);
        positions.push_back(npos+box.sizes.x*box.a+box.sizes.y*box.b);
        positions.push_back(npos+box.sizes.y*box.b);

        float d = -dot(box.c,npos);
        planeCoef = vec4(box.c.x,box.c.y,box.c.z,d);
    }
}