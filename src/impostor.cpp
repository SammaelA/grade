#include "impostor.h"
#include "texture_manager.h"

using glm::vec3;
using glm::vec4;
using glm::mat4;

void ImpostorBaker::prepare(Tree &t, int branch_level, std::vector<Clusterizer::Cluster> &clusters, std::list<int> &numbers,
                 ImpostorsData *data)
{
    //BillboardCloudRaw::prepare(t,branch_level,clusters,numbers,data);
    if (!data)
        return;
    std::map<int, InstanceDataArrays> all_transforms;
    std::vector<Branch> base_branches;
    BranchHeap heap;
    LeafHeap l_heap;
    for (int i : numbers)
    {
        InstanceDataArrays IDA;
        Branch *b = clusters[i].prepare_to_replace(IDA);
        if (!b)
            continue;
        all_transforms.emplace(i, IDA);
        base_branches.push_back(Branch());
        base_branches.back().deep_copy(b, heap, &l_heap);
        base_branches.back().base_seg_n = i;
    }
    int slices_n = 8;
    std::map<int,int> proj;
    if (data)
    {
        data->level = branch_level;
        data->valid = true;
        data->atlas = atlas;
        data->impostors.clear();
        for (auto it = all_transforms.begin(); it != all_transforms.end(); it++)
        {
            proj.emplace(it->first,data->impostors.size());
            data->impostors.push_back(Impostor());
            data->impostors.back().IDA = it->second;
        }
    }
    atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    glm::ivec4 sizes = atlas.get_sizes();
    int cnt = ceil(sqrt((base_branches.size()*(slices_n + 1))/atlas.layers_count() + 1));
    int tex_size = (sizes.x) / 3;
    atlas.set_grid(tex_size, tex_size);
    atlas.clear();

    for (Branch &b : base_branches)
    {
        //data->impostors.push_back(Impostor());
        make_impostor(b,data->impostors[proj.at(b.base_seg_n)],slices_n); 
    }
    data->valid = !data->impostors.empty();
    data->level = 0;
    data->atlas = atlas;
}
void ImpostorBaker::make_impostor(Branch &br, Impostor &imp, int slices_n)
{
    BBox bbox = get_bbox(&br,glm::vec3(1,0,0),glm::vec3(0,1,0),glm::vec3(0,0,1));
    imp.bcyl.center = bbox.position + 0.5f*bbox.sizes;
    imp.bcyl.r = 0.5*sqrt(SQR(bbox.sizes.x) + SQR(bbox.sizes.z));
    imp.bcyl.h_2 = 0.5*bbox.sizes.y;

    glm::mat4 rot = glm::rotate(glm::mat4(1.0f),2*PI/slices_n,glm::vec3(0,1,0));
    glm::vec4 a = glm::vec4(imp.bcyl.r,0,0,1);
    glm::vec3 b = glm::vec3(0,imp.bcyl.h_2,0);
    glm::vec4 c = glm::vec4(0,0,imp.bcyl.r,1);

    BBox cur;
    cur.position = imp.bcyl.center - glm::vec3(a) - b - glm::vec3(c);

    cur.sizes = glm::vec3(2 * length(c), 2 * length(a), 2 * length(b));
    cur.a = glm::normalize(glm::vec3(2.0f*c));
    cur.b = glm::normalize(glm::vec3(2.0f*a));
    cur.c = glm::normalize(glm::vec3(2.0f*b));

    mat4 rot_inv(vec4(cur.a, 0), vec4(cur.b, 0), vec4(cur.c, 0), vec4(0, 0, 0, 1));
    mat4 in_rot = inverse(rot_inv);
    cur.position = in_rot * vec4(cur.position,1);

    int num = atlas.add_tex();
    vec3 base_joint = vec3(0, 0, 0);
    Visualizer tg(ttd[br.type_id].wood, ttd[br.type_id].leaf, nullptr);
    Billboard bill(cur, num, br.base_seg_n, 0, base_joint);
    create_billboard(ttd[br.type_id], &br, cur, tg, num, bill, 0.67);

    bill.positions[0] = imp.bcyl.center - glm::vec3(a) + glm::vec3(c);
    bill.positions[1] = imp.bcyl.center + glm::vec3(a) + glm::vec3(c);
    bill.positions[2] = imp.bcyl.center - glm::vec3(a) - glm::vec3(c);
    bill.positions[3] = imp.bcyl.center + glm::vec3(a) - glm::vec3(c);

    imp.top_slice = bill;

    for (int i=0;i<slices_n;i++)
    {
        cur.position = imp.bcyl.center - glm::vec3(a) - b - glm::vec3(c);
        cur.sizes = glm::vec3(2*length(a),2*length(b),2*length(c));
        cur.a = glm::normalize(glm::vec3(2.0f*a));
        cur.b = glm::normalize(glm::vec3(2.0f*b));
        cur.c = glm::normalize(glm::vec3(2.0f*c));
        rot_inv = mat4(vec4(cur.a, 0), vec4(cur.b, 0), vec4(cur.c, 0), vec4(0, 0, 0, 1));
        in_rot = inverse(rot_inv);
        cur.position = in_rot * vec4(cur.position,1);
        num = atlas.add_tex();

        bill = Billboard(cur, num, br.base_seg_n, 0, base_joint);
        create_billboard(ttd[br.type_id], &br, cur, tg, num, bill, 0.67);

        bill.positions[0] = imp.bcyl.center - glm::vec3(a) - b;
        bill.positions[1] = imp.bcyl.center + glm::vec3(a) - b;
        bill.positions[2] = imp.bcyl.center - glm::vec3(a) + b;
        bill.positions[3] = imp.bcyl.center + glm::vec3(a) + b;

        imp.slices.push_back(bill);

        a = rot * a;
        c = rot * c;
    }
    glGenerateTextureMipmap(atlas.tex().texture);
}

void ImpostorBaker::prepare_all_grove(Tree &t, GroveGenerationData &ggd, int branch_level, std::vector<Clusterizer::Cluster> &clusters,
                           std::list<int> &numbers, ImpostorsData *data)
{
    if (!data)
        return;
    
    data->level = branch_level;
    data->valid = true;
    data->atlas = atlas;
    data->impostors.clear();
    data->impostors.push_back(Impostor());
    data->impostors.back().id = 0;

    InstanceDataArrays ida;
    ida.centers_par = {ggd.pos};
    ida.centers_self = {ggd.pos};
    ida.transforms = {glm::mat4(1.0f)};

    data->impostors.back().IDA = ida;

    std::map<int, InstanceDataArrays> all_transforms;
    std::vector<Branch> base_branches;
    BranchHeap heap;
    LeafHeap l_heap;
        for (int i : numbers)
    {
        InstanceDataArrays IDA;
        Branch *b = clusters[i].prepare_to_replace(IDA);
        if (!b)
            continue;
        all_transforms.emplace(i, IDA);
        base_branches.push_back(Branch());
        base_branches.back().deep_copy(b, heap, &l_heap);
        base_branches.back().base_seg_n = i;
    }
    int slices_n = 8;

    atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    glm::ivec4 sizes = atlas.get_sizes();
    int tex_size = (sizes.x)/3;
    atlas.set_grid(tex_size, tex_size);
    atlas.clear();

    BBox bbox = BBox();
    bbox.a = glm::vec3(1,0,0);
    bbox.b = glm::vec3(0,1,0);
    bbox.c = glm::vec3(0,0,1);
    float y_size = MAX(MAX(2.0f*ggd.size.x, ggd.size.y), 2.0f*ggd.size.z);
    float y_thr = 10;
    bbox.sizes = vec3(2.0f*ggd.size.x, y_size + y_thr, 2.0f*ggd.size.z);
    bbox.position = ggd.pos - vec3(ggd.size.x, -y_thr, ggd.size.z);
    Impostor &imp = data->impostors.back();

    imp.bcyl.center = bbox.position + 0.5f*bbox.sizes;
    imp.bcyl.r = 0.5*sqrt(SQR(bbox.sizes.x) + SQR(bbox.sizes.z));
    imp.bcyl.h_2 = 0.5*bbox.sizes.y;

    glm::mat4 rot = glm::rotate(glm::mat4(1.0f),2*PI/slices_n,glm::vec3(0,1,0));
    glm::vec4 a = glm::vec4(imp.bcyl.r,0,0,1);
    glm::vec3 b = glm::vec3(0,imp.bcyl.h_2,0);
    glm::vec4 c = glm::vec4(0,0,imp.bcyl.r,1);

    BBox cur;
    cur.position = imp.bcyl.center - glm::vec3(a) - b - glm::vec3(c);

    cur.sizes = glm::vec3(2 * length(c), 2 * length(a), 2 * length(b));
    cur.a = glm::normalize(glm::vec3(2.0f*c));
    cur.b = glm::normalize(glm::vec3(2.0f*a));
    cur.c = glm::normalize(glm::vec3(2.0f*b));

    mat4 rot_inv(vec4(cur.a, 0), vec4(cur.b, 0), vec4(cur.c, 0), vec4(0, 0, 0, 1));
    mat4 in_rot = inverse(rot_inv);
    cur.position = in_rot * vec4(cur.position,1);

    int num = atlas.add_tex();
    vec3 base_joint = vec3(0, 0, 0);
    Visualizer tg(ttd[base_branches.front().type_id].wood, ttd[base_branches.front().type_id].leaf, nullptr);
    Billboard bill(cur, num, 0, 0, base_joint);
    create_billboard(ttd, all_transforms, base_branches, cur, tg, num, bill, 0.8);

    bill.positions[0] = imp.bcyl.center - glm::vec3(a) + glm::vec3(c);
    bill.positions[1] = imp.bcyl.center + glm::vec3(a) + glm::vec3(c);
    bill.positions[2] = imp.bcyl.center - glm::vec3(a) - glm::vec3(c);
    bill.positions[3] = imp.bcyl.center + glm::vec3(a) - glm::vec3(c);

    imp.top_slice = bill;

    for (int i=0;i<slices_n;i++)
    {
        cur.position = imp.bcyl.center - glm::vec3(a) - b - glm::vec3(c);
        cur.sizes = glm::vec3(2*length(a),2*length(b),2*length(c));
        cur.a = glm::normalize(glm::vec3(2.0f*a));
        cur.b = glm::normalize(glm::vec3(2.0f*b));
        cur.c = glm::normalize(glm::vec3(2.0f*c));
        rot_inv = mat4(vec4(cur.a, 0), vec4(cur.b, 0), vec4(cur.c, 0), vec4(0, 0, 0, 1));
        in_rot = inverse(rot_inv);
        cur.position = in_rot * vec4(cur.position,1);
        num = atlas.add_tex();

        bill = Billboard(cur, num, 0, 0, base_joint);
        create_billboard(ttd, all_transforms, base_branches, cur, tg, num, bill, 0.8);

        bill.positions[0] = imp.bcyl.center - glm::vec3(a) - b;
        bill.positions[1] = imp.bcyl.center + glm::vec3(a) - b;
        bill.positions[2] = imp.bcyl.center - glm::vec3(a) + b;
        bill.positions[3] = imp.bcyl.center + glm::vec3(a) + b;

        imp.slices.push_back(bill);

        a = rot * a;
        c = rot * c;
    }
    glGenerateTextureMipmap(atlas.tex().texture);

    
    data->valid = !data->impostors.empty();
    data->level = 0;
    data->atlas = atlas;
}

ImpostorRenderer::ImpostorRenderer(ImpostorsData *data):
impostorRenderer({"impostor_render.vs", "impostor_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
impostorRendererInstancing({"impostor_render_instancing.vs", "impostor_render_instancing.fs"}, {"in_Position", "in_Normal", "in_Tex"})
//impostorRendererInstancing({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{
    if (!data || !data->valid)
    {
        logerr("empty or corrupted impostor data");
        return;
    }
    this->data = data;
    if (data->impostors.size() == 1)
    {
        hth = 0.1;
        vth = 0.5;
    }
    else
    {
        hth = 0;
        vth = 0.1;
    }
    std::vector<float> s_verts;
    for (Impostor &imp : data->impostors)
    {
        offsets.push_back(s_verts.size()/8);
        std::vector<glm::vec3> tcs = imp.top_slice.get_tc(data->atlas);
        for (int i = 0;i<MIN(imp.top_slice.positions.size(),tcs.size());i++)
        {
            vec3 &pos = imp.top_slice.positions[i];
            
            s_verts.push_back(pos.x);
            s_verts.push_back(pos.y);
            s_verts.push_back(pos.z);
            s_verts.push_back(1);

            s_verts.push_back(tcs[i].x);
            s_verts.push_back(tcs[i].y);
            s_verts.push_back(tcs[i].z);
            s_verts.push_back(1);
        }
        for (Billboard &b : imp.slices)
        {
            tcs = b.get_tc(data->atlas);
            for (int i = 0;i<MIN(b.positions.size(),tcs.size());i++)
            {
                vec3 &pos = b.positions[i];

                s_verts.push_back(pos.x);
                s_verts.push_back(pos.y);
                s_verts.push_back(pos.z);
                s_verts.push_back(1);
                            
                s_verts.push_back(tcs[i].x);
                s_verts.push_back(tcs[i].y);
                s_verts.push_back(tcs[i].z);
                s_verts.push_back(1);
            }
        }
        Model *bm = new Model();
        std::vector<float> vertexes = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
        std::vector<float> tc = {0,0,0,0, 1,0,0,0, 0,1,0,0, 1,1,0,0};
        std::vector<GLuint> indices = {0, 1, 3, 2, 0, 3};

        std::function<void(Model *)> _c_mip = [&](Model *h) {
            bm->positions = vertexes;
            bm->colors = tc;
            bm->indices = indices;
        };
        bm->construct(_c_mip);
        models.push_back(bm);
    }
    glGenBuffers(1, &slicesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slicesBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*s_verts.size(), s_verts.data(), GL_DYNAMIC_DRAW);


}
void ImpostorRenderer::render(std::vector<uint> &counts, glm::mat4 &projectionCamera, glm::vec3 camera_pos, glm::vec2 LOD_min_max,
                              glm::vec4 screen_size)
{
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slicesBuffer);
    impostorRendererInstancing.use();
    impostorRendererInstancing.uniform("projectionCamera", projectionCamera);
    impostorRendererInstancing.texture("tex", data->atlas.tex());
    impostorRendererInstancing.uniform("camera_pos", camera_pos);
    impostorRendererInstancing.uniform("screen_size", screen_size);
    impostorRendererInstancing.texture("noise", textureManager.get("noise"));
    impostorRendererInstancing.uniform("hor_vert_transition_thr", glm::vec2(hth,vth));

    for (int i=0;i<models.size();i++)
    {
        float a_step = (2*PI)/(data->impostors[i].slices.size());

        impostorRendererInstancing.uniform("slice_offset",(int)offsets[i]);
        impostorRendererInstancing.uniform("slice_verts",(int)data->impostors[i].slices[0].positions.size());
        impostorRendererInstancing.uniform("slice_count", (int)(data->impostors[i].slices.size()));
        impostorRendererInstancing.uniform("id",(uint)data->impostors[i].id);
        impostorRendererInstancing.uniform("imp_center", data->impostors[i].bcyl.center);
        impostorRendererInstancing.uniform("angle_step",a_step);
        impostorRendererInstancing.uniform("delta",0.5f*a_step);

        if (data->impostors[i].id >= 0 && data->impostors[i].id < counts.size() && counts[data->impostors[i].id] > 0)
        {
            Model *m = models[i];
            if(m && m->indexed)
            {
                glBindVertexArray(m->vao);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m->ibo);
                glDrawElementsInstanced(GL_TRIANGLES, m->SIZE, GL_UNSIGNED_INT, 0, counts[data->impostors[i].id]);
            }
        }
    }
}