#include "impostor.h"

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

    atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    glm::ivec4 sizes = atlas.get_sizes();
    int cnt = ceil(sqrt((base_branches.size()*(slices_n + 1))/atlas.layers_count() + 1));
    int tex_size = (sizes.x) / 3;
    atlas.set_grid(tex_size, tex_size);
    atlas.clear();

    for (Branch &b : base_branches)
    {
        data->impostors.push_back(Impostor());
        make_impostor(b,data->impostors.back(),slices_n); 
    }
    data->valid = !data->impostors.empty();
    data->level = 0;
    data->atlas = atlas;
}
void ImpostorBaker::make_impostor(Branch &br, Impostor &imp, int slices_n)
{
    logerr("!!!! %d",slices_n);
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
    create_billboard(ttd[br.type_id], &br, cur, tg, num, bill);
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
        logerr("!!!! %d",num);
        bill = Billboard(cur, num, br.base_seg_n, 0, base_joint);
        create_billboard(ttd[br.type_id], &br, cur, tg, num, bill);

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
ImpostorRenderer::ImpostorRenderer(ImpostorsData *data):
impostorRenderer({"impostor_render.vs", "impostor_render.fs"}, {"in_Position", "in_Normal", "in_Tex"})
//impostorRendererInstancing({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{
    if (!data || !data->valid)
    {
        logerr("empty or corrupted impostor data");
        return;
    }
    this->data = data;

    std::vector<float> s_verts;
    for (Impostor &imp : data->impostors)
    {
        offsets.push_back(s_verts.size()/4);
        for (vec3 &pos : imp.top_slice.positions)
        {
            s_verts.push_back(pos.x);
            s_verts.push_back(pos.y);
            s_verts.push_back(pos.z);
            s_verts.push_back(1);
        }
        for (Billboard &b : imp.slices)
        {
            for (vec3 &pos : b.positions)
            {
                s_verts.push_back(pos.x);
                s_verts.push_back(pos.y);
                s_verts.push_back(pos.z);
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
void ImpostorRenderer::render(glm::mat4 &projectionCamera, glm::vec3 camera_pos, glm::vec2 LOD_min_max,
                              glm::vec4 screen_size)
{
    static float alpha = 0;
    alpha += 0.01;
    impostorRenderer.use();
    impostorRenderer.uniform("projectionCamera", projectionCamera);
    impostorRenderer.texture("tex", data->atlas.tex());
    for (int i=0;i<models.size();i++)
    {
        impostorRenderer.uniform("slice_offset",(int)offsets[i]);
        impostorRenderer.uniform("slice_verts",(int)data->impostors[i].slices[0].positions.size());
        //impostorRenderer.uniform("slice_n",j);
        vec3 center = data->impostors[i].bcyl.center;
        vec3 viewdir = center - camera_pos;
        viewdir.y = 0;
        viewdir = glm::normalize(viewdir);
        //viewdir = vec3(cos(alpha),0,sin(alpha));
        float phi = acos(dot(vec3(0,0,-1),viewdir));
        if (viewdir.x > 0)
            phi = 2*PI - phi;
        int slice_count = data->impostors[i].slices.size();
        int slice_n = (int)round(phi*slice_count/(2*PI)) % 8; 
        phi = phi - slice_n*2*PI/slice_count;
        //slogerr("phi %f d_phi %f (%f) slice_n %d %d (%f %f %f)",phi + slice_n*2*PI/slice_count, phi, 
        //dot(vec3(0,0,-1),viewdir), slice_n, slice_count, viewdir.x,viewdir.y,viewdir.z);
        mat4 rot = glm::translate(glm::rotate(glm::translate(mat4(1.0f),center), phi,vec3(0,1,0)),-center);
        impostorRenderer.uniform("rot", rot);
        for (int j=1;j<9;j++)
        {
            impostorRenderer.uniform("slice_n", slice_n + 1);
            models[i]->render(GL_TRIANGLES);
        }
    }
}