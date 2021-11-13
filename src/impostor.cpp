#include "impostor.h"
#include "graphics_utils/texture_manager.h"

using glm::vec3;
using glm::vec4;
using glm::mat4;

void ImpostorBaker::prepare(Quality _quality, int branch_level, ClusterData &cluster, 
                            std::vector<TreeTypeData> &_ttd, ImpostorsData *data, std::list<Impostor>::iterator &impostor)
{
    ImpostorGenerationParams params;
    params.quality = _quality;
    prepare(params, branch_level, cluster, _ttd, data, impostor);
}

void ImpostorBaker::prepare(ImpostorGenerationParams params, int branch_level, ClusterData &cluster, std::vector<TreeTypeData> &_ttd,
                            ImpostorsData *data, std::list<Impostor>::iterator &impostor)
{
    quality = params.quality;
    static const int slices_n = 8;
    if (!data)
    {
        logerr("empty data cannot create impostors");
    }
    if (!data->atlas.is_valid())
    {
        AtlasParams params = set_atlas_params(quality, (slices_n + 1)*4);
        int atlas_capacity = (params.x/params.grid_x)*(params.y/params.grid_y)*params.layers;
        TextureAtlas a = TextureAtlas(params.x,params.y,params.layers);
        data->atlas = a;
        atlas = &(data->atlas);
        atlas->set_grid(params.grid_x,params.grid_y);
        atlas->set_clear_color(glm::vec4(0, 0, 0, 0));
    }
    else
    {
        atlas = &(data->atlas);
    }

    ttd = _ttd;
    std::vector<ClusterData> clusters = {cluster};
    prepare(branch_level, clusters, params, data);
    auto it = data->impostors.end();
    it--;
    impostor = it;
    atlas = nullptr;
}

void ImpostorBaker::prepare(int branch_level, std::vector<ClusterData> &clusters, ImpostorGenerationParams &params,
                            ImpostorsData *data)
{
    if (!data)
        return;
    std::map<int, InstanceDataArrays> all_transforms;
    std::vector<Branch> base_branches;
    BranchHeap heap;
    LeafHeap l_heap;
    for (int i = 0; i < clusters.size(); i++)
    {
        InstanceDataArrays IDA = clusters[i].IDA;
        Branch *b = clusters[i].base;
        if (IDA.transforms.size() == 1)
            IDA.transforms.front() = glm::mat4(1.0f);

        if (!b)
            continue;
        all_transforms.emplace(i, IDA);
        base_branches.push_back(Branch());
        base_branches.back().deep_copy(b, heap, &l_heap);
        base_branches.back().mark_A = i;
    }
    int slices_n = 8;
    std::map<int,int> proj;
    if (data)
    {
        data->level = branch_level;
        data->valid = true;
        for (auto it = all_transforms.begin(); it != all_transforms.end(); it++)
        {
            proj.emplace(it->first,data->impostors.size());
            data->impostors.push_back(Impostor());
            data->impostors.back().IDA = it->second;
        }
    }
    if (atlas)
    {
        /*
        atlas->set_clear_color(glm::vec4(0, 0, 0, 0));
        glm::ivec4 sizes = atlas->get_sizes();
        int cnt = ceil(sqrt((base_branches.size()*(slices_n + 1))/atlas->layers_count() + 1));
        int tex_size = (sizes.x) / cnt;
        atlas->set_grid(tex_size, tex_size);
        atlas->clear();
        */
    }
    else
    {
        AtlasParams params = set_atlas_params(quality, (slices_n + 1)*base_branches.size());
        int atlas_capacity = (params.x/params.grid_x)*(params.y/params.grid_y)*params.layers;
        atlas = new TextureAtlas(params.x,params.y,params.layers);
        atlas->set_grid(params.grid_x,params.grid_y);
        atlas->set_clear_color(glm::vec4(0, 0, 0, 0));
        atlas->clear();
    }
    std::vector<std::list<Impostor>::iterator> its;
    std::list<Impostor>::iterator it = data->impostors.begin();
    while (it != data->impostors.end())
    {
        its.push_back(it);
        it++;
    }
    data->atlas = *atlas;
    atlas = nullptr;
    for (Branch &b : base_branches)
    {
        BBox bbox = get_bbox(&b,glm::vec3(1,0,0),glm::vec3(0,1,0),glm::vec3(0,0,1));
        make_impostor(b, ttd[b.type_id], *(its[proj.at(b.mark_A)]),params, data->atlas, bbox); 
    }

    data->valid = !data->impostors.empty();
    data->level = 0;
}
void ImpostorBaker::make_impostor(Branch &br, TreeTypeData &tree_type, Impostor &imp, ImpostorGenerationParams &params, 
                                  TextureAtlas &atl, BBox &bbox)
{
    imp.bcyl.center = bbox.position + 0.5f*bbox.sizes;
    imp.bcyl.r = 0.5*sqrt(SQR(bbox.sizes.x) + SQR(bbox.sizes.z));
    imp.bcyl.h_2 = 0.5*bbox.sizes.y;

    glm::mat4 rot = glm::rotate(glm::mat4(1.0f),2*PI/params.slices_n,glm::vec3(0,1,0));
    glm::mat4 base_rot = glm::rotate(glm::mat4(1.0f),params.add_rotation_y,glm::vec3(0,1,0));
    glm::vec4 a = base_rot*glm::vec4(imp.bcyl.r,0,0,1);
    glm::vec3 b = glm::vec3(0,imp.bcyl.h_2,0);
    glm::vec4 c = base_rot*glm::vec4(0,0,imp.bcyl.r,1);

    BBox cur;
    mat4 rot_inv,in_rot;
    int num;
    vec3 base_joint;
    Billboard bill;
    Visualizer tg;
    if (params.need_top_view)
    {
        cur.position = imp.bcyl.center - glm::vec3(a) - b - glm::vec3(c);

        cur.sizes = glm::vec3(2 * length(c), 2 * length(a), 2 * length(b));
        cur.a = glm::normalize(glm::vec3(2.0f*c));
        cur.b = glm::normalize(glm::vec3(2.0f*a));
        cur.c = glm::normalize(glm::vec3(2.0f*b));

        rot_inv = mat4(vec4(cur.a, 0), vec4(cur.b, 0), vec4(cur.c, 0), vec4(0, 0, 0, 1));
        in_rot = inverse(rot_inv);
        cur.position = in_rot * vec4(cur.position,1);

        num = atl.add_tex();
        base_joint = vec3(0, 0, 0);
        tg = Visualizer(tree_type.wood, tree_type.leaf, nullptr);
        bill = Billboard(cur, num, br.mark_A, 0, base_joint);
        create_billboard(tree_type, &br, cur, tg, num, bill, atl, params);

        bill.positions[0] = imp.bcyl.center - glm::vec3(a) + glm::vec3(c);
        bill.positions[1] = imp.bcyl.center + glm::vec3(a) + glm::vec3(c);
        bill.positions[2] = imp.bcyl.center - glm::vec3(a) - glm::vec3(c);
        bill.positions[3] = imp.bcyl.center + glm::vec3(a) - glm::vec3(c);

        imp.top_slice = bill;
    }

    for (int i=0;i<params.slices_n;i++)
    {
        cur.position = imp.bcyl.center - glm::vec3(a) - b - glm::vec3(c);
        cur.sizes = glm::vec3(2*length(a),2*length(b),2*length(c));
        cur.a = glm::normalize(glm::vec3(2.0f*a));
        cur.b = glm::normalize(glm::vec3(2.0f*b));
        cur.c = glm::normalize(glm::vec3(2.0f*c));
        rot_inv = mat4(vec4(cur.a, 0), vec4(cur.b, 0), vec4(cur.c, 0), vec4(0, 0, 0, 1));
        in_rot = inverse(rot_inv);
        cur.position = in_rot * vec4(cur.position,1);
        num = atl.add_tex();

        bill = Billboard(cur, num, br.mark_A, 0, base_joint);
        create_billboard(tree_type, &br, cur, tg, num, bill, atl, params);

        bill.positions[0] = imp.bcyl.center - glm::vec3(a) - b;
        bill.positions[1] = imp.bcyl.center + glm::vec3(a) - b;
        bill.positions[2] = imp.bcyl.center - glm::vec3(a) + b;
        bill.positions[3] = imp.bcyl.center + glm::vec3(a) + b;

        imp.slices.push_back(bill);

        a = rot * a;
        c = rot * c;
    }
}

void ImpostorBaker::prepare_all_grove(GroveGenerationData &ggd, int branch_level, std::vector<ClusterData> &clusters,
                                      ImpostorsData *data)
{
    if (!data)
        return;

    data->level = branch_level;
    data->valid = true;
    data->atlas = *atlas;
    data->impostors.clear();
    data->impostors.push_back(Impostor());
    data->impostors.back().id = 0;

    InstanceDataArrays ida;
    ida.centers_par = {ggd.pos};
    ida.centers_self = {ggd.pos};
    ida.transforms = {glm::mat4(1.0f)};

    data->impostors.back().IDA = ida;

    return;
}
ImpostorRenderer::~ImpostorRenderer()
{
    #define DELBUF(a) if (a) { glDeleteBuffers(1, &(a)); a = 0;}

    DELBUF(slicesBuffer);
    DELBUF(impostorsDataBuffer);

    for (int i=0;i<models.size();i++)
    {
        if (models[i])
        {
            delete models[i];
            models[i] = nullptr;
        }
    }
}
ImpostorRenderer::ImpostorRenderer(ImpostorsData *data):
impostorRenderer({"impostor_render.vs", "impostor_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
impostorRendererInstancing({"impostor_render_instancing.vs", "impostor_render_instancing.fs"}, {"in_Position", "in_Normal", "in_Tex"})
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

        std::function<void(Model *)> _c_mip = [&](Model *h) 
        {
            bm->positions = vertexes;
            bm->colors = tc;
            bm->indices = indices;
        };
        bm->construct(_c_mip);
        models.push_back(bm);
    }
    std::vector<ImpostorData> imp_data_buffer;
    auto it = data->impostors.begin();
    for (int i=0;i<models.size();i++)
    {
        ImpostorData dat;
        dat.slice_count = (it->slices.size());
        dat.slice_offset = offsets[i];
        dat.slice_verts = it->slices[0].positions.size();
        dat.imp_center = glm::vec4(it->bcyl.center,1);
        imp_data_buffer.push_back(dat);
        it++;
    }
    glGenBuffers(1, &impostorsDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, impostorsDataBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(ImpostorData)*imp_data_buffer.size(), imp_data_buffer.data(), GL_STATIC_DRAW);
    glGenBuffers(1, &slicesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slicesBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*s_verts.size(), s_verts.data(), GL_STATIC_DRAW);


}
void ImpostorRenderer::render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light,
                              glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                              glm::vec4 screen_size, bool to_shadow, GroveRendererDebugParams dbgpar)
{
    if (to_shadow)//we do not render impostors to shadow
        return;
    if (!data || !data->valid)
        return;
    
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slicesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, impostorsDataBuffer);

    impostorRendererInstancing.use();
    impostorRendererInstancing.uniform("projection", projection);
    impostorRendererInstancing.uniform("view", view);
    impostorRendererInstancing.texture("color_tex", data->atlas.tex(0));
    impostorRendererInstancing.texture("normal_tex", data->atlas.tex(1));
    impostorRendererInstancing.uniform("camera_pos", camera_pos);
    impostorRendererInstancing.uniform("screen_size", screen_size);
    impostorRendererInstancing.texture("noise", textureManager.get("noise"));
    impostorRendererInstancing.uniform("hor_vert_transition_thr", glm::vec2(hth,vth));
    impostorRendererInstancing.uniform("delta",0.5f);
    impostorRendererInstancing.uniform("type_id",(uint)mdrd.type_id);
    impostorRendererInstancing.uniform("vertex_id_offset",mdrd.base_vertex_id);
    impostorRendererInstancing.uniform("debug_model_id",dbgpar.need_focus_model ? dbgpar.model_focused : -1);
    glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                        mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
}