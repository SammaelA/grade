#include "graphics_utils/impostor.h"
#include "graphics_utils/texture_manager.h"
#include "generation/generation_settings.h"

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
        int mult = params.monochrome ? 2 : 4;
        AtlasParams params = set_atlas_params(quality, mult*(slices_n + 1));
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
    Model br_m, l_m;
    create_models(&br, tg, params, br_m, l_m);
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
        create_billboard_model(tree_type, &br, cur, tg, num, bill, atl, params, br_m, l_m);

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
        create_billboard_model(tree_type, &br, cur, tg, num, bill, atl, params, br_m, l_m);

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