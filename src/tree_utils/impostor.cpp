#include "tree_utils/impostor.h"
#include "tinyEngine/engine.h"
#include "generation/generation_task.h"


void ImpostorBaker::prepare(Quality _quality, int branch_level, ClusterData &cluster, 
                            const std::vector<TreeTypeData> &_ttd, ImpostorsData *data, std::list<Impostor>::iterator &impostor)
{
    ImpostorGenerationParams params;
    params.quality = _quality;
    prepare(params, branch_level, cluster, _ttd, data, impostor);
}

void ImpostorBaker::prepare(ImpostorGenerationParams params, int branch_level, ClusterData &cluster, const std::vector<TreeTypeData> &_ttd,
                            ImpostorsData *data, std::list<Impostor>::iterator &impostor, int clusters_expected)
{
    quality = params.quality;
    static const int slices_n = 8;
    if (!data)
    {
        logerr("empty data cannot create impostors");
    }
    if (!data->atlas.is_valid())
    {
        int mult = clusters_expected + 1;
        //logerr("created impostors for %d clusters", clusters_expected);
        AtlasParams params = set_atlas_params(quality, mult*(slices_n + 1));
        int atlas_capacity = (params.x/params.grid_x)*(params.y/params.grid_y)*params.layers;
        TextureAtlas a = TextureAtlas(params.x,params.y,params.layers);
        data->atlas = a;
        atlas = &(data->atlas);
        atlas->set_grid(params.grid_x,params.grid_y);
        atlas->set_clear_color(float4(0, 0, 0, 0));
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
            IDA.transforms.front() = float4x4();

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
        atlas->set_clear_color(float4(0, 0, 0, 0));
        int4 sizes = atlas->get_sizes();
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
        atlas->set_clear_color(float4(0, 0, 0, 0));
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
        BBox bbox = get_bbox(&b,float3(1,0,0),float3(0,1,0),float3(0,0,1));
        int type_pos = 0;
        for (auto &tt : ttd)
        {
            if (tt.type_id == b.type_id)
            {
                break;
            }
            type_pos++;
        }
        //logerr("type id %d %d", b.type_id, type_pos);
        make_impostor(b, ttd[type_pos], *(its[proj.at(b.mark_A)]),params, data->atlas, bbox); 
    }

    data->valid = !data->impostors.empty();
    data->level = 0;
}
void ImpostorBaker::make_impostor(Branch &br, const TreeTypeData &tree_type, Impostor &imp, ImpostorGenerationParams &params, 
                                  TextureAtlas &atl, BBox &bbox)
{
    imp.bcyl.center = bbox.position + 0.5f*bbox.sizes;
    imp.bcyl.r = 0.5*sqrt(SQR(bbox.sizes.x) + SQR(bbox.sizes.z));
    imp.bcyl.h_2 = 0.5*bbox.sizes.y;

    float4x4 rot = LiteMath::rotate(float4x4(),2*PI/params.slices_n,float3(0,1,0));
    float4x4 base_rot = LiteMath::rotate(float4x4(),params.add_rotation_y,float3(0,1,0));
    float4 a = base_rot*float4(imp.bcyl.r,0,0,1);
    float3 b = float3(0,imp.bcyl.h_2,0);
    float4 c = base_rot*float4(0,0,imp.bcyl.r,1);

    BBox cur;
    float4x4 rot_inv,in_rot;
    int num;
    float3 base_joint;
    Billboard bill;
    Model br_m, l_m;
    create_models(&br, params, br_m, l_m);
    if (params.need_top_view)
    {
        cur.position = imp.bcyl.center - to_float3(a) - b - to_float3(c);

        cur.sizes = float3(2 * length(c), 2 * length(a), 2 * length(b));
        cur.a = normalize(to_float3(2.0f*c));
        cur.b = normalize(to_float3(2.0f*a));
        cur.c = normalize(2.0f*b);

        rot_inv = to_float4x4(to_float4(cur.a, 0), to_float4(cur.b, 0), to_float4(cur.c, 0), float4(0, 0, 0, 1));
        in_rot = inverse4x4(rot_inv);
        cur.position = to_float3(in_rot * to_float4(cur.position,1));

        num = atl.add_tex();
        base_joint = float3(0, 0, 0);
        bill = Billboard(cur, num, br.mark_A, 0, base_joint);
        create_billboard_model(tree_type, &br, cur, num, bill, atl, params, br_m, l_m);

        bill.positions[0] = imp.bcyl.center - to_float3(a) + to_float3(c);
        bill.positions[1] = imp.bcyl.center + to_float3(a) + to_float3(c);
        bill.positions[2] = imp.bcyl.center - to_float3(a) - to_float3(c);
        bill.positions[3] = imp.bcyl.center + to_float3(a) - to_float3(c);

        imp.top_slice = bill;
    }

    for (int i=0;i<params.slices_n;i++)
    {
        cur.position = imp.bcyl.center - to_float3(a) - b - to_float3(c);
        cur.sizes = float3(2*length(a),2*length(b),2*length(c));
        cur.a = normalize(to_float3(2.0f*a));
        cur.b = normalize(2.0f*b);
        cur.c = normalize(to_float3(2.0f*c));
        rot_inv = to_float4x4(to_float4(cur.a, 0), to_float4(cur.b, 0), to_float4(cur.c, 0), float4(0, 0, 0, 1));
        in_rot = inverse4x4(rot_inv);
        cur.position = to_float3(in_rot * to_float4(cur.position,1));
        num = atl.add_tex();
        //logerr("imp scale %f %f %f -- from %f %f %f", cur.sizes.x, cur.sizes.y, cur.sizes.z, bbox.sizes.x, bbox.sizes.y,
        //       bbox.sizes.z);
        bill = Billboard(cur, num, br.mark_A, 0, base_joint);
        create_billboard_model(tree_type, &br, cur, num, bill, atl, params, br_m, l_m);

        bill.positions[0] = imp.bcyl.center - to_float3(a) - b;
        bill.positions[1] = imp.bcyl.center + to_float3(a) - b;
        bill.positions[2] = imp.bcyl.center - to_float3(a) + b;
        bill.positions[3] = imp.bcyl.center + to_float3(a) + b;

        imp.slices.push_back(bill);

        a = rot * a;
        c = rot * c;
    }
}