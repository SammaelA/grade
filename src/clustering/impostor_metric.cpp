#include "impostor_metric.h"
#include "dist_data_table.h"
#include "graphics_utils/billboard_cloud.h"
#include "graphics_utils/texture_manager.h"
#include "impostor_similarity_params.h"

ImpostorSimilarityParams isimParams;

void BranchClusteringDataImpostor::clear()
{

}
void ImpostorClusteringHelper::clear_branch_data(BranchClusteringData *base, ClusteringContext *ictx)
{
    BranchClusteringDataImpostor *id = dynamic_cast<BranchClusteringDataImpostor *>(base);

    if (id && ictx)
    {
        for (auto &it : id->self_impostor->slices)
        {
            ictx->self_impostors_data->atlas.remove_tex(it.id);
        }   
        ictx->self_impostors_data->impostors.erase(id->self_impostor);
    }

    delete base;
}

BranchClusteringData *ImpostorClusteringHelper::convert_branch(Block &settings, Branch *base, ClusteringContext *ictx, 
                                                               BaseBranchClusteringData &data)
{
    isimParams.load(&settings);
    BranchClusteringDataImpostor *id = new BranchClusteringDataImpostor();

    if (!ictx)
    {
        logerr("ImpostorClusteringHelper given clustering context with wrong type");
        return nullptr;
    }
    if (!ictx->self_impostors_data)
    {
        ictx->self_impostors_data = new ImpostorsData();
        int mips = log2(isimParams.impostor_texture_size) + 2;
        TextureAtlas a = TextureAtlas(8*isimParams.impostor_texture_size,8*isimParams.impostor_texture_size, 2, mips);

        ictx->self_impostors_data->atlas = a;
        ictx->self_impostors_data->atlas.set_grid(isimParams.impostor_texture_size, isimParams.impostor_texture_size);
        ictx->self_impostors_data->atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    }
    ImpostorBaker ib;
    ImpostorBaker::ImpostorGenerationParams params;
    params.monochrome = true;
    params.leaf_scale = isimParams.leaf_size_mult;
    params.wood_scale = isimParams.wood_size_mult;
    params.need_top_view = false;
    params.quality = isimParams.impostor_texture_size;
    params.slices_n = isimParams.impostor_similarity_slices;
    params.level_from = isimParams.impostor_metric_level_from;
    params.level_to = isimParams.impostor_metric_level_to;
    params.leaf_opacity = isimParams.leaves_opacity;

    BranchHeap bh;
    LeafHeap lh;
    Branch *tmp_b = bh.new_branch();
    tmp_b->deep_copy(base, bh, &lh);
    glm::mat4 tr = glm::rotate(glm::mat4(1.0f), PI/2, glm::vec3(0,0,1))*glm::inverse(data.transform);
    //when we render impostor we assume that the main axis is y,
    //while after glm::inverse(data.transform) the main axis is x
    tmp_b->transform(tr, 1);

            ClusterData cd;
            cd.base = tmp_b;
            cd.base_pos = 0;
            cd.IDA.type_ids.push_back(base->type_id);
            cd.IDA.tree_ids.push_back(0);
            cd.IDA.centers_par.push_back(tmp_b->center_par);
            cd.IDA.centers_self.push_back(tmp_b->center_self);
            cd.IDA.transforms.push_back(glm::mat4(1.0f));
            cd.ACDA.originals.push_back(nullptr);
            //since we delete full trees right after packing the in clusters originals now mean nothing
            cd.ACDA.ids.push_back(base->self_id);
            cd.ACDA.rotations.push_back(0);
            cd.ACDA.clustering_data.push_back(nullptr);

            ictx->self_impostors_data->impostors.emplace_back();
            BBox bbox = BillboardCloudRaw::get_bbox(tmp_b,glm::vec3(1,0,0),glm::vec3(0,1,0),glm::vec3(0,0,1));
            ib.make_impostor(*tmp_b, (*ictx->types)[base->type_id], ictx->self_impostors_data->impostors.back(), params, 
                             ictx->self_impostors_data->atlas, bbox);
            id->self_impostor = ictx->self_impostors_data->impostors.end();
            id->self_impostor--;
            return id;
}

IntermediateClusteringData *ImpostorClusteringHelper::prepare_intermediate_data(Block &settings, 
                                                                                std::vector<BranchClusteringData *> branches,
                                                                                ClusteringContext *ictx)
{
    isimParams.load(&settings);
    IntermediateClusteringDataDDT *data = new IntermediateClusteringDataDDT();
    if (!ictx)
    {
        logerr("ImpostorClusteringHelper given clustering context with wrong type");
        return nullptr;
    }
    data->branches = branches;
    data->ddt.create(branches.size());
    data->elements_count = branches.size();
    std::vector<BranchClusteringDataImpostor *> real_branches;
    for (int i = 0; i < branches.size(); i++)
    {
        BranchClusteringDataImpostor *imp_dt = dynamic_cast<BranchClusteringDataImpostor *>(branches[i]);
        if (!imp_dt)
        {
            logerr("ImpostorClusteringHelper error: wrong type of BranchClusteringData");
            return data;
        }
        real_branches.push_back(imp_dt);
    }
    ictx->self_impostors_raw_atlas = new TextureAtlasRawData(ictx->self_impostors_data->atlas);
    /*
    if (current_clustering_step == ClusteringStep::BRANCHES)
    {
        textureManager.save_png_raw(ictx->self_impostors_raw_atlas->get_raw_data(),
                                    ictx->self_impostors_raw_atlas->get_w(),
                                    ictx->self_impostors_raw_atlas->get_h(),
                                    4,
                                    "raw png");
    }
    */
    ProgressBar pb = ProgressBar("Preparing DDT with imposter metric", SQR(real_branches.size()),"iterations", true);
    int cnt = 0;
    for (int i = 0; i < real_branches.size(); i++)
    {
        for (int j = 0; j < real_branches.size(); j++)
        {
            Answer a;
            DistData d;
            auto p = data->ddt.get(i,j);
            a = p.first;
            d = p.second;
            if (i == j)
            {
                a.exact = true;
                a.from = 0;
                a.to = 0;
                d.dist = 0;
                d.rotation = 0;
            }
            else if (j < i)
            {
                p = data->ddt.get(j,i);
                a = p.first;
                d = p.second;
            }
            else
            {
                a = dist_impostor(*(real_branches[i]),*(real_branches[j]),ictx,0,1,&d);
            }
            data->ddt.set(i,j,a,d);
            cnt++;
            if (cnt % 1000 == 0)
                pb.iter(cnt);
        }
    }
    /*for (int i = 0; i < real_branches.size(); i++)
    {
        for (int j = 0; j < real_branches.size(); j++)
        {
            debug("%.3f ", data->ddt.get(i,j).first.from);
        }
        debugnl();
    }*/
    delete ictx->self_impostors_raw_atlas;
    return data;
}
int ccnt = 0;

float imp_dist(int w, int h, Billboard &b1, Billboard &b2, TextureAtlasRawData *raw)
{

    float diff = 0;
    float sum = 0;
    unsigned char *data = nullptr;
    if (ccnt<0)
    {
        data = new unsigned char[4*w*h];
        ccnt++;
    }
    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < h; j++)
        {
            if (data)
            {
                data[4*(w*i + j)] = 0;
                data[4*(w*i + j) + 1] = 0;
                data[4*(w*i + j) + 2] = 0;
                data[4*(w*i + j) + 3] = 255;
            }
            unsigned char a1 = raw->get_pixel_uc(i, j, Channel::A, b1.id);
            unsigned char a2 = raw->get_pixel_uc(i, j, Channel::A, b2.id);
            if (a1 == 0 && a2 == 0)
                continue;
            else 
            {
                unsigned char r1 = raw->get_pixel_uc(i, j, Channel::R, b1.id);
                unsigned char r2 = raw->get_pixel_uc(i, j, Channel::R, b2.id);
                unsigned char g1 = raw->get_pixel_uc(i, j, Channel::G, b1.id);
                unsigned char g2 = raw->get_pixel_uc(i, j, Channel::G, b2.id);
                
                float f = (abs(r1-r2) + abs(g1-g2))/255.0;

                diff += f;
                sum += MAX(1,f);

                if (data)
                {
                    data[4*(w*i + j)] = abs(r1-r2);
                    data[4*(w*i + j) + 1] = abs(g1-g2);
                    data[4*(w*i + j) + 2] = 0;
                }
            }
        }
    }
    if (data)
    {
        textureManager.save_png_raw(data, w, h, 4, "debug_"+std::to_string(ccnt));
        delete data;
    }
    return (double)diff/(sum+1);
}
Answer ImpostorClusteringHelper::dist_impostor(BranchClusteringDataImpostor &bwd1, BranchClusteringDataImpostor &bwd2, 
                                               ClusteringContext *current_data, float min, float max, DistData *data)
{
    if (!current_data->self_impostors_raw_atlas || !current_data->self_impostors_data)
        return Answer(true,1000,1000);
    glm::ivec4 sizes = current_data->self_impostors_data->atlas.get_sizes();
    int w = sizes.x/sizes.z;
    int h = sizes.y/sizes.w; 
    
    if (bwd1.self_impostor->slices.empty() || bwd1.self_impostor->slices.size() != bwd2.self_impostor->slices.size())
        return Answer(true,1000,1000);
    float min_av_dist = 1;
    int best_rot = 0;
    int sz = bwd1.self_impostor->slices.size();
    bool simple_comp = isimParams.impostor_metric_level_to > 500;
    if (simple_comp)
    {
        float best_d = 1;
        int best_r = 0;
        for (int r=0;r<sz;r++)
        {
            float dts = imp_dist(w,h,bwd1.self_impostor->slices[0],bwd2.self_impostor->slices[r],
                              current_data->self_impostors_raw_atlas);
            if (dts < best_d)
            {
                best_d = dts;
                best_r = r;
            }
        }
        float av_dst = 0;
        for (int i=0;i<sz;i++)
        {
            av_dst += imp_dist(w,h,bwd1.self_impostor->slices[i],bwd2.self_impostor->slices[(i + best_r)%sz],
                              current_data->self_impostors_raw_atlas);
        }
        av_dst /= sz;
        min_av_dist = av_dst;
    }
    else
    {
        for (int r=0;r<sz;r++)
        {
            float av_dst = 0;
            for (int i=0;i<sz;i++)
            {
                av_dst += imp_dist(w,h,bwd1.self_impostor->slices[i],bwd2.self_impostor->slices[(i + r)%sz],
                                current_data->self_impostors_raw_atlas);
            }
            av_dst /= sz;
            if (av_dst < min_av_dist)
            {
                min_av_dist = av_dst;
                best_rot = r;
            }
        }
    }

    data->rotation = (2*PI*best_rot)/sz;
    #define SZ_DIFF(a,b) pow(MAX(1, MAX(a,b)/MIN(a,b) - isimParams.size_diff_tolerance), isimParams.size_diff_factor)
    glm::vec3 &s1 = bwd1.sizes;
    glm::vec3 &s2 = bwd2.sizes;
    float dist_discriminator = SZ_DIFF(s1.x, s2.x) *
                               SZ_DIFF(sqrt(SQR(s1.y) + SQR(s1.z)), sqrt(SQR(s2.y) + SQR(s2.z)));
     min_av_dist += dist_discriminator - 1;
    return Answer(true,min_av_dist,min_av_dist);
}