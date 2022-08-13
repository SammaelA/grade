#include "hasing.h"
#include "impostor_similarity_params.h"
#include "graphics_utils/impostor.h"
#include "graphics_utils/texture_manager.h"

void HashBasedClusteringHelper::create_impostor_temp(Block &settings, Branch *base, ClusteringContext *ctx, 
                          BaseBranchClusteringData &data, ImpostorSimilarityParams &isimParams, 
                          TextureAtlas &temp_atlas, Impostor &imp)
{

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
    glm::mat4 tr = glm::rotate(glm::mat4(1.0f), PI / 2, glm::vec3(0, 0, 1)) * glm::inverse(data.transform);
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
    cd.ACDA.clustering_data.push_back(nullptr);
    
    BBox bbox = BillboardCloudRaw::get_bbox(tmp_b,glm::vec3(1,0,0),glm::vec3(0,1,0),glm::vec3(0,0,1));
    ib.make_impostor(*tmp_b, ctx->types[base->type_id], imp, params, temp_atlas, bbox);
}

BranchClusteringData *HashBasedClusteringHelper::convert_branch_impostor(Block &settings, Branch *base, ClusteringContext *ctx,
                                                                         BaseBranchClusteringData &data)
{
    ImpostorSimilarityParams isimParams;
    isimParams.load(&settings);
    BranchHash *id = new BranchHash();
    Impostor imp;
    
    TextureAtlas atl = TextureAtlas(isimParams.impostor_similarity_slices*isimParams.impostor_texture_size, isimParams.impostor_texture_size, 1);
    atl.set_grid(isimParams.impostor_texture_size, isimParams.impostor_texture_size);
    atl.set_clear_color(glm::vec4(0, 0, 0, 0));
    
    create_impostor_temp(settings, base, ctx, data, isimParams, atl, imp);
    TextureAtlasRawData rawAtlas = TextureAtlasRawData(atl);

    int hash_size = get_default_block().get_int("impostor_hash_size", 8);
    hash_size = settings.get_int("impostor_hash_size", hash_size);
    hash_size = MIN(hash_size, (int)(isimParams.impostor_texture_size));

    bool separate_colors = get_default_block().get_bool("separate_colors", true);
    separate_colors = settings.get_bool("separate_colors", separate_colors);
    
    bool relative_to_average = get_default_block().get_bool("relative_to_average", false);
    relative_to_average = settings.get_bool("relative_to_average", relative_to_average);

    int step = (int)(isimParams.impostor_texture_size)/hash_size;
    if (step*hash_size != (int)(isimParams.impostor_texture_size))
    {
        debug("warning: impostor_hash_size = %d is not a divider of impostor size = %d", hash_size, (int)(isimParams.impostor_texture_size));
    }
    for (int hash_n = 0;hash_n < isimParams.impostor_similarity_slices;hash_n++)
    {
        auto &bill = imp.slices[hash_n];
        id->hashes.emplace_back();
        auto &H = id->hashes.back().data;
        if (separate_colors)
        {
            float av_r = 0;
            float av_g = 0;
            H = std::vector<float>(2*hash_size*hash_size, 0);
            int cnt = 0;
            logerr("hash size %d %d", hash_size, step);
            for (int i=0;i<hash_size;i++)
            {
                for (int j=0;j<hash_size;j++)
                {
                    float res_r = 0;
                    float res_g = 0;
                    for (int x = i*step;x<(i+1)*step;x++)
                    {
                        for (int y = i*step;y<(i+1)*step;y++)
                        {
                            res_r += rawAtlas.get_pixel_uc(x,y,Channel::R,bill.id);
                            res_g += rawAtlas.get_pixel_uc(x,y,Channel::G,bill.id);
                        }
                    }
                    H[cnt] = res_r/(255*step*step);
                    H[cnt+1] = res_g/(255*step*step);
                    av_r += H[cnt];
                    av_g += H[cnt+1];

                    cnt+=2;
                }
            }
            if (relative_to_average)
            {
                av_r = av_r/(hash_size*hash_size);
                av_g = av_g/(hash_size*hash_size);

                for (int i=0;i<2*hash_size*hash_size;i+=2)
                {
                    H[i] = (H[i] > av_r) ? 1 : 0;
                    H[i + 1] = (H[i] > av_g) ? 1 : 0;
                }
            }
        }
        else
        {
            float av = 0;
            id->hashes.back().data = std::vector<float>(hash_size*hash_size, 0);
            int cnt = 0;
            for (int i=0;i<hash_size;i++)
            {
                for (int j=0;j<hash_size;j++)
                {
                    float res = 0;
                    for (int x = i*step;x<(i+1)*step;x++)
                    {
                        for (int y = i*step;y<(i+1)*step;y++)
                        {
                            res += rawAtlas.get_pixel_uc(x,y,Channel::R,bill.id);
                            res += rawAtlas.get_pixel_uc(x,y,Channel::G,bill.id);
                        }
                    }
                    H[cnt] = res/(2*255*step*step);
                    av += H[cnt];
                    cnt++;
                }
            }
            if (relative_to_average)
            {
                av = av/(hash_size*hash_size);

                for (int i=0;i<hash_size*hash_size;i++)
                {
                    H[i] = (H[i] > av) ? 1 : 0;
                }
            }          
        }
        for (auto &f : H)
        {
            debug("%f ",f);
        }
        debugnl();
    }

    return id;
}

void dct(float *values, float *dcts, int N, int M)
{
    for (int u = 0;u<N;u++)
    {
        for (int v=0;v<M;v++)
        {
            dcts[u*M + v] = 0;
            for (int i = 0;i<N;i++)
            {
                for (int j=0;j<M;j++)
                {
                    dcts[u*M + v] += values[i*M + j]*cos(M_PI/((float)N)*(i+1./2.)*u)*cos(M_PI/((float)M)*(j+1./2.)*v);          
                }
            }
        }
    }
}

BranchClusteringData *HashBasedClusteringHelper::convert_branch_impostor_dct(Block &settings, Branch *base, 
                                                                             ClusteringContext *ctx, 
                                                                             BaseBranchClusteringData &data)
{
    ImpostorSimilarityParams isimParams;
    isimParams.load(&settings);
    BranchHash *id = new BranchHash();
    Impostor imp;

    TextureAtlas atl = TextureAtlas(isimParams.impostor_similarity_slices*isimParams.impostor_texture_size, isimParams.impostor_texture_size, 1);
    atl.set_grid(isimParams.impostor_texture_size, isimParams.impostor_texture_size);
    atl.set_clear_color(glm::vec4(0, 0, 0, 0));

    create_impostor_temp(settings, base, ctx, data, isimParams, atl, imp);
    TextureAtlasRawData rawAtlas = TextureAtlasRawData(atl);

    int dct_block_size = get_default_block().get_int("dct_block_size", 32);
    dct_block_size = settings.get_int("dct_block_size", dct_block_size);
    dct_block_size = MIN(dct_block_size, (int)(isimParams.impostor_texture_size));

    int dct_use_size = get_default_block().get_int("dct_use_size", 8);
    dct_use_size = settings.get_int("dct_use_size", dct_use_size);
    dct_use_size = MIN(dct_use_size, dct_block_size);

    bool separate_colors = get_default_block().get_bool("separate_colors", true);
    separate_colors = settings.get_bool("separate_colors", separate_colors);

    int step = (int)(isimParams.impostor_texture_size)/dct_block_size;
    for (int hash_n = 0;hash_n < isimParams.impostor_similarity_slices;hash_n++)
    {
        auto &bill = imp.slices[hash_n];
        id->hashes.emplace_back();
        auto &H = id->hashes.back().data;
        H = std::vector<float>((1 + (int)(separate_colors))*dct_use_size*dct_use_size -1, 0);
        for (int channel=0;channel<1 + (int)(separate_colors);channel++)
        {
            float *values = new float[dct_block_size*dct_block_size];
            float *dcts = new float[dct_block_size*dct_block_size];

            for (int i = 0; i < dct_block_size; i++)
            {
                for (int j = 0; j < dct_block_size; j++)
                {
                    float res = 0;
                    for (int x = i * step; x < (i + 1) * step; x++)
                    {
                        for (int y = i * step; y < (i + 1) * step; y++)
                        {
                            res += rawAtlas.get_pixel_uc(x, y, Channel(channel), bill.id);
                            if (!separate_colors)
                                res += rawAtlas.get_pixel_uc(x, y, Channel::G, bill.id);
                        }
                    }
                    values[i*dct_block_size + j] = res/(step*step*(2-(int)separate_colors));
                }
            }

            dct(values, dcts, dct_block_size, dct_block_size);
            int cnt = 0;
            for (int i = 0; i < dct_block_size; i++)
            {
                for (int j = 0; j < dct_block_size; j++)
                {
                    if (i+j)
                    {
                        H[cnt] = dcts[i*dct_block_size + j];
                    }
                }
            }
        }
    }

    return id;
}