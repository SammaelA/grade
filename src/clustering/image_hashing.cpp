#include "hasing.h"
#include "impostor_similarity_params.h"
#include "../impostor.h"
#include "../texture_manager.h"

BranchClusteringData *HashBasedClusteringHelper::convert_branch_impostor(Block &settings, Branch *base, ClusteringContext *ctx,
                                                                         BaseBranchClusteringData &data)
{
    ImpostorSimilarityParams isimParams;
    isimParams.load(&settings);
    BranchHash *id = new BranchHash();

    ImpostorBaker ib;
    ImpostorBaker::ImpostorGenerationParams params;
    params.fixed_colors = true;
    params.leaf_size_mult = isimParams.leaf_size_mult;
    params.wood_size_mult = isimParams.wood_size_mult;
    params.need_top_view = false;
    params.quality = Quality::LOW_AS_F;
    params.slices_n = isimParams.impostor_similarity_slices;
    params.level_from = isimParams.impostor_metric_level_from;
    params.level_to = isimParams.impostor_metric_level_to;

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
    cd.ACDA.originals.push_back(nullptr);
    cd.ACDA.ids.push_back(base->self_id);
    cd.ACDA.rotations.push_back(0);
    cd.ACDA.clustering_data.push_back(nullptr);

    std::list<Impostor>::iterator imp_iter;
    ImpostorsData impData;
    TextureAtlas atl = TextureAtlas(isimParams.impostor_similarity_slices*Quality::LOW_AS_F, Quality::LOW_AS_F, 1);
    impData.atlas = atl;
    impData.atlas.set_grid(Quality::LOW_AS_F, Quality::LOW_AS_F);
    impData.atlas.set_clear_color(glm::vec4(0, 0, 0, 0));

    ib.prepare(params, 1, cd, *(ctx->types), &impData, imp_iter);

    TextureAtlasRawData rawAtlas = TextureAtlasRawData(impData.atlas);
    int hash_size = get_default_block().get_int("impostor_hash_size", 8);
    hash_size = settings.get_int("impostor_hash_size", hash_size);
    
    bool separate_colors = get_default_block().get_bool("separate_colors", true);
    separate_colors = settings.get_bool("separate_colors", separate_colors);
    
    int step = (int)(Quality::LOW_AS_F)/hash_size;
    if (step*hash_size != (int)(Quality::LOW_AS_F))
    {
        debug("warning: impostor_hash_size = %d is not a divider of impostor size = %d", hash_size, (int)(Quality::LOW_AS_F));
    }
    for (int hash_n = 0;hash_n < isimParams.impostor_similarity_slices;hash_n++)
    {
        auto &bill = imp_iter->slices[hash_n];
        id->hashes.emplace_back();
        if (separate_colors)
        {
            id->hashes.back().data = std::vector<float>(2*hash_size*hash_size, 0);
            int cnt = 0;
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
                    id->hashes.back().data[cnt] = res_r/(255*step*step);
                    id->hashes.back().data[cnt+1] = res_g/(255*step*step);
                    //debug("%d %d", (int)(255*id->hashes.back().data[cnt]), (int)(255*id->hashes.back().data[cnt+1]));
                    cnt+=2;
                }
            }
        }
        else
        {
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
                    id->hashes.back().data[cnt] = res/(2*255*step*step);
                    cnt++;
                }
            }            
        }
    }

    return id;
}
