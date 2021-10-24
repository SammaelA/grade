#include "experimental_hashing.h"
#include "impostor_similarity_params.h"
#include "../texture_manager.h"

BranchClusteringData *ImpostorHashClusteringHelper2::convert_branch(Block &settings, Branch *base,
                                                                    ClusteringContext *ctx,
                                                                    BaseBranchClusteringData &data)
{
    BranchHash *branchHash = new BranchHash();
    res_branches.push_back(branchHash);
    src_branches.push_back(base);
    base_branches.push_back(data);

    return branchHash;
}

IntermediateClusteringData *ImpostorHashClusteringHelper2::prepare_intermediate_data(Block &settings, std::vector<BranchClusteringData *> branches,
                                                                                     ClusteringContext *ctx)
{
    return prepare_intermediate_data_ddt(settings, branches, ctx);
}
//BranchClusteringData that is returned after convert_branch could be not fully initialized.
//call this function before using it.
void ImpostorHashClusteringHelper2::branch_conversion_flush(Block &settings, ClusteringContext *ctx)
{
    ImpostorSimilarityParams isimParams;
    isimParams.load(&settings);
    int b_count = res_branches.size();

    TextureAtlas atl = TextureAtlas(isimParams.impostor_similarity_slices * isimParams.impostor_texture_size,
                                    isimParams.impostor_texture_size, 1);
    atl.set_grid(isimParams.impostor_texture_size, isimParams.impostor_texture_size);
    atl.set_clear_color(glm::vec4(0, 0, 0, 0));
    std::vector<Impostor> imp_iters;

    for (int br_n = 0; br_n < b_count; br_n++)
    {
        BranchHash *id = res_branches[br_n];
        Branch *base = src_branches[br_n];
        BaseBranchClusteringData &data = base_branches[br_n];
        imp_iters.emplace_back();

        create_impostor_temp(settings, base, ctx, data, isimParams, atl, imp_iters.back());
    }

    TextureAtlasRawData rawAtlas = TextureAtlasRawData(atl);

        int hash_size = get_default_block().get_int("impostor_hash_size", 8);
    hash_size = settings.get_int("impostor_hash_size", hash_size);
    hash_size = MIN(hash_size, (int)(isimParams.impostor_texture_size));

    int step = (int)(isimParams.impostor_texture_size)/hash_size;
    if (step*hash_size != (int)(isimParams.impostor_texture_size))
    {
        debug("warning: impostor_hash_size = %d is not a divider of impostor size = %d", hash_size, (int)(isimParams.impostor_texture_size));
    }
    glm::vec2 averige_sizes = glm::vec2(0,0);
    for (int i=0;i<b_count;i++)
    {
        averige_sizes.x += res_branches[i]->sizes.x;
        averige_sizes.y += sqrt(SQR(res_branches[i]->sizes.y)*SQR(res_branches[i]->sizes.z));
    }
    averige_sizes = (1.0f/b_count)*averige_sizes;
    for (int i=0;i<b_count;i++)
    {
        BranchHash *id = res_branches[i];
        id->hashes.emplace_back();
        auto &H = id->hashes.back().data;
        int cnt = 0;
        H = std::vector<float>(2*isimParams.impostor_similarity_slices*hash_size*hash_size + 2, 0);
        float sum = 0;
        for (int hash_n = 0; hash_n < isimParams.impostor_similarity_slices; hash_n++)
        {
            auto &bill = imp_iters[i].slices[hash_n];
            int cnt = 0;

            for (int i = 0; i < hash_size; i++)
            {
                for (int j = 0; j < hash_size; j++)
                {
                    float res_r = 0;
                    float res_g = 0;
                    for (int x = i * step; x < (i + 1) * step; x++)
                    {
                        for (int y = i * step; y < (i + 1) * step; y++)
                        {
                            res_r += rawAtlas.get_pixel_uc(x, y, Channel::R, bill.id);
                            res_g += rawAtlas.get_pixel_uc(x, y, Channel::G, bill.id);
                        }
                    }
                    H[cnt] = res_r / (255 * step * step);
                    sum += H[cnt];
                    H[cnt + 1] = res_g / (255 * step * step);
                    sum += H[cnt+1];
                    cnt += 2;
                }
            }
        }
        float diff_x = (id->sizes.x > averige_sizes.x) ? (id->sizes.x/averige_sizes.x - 1) : (1 - averige_sizes.x/id->sizes.x);
        float av_sz = sqrt(SQR(id->sizes.y)*SQR(id->sizes.z));
        float diff_y = (av_sz > averige_sizes.y) ? (av_sz/averige_sizes.y - 1) : (1 - averige_sizes.y/av_sz);

        H[cnt] = sum*diff_x;
        H[cnt + 1] = sum*diff_y;
        logerr("%f %f",H[cnt],H[cnt+1]);
    }
    
    rawAtlas.clear();
}