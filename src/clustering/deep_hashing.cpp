#include "deep_hashing.h"
#include "impostor_similarity_params.h"
#include "../impostor.h"
#include "../texture_manager.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

BranchClusteringData *DeepHashBasedClusteringHelper::convert_branch_deep_hash(Block &settings, Branch *base, 
                                                      ClusteringContext *ctx, 
                                                      BaseBranchClusteringData &data)
{
    BranchHash *branchHash = new BranchHash();
    res_branches.push_back(branchHash);
    src_branches.push_back(base);
    base_branches.push_back(data);

    return branchHash;
}

void DeepHashBasedClusteringHelper::branch_conversion_flush(Block &settings, ClusteringContext *ctx)
{
    ImpostorSimilarityParams isimParams;
    isimParams.load(&settings);
    int b_count = res_branches.size();
        
    ImpostorsData impData;
    TextureAtlas atl = TextureAtlas(isimParams.impostor_similarity_slices*Quality::LOW_AS_F, Quality::LOW_AS_F, 1);
    impData.atlas = atl;
    impData.atlas.set_grid(Quality::LOW_AS_F, Quality::LOW_AS_F);
    impData.atlas.set_clear_color(glm::vec4(0, 0, 0, 0));
    std::vector<std::list<Impostor>::iterator> imp_iters;
    for (int br_n = 0; br_n<b_count;br_n++)
    {
        BranchHash *id = res_branches[br_n];
        Branch *base = src_branches[br_n];
        BaseBranchClusteringData &data = base_branches[br_n];
        std::list<Impostor>::iterator imp_iter;
        
        create_impostor_temp(settings, base, ctx, data, isimParams, imp_iter, impData);

        imp_iters.push_back(imp_iter);
    }
    
    TextureAtlasRawData rawAtlas = TextureAtlasRawData(impData.atlas);
    //textureManager.save_bmp(impData.atlas.tex(0),"flush_atlas");

        int hash_size = get_default_block().get_int("impostor_hash_size", 8);
        hash_size = settings.get_int("impostor_hash_size", hash_size);
        
        bool separate_colors = get_default_block().get_bool("separate_colors", true);
        separate_colors = settings.get_bool("separate_colors", separate_colors);
        
        bool relative_to_average = get_default_block().get_bool("relative_to_average", false);
        relative_to_average = settings.get_bool("relative_to_average", relative_to_average);

        int step = (int)(Quality::LOW_AS_F)/hash_size;
        if (step*hash_size != (int)(Quality::LOW_AS_F))
        {
            debug("warning: impostor_hash_size = %d is not a divider of impostor size = %d", hash_size, (int)(Quality::LOW_AS_F));
        }
    std::string dir_path = "./tmp";
    std::string cluster_labels = " 1";
    for (int k =0;k<32;k++)
    {
        cluster_labels = cluster_labels + " 0";
    }
    std::string database_file;
    std::string database_name = "database.txt";
    int sl = 0, ww = 0, hh = 0;
    unsigned char *sl_data = safe_new<unsigned char>(rawAtlas.get_slice_size(0), "sl_data");

    for (int br_n = 0; br_n<b_count;br_n++)
    {
        std::list<Impostor>::iterator imp_iter = imp_iters[br_n];
        BranchHash *id = res_branches[br_n];

        for (int hash_n = 0;hash_n < isimParams.impostor_similarity_slices;hash_n++)
        {
            std::string file_path = dir_path + "/" + std::to_string(sl)+".bmp";
            auto &bill = imp_iter->slices[hash_n];
            rawAtlas.get_slice(bill.id, sl_data, &ww, &hh);
            textureManager.save_bmp_raw_directly(sl_data, ww, hh, 4, file_path);
            std::string record = "./" + std::to_string(sl)+".bmp" + cluster_labels + "\n";
            database_file += record;
            sl++;

            id->hashes.emplace_back();
        }

    }
    
    std::ofstream database_ofs;
    database_ofs.open(dir_path+"/prepare.txt");
    database_ofs << database_file;
    database_ofs.close();
    

            double *data = nullptr;
            int hash_len = 0;
            int hash_count = 0;
            int expected_hash_len = 64;
            std::string script_name = settings.get_string("script", "get_hashes");
            std::string script_args = script_name;
            std::string model_path = "./models/tm_1_dch.npy";
            model_path = settings.get_string("model",model_path);
            script_args = script_args + " --model-weights = " + model_path;
            
            PythonHelper ph;

            ph.init("./scripts/deep_hashing");

            ph.run_script(script_name, script_args);
            bool res = ph.get_numpy_2d_array_double("arr",&hash_count,  &hash_len, &data);
            ph.finish_script();
            //ph.run_script(script_name, script_args);
            //bool res = ph.get_numpy_2d_array_double("arr",&hash_count,  &hash_len, &data);
            //ph.finish_script();
            //ph.finish();
            if (res)
            {
                if (hash_len != expected_hash_len || hash_count < b_count*isimParams.impostor_similarity_slices)
                {
                    logerr("deep hashing error data got from python script %d %d %d %d",
                           hash_len, expected_hash_len, hash_count, b_count*isimParams.impostor_similarity_slices);
                }
                else
                {
                    int cnt = 0;
                    for (int br_n = 0; br_n<b_count;br_n++)
                    {
                        BranchHash *id = res_branches[br_n];

                        for (int hash_n = 0;hash_n < isimParams.impostor_similarity_slices;hash_n++)
                        {
                            id->hashes[hash_n].data = std::vector<float>(hash_len, 0);
                            for (int i=0;i<hash_len;i++)
                            {
                                id->hashes[hash_n].data[i] = data[hash_len*cnt + i];
                            }
                            cnt++;
                        }
                    }

                }
                delete data;
            }
            else
            {
                logerr("deep hashing failed");
            }

    safe_delete<unsigned char>(sl_data, "sl_data");
    rawAtlas.clear();
}