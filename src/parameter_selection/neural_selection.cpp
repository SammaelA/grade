#include "neural_selection.h"
#include "tree_generators/all_generators.h"
#include "generation/grove_packer.h"
#include "graphics_utils/texture_manager.h"
#include "parameter_selection_utils.h"
#include "save_utils/csv.h"
#include <vector>
#include <thread>

void NeuralEstimator::prepare_dataset(ParameterList &param_list, GroveGenerationData &tree_ggd, 
                                      std::string save_path, int impostor_size, int impostor_slices,
                                      int batch_size, int batches_cnt, bool save_metainfo,
                                      int joints_limit)
{
    AbstractTreeGenerator::set_joints_limit(joints_limit);
    Heightmap hmap = Heightmap(glm::vec3(0,0,0), glm::vec2(100,100), 10);
    std::vector<float> par_list_normalized;
    param_list.to_simple_list(par_list_normalized, true, true);
    if (par_list_normalized.empty())
    {
        logerr("empty parameter list for prepare_dataset");
        return;
    }
    std::vector<ParameterList> params(batch_size, param_list);
    std::vector<std::vector<float>> params_normalzied(batch_size);
    tree_ggd.impostor_generation_params.slices_n = impostor_slices;
    tree_ggd.impostor_generation_params.quality = impostor_size;
    
    int total_images_count = 0;
    int params_cnt = par_list_normalized.size();
    logerr("generator has %d parameters", params_cnt);
    std::vector<std::string> columns;
    columns.push_back("image_dir");
    for (int i=0;i<params_cnt;i++)
        columns.push_back("par_"+std::to_string(i));
    CSVData dataset(columns);
    unsigned char *sl_data = new unsigned char[4*impostor_size*impostor_size];

    for (int batch = 0; batch < batches_cnt; batch++)
    {
        for (int i=0;i<batch_size;i++)
        {
            for (auto &p : par_list_normalized)
                p = urand(0,1);
            params_normalzied[i] = par_list_normalized;
            params[i].from_simple_list(par_list_normalized, true, true);
            
        }

        textureManager.set_textures_tag(1);
        GrovePacker packer;
        Tree *trees = new Tree[params.size()];
        GrovePacked tmp_g;
        for (int i = tree_ggd.types.size(); i < params.size(); i++)
        {
            tree_ggd.types.push_back(tree_ggd.types[0]);
        }
        tree_ggd.trees_count = params.size();
        for (int i = 0; i < params.size(); i++)
        {
            tree_ggd.types[i].get_params()->read_parameter_list(params[i]);
        }

        int num_threads = MIN(16, params.size());
        int step = ceil(params.size() / (float)num_threads);
        LightVoxelsCube **thr_voxels = new LightVoxelsCube *[num_threads];
        bool voxels_needed = false;
        glm::vec3 max_size = glm::vec3(100,100,100);
        float min_base_size = 1000;
        for (auto &t : tree_ggd.types)
        {
            auto *gen = get_generator(t.generator_name);
            if (gen->use_voxels_for_generation())
                voxels_needed = true;
            max_size = max(max_size, t.get_params()->get_tree_max_size());
            min_base_size = MIN(min_base_size, t.get_params()->get_scale_factor());
            delete gen;
        }
        for (int i = 0; i < num_threads; i++)
        {
            if (voxels_needed)
                thr_voxels[i] = new LightVoxelsCube(glm::vec3(0,0,0), max_size, min_base_size, 1, 1, 2);
            else
                thr_voxels[i] = new LightVoxelsCube(glm::vec3(0,0,0),glm::ivec3(1,1,1),1,1,2);
        }
        std::vector<std::thread> threads;
        for (int i = 0; i < num_threads; i++)
        {
            int start_n = step * i;
            int stop_n = MIN(step * (i + 1), params.size());
            threads.push_back(std::thread(&ps_utils::gen_tree_task, start_n, stop_n, thr_voxels[i], &(tree_ggd.types), trees));
        }
        for (auto &t : threads)
            t.join();
        
        std::vector<int> valid_ids;
        int valid_cnt = 0;
        for (int i=0;i<params.size();i++)
        {
            if (GrovePacker::is_valid_tree(trees[i]))
            {
                valid_cnt++;
            }
            else
            {
                valid_ids.push_back(i);
            }
        }
        logerr("valid %d/%d", valid_cnt, batch_size);
        packer.add_trees_to_grove(tree_ggd, tmp_g, trees, &hmap, false);
        TextureAtlasRawData raw_atlas(tmp_g.impostors[1].atlas);
        int imp_n = 0;
        for (auto &imp : tmp_g.impostors[1].impostors)
        {
            for (auto &slice : imp.slices)
            {
                raw_atlas.get_slice(slice.id, sl_data);
                std::string dir = "images/" + std::to_string(total_images_count) + ".png";
                std::vector<std::string> values;
                values.push_back(dir);
                for (int i=0;i<params_cnt;i++)
                {
                    values.push_back(std::to_string(params_normalzied[imp_n][i]));
                }
                dataset.add_row(values);
                textureManager.save_png_raw_directly(sl_data, impostor_size, impostor_size, 4, save_path + "/" + dir);
                total_images_count++;
            }
            imp_n++;
        }

        delete[] trees;
        for (int i = 0; i < num_threads; i++)
        {
            if (thr_voxels[i])
                delete thr_voxels[i];
        }
        delete[] thr_voxels;
        for (auto &imp : tmp_g.impostors)
        {
            imp.atlas.destroy();
        }
        textureManager.clear_unnamed_with_tag(1);
    }

    CSVSaver saver;
    saver.save_csv_in_file(dataset, save_path + "/dataset.csv");
    delete[] sl_data;
}