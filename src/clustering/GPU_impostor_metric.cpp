#include "GPU_impostor_metric.h"
#include "dist_data_table.h"
#include "../billboard_cloud.h"
#include "../texture_manager.h"
#include "../tinyEngine/TinyEngine.h"

IntermediateClusteringData *GPUImpostorClusteringHelper::prepare_intermediate_data(Block &settings, 
                                                                                   std::vector<BranchClusteringData *> branches,
                                                                                   ClusteringContext *ictx)
{
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
    if (current_clustering_step == ClusteringStep::BRANCHES)
    {
        textureManager.save_bmp_raw(ictx->self_impostors_raw_atlas->get_raw_data(),
                                    ictx->self_impostors_raw_atlas->get_w(),
                                    ictx->self_impostors_raw_atlas->get_h(),
                                    4,
                                    "raw bmp");
    }
    int tasks_cnt_per_impostor = real_branches[0]->self_impostor->slices.size();
    int sz = real_branches.size();
    int tasks_cnt = tasks_cnt_per_impostor*tasks_cnt_per_impostor*sz*sz;
    Task *tasks = new Task[tasks_cnt];
    float *results = new float[tasks_cnt];
    int p = 0;
    auto &atlas = ictx->self_impostors_data->atlas;
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
        {
            for (int k = 0; k<tasks_cnt_per_impostor;k++)
            {
                for (int l = 0; l<tasks_cnt_per_impostor;l++)
                {
                    if (j>i)
                    {
                        Billboard &b1 = real_branches[i]->self_impostor->slices[k];
                        Billboard &b2 = real_branches[j]->self_impostor->slices[l];
                        atlas.pixel_offsets(b1.id,tasks[p].from);
                        atlas.pixel_offsets(b2.id,tasks[p].to);
                    }
                    else
                    {
                        tasks[p].from = glm::uvec4(0,0,0,0);
                        tasks[p].to = glm::uvec4(0,0,0,0);
                    }
                    p++;
                }
            }
        }
    }
        glm::ivec4 sizes = atlas.get_sizes();
        int layers = atlas.layers_count();
        glGenBuffers(1, &tasks_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, tasks_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(Task)*tasks_cnt, tasks, GL_STATIC_DRAW);

        glGenBuffers(1, &results_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, results_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*tasks_cnt, results, GL_STATIC_DRAW);

        glGenBuffers(1, &raw_data_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, raw_data_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(unsigned)*sizes.x*sizes.y*layers, 
                     ictx->self_impostors_raw_atlas->get_raw_data(), GL_STATIC_DRAW);

        //glBindImageTexture(0, atlas.tex(0).texture, 0, true, 0, GL_READ_ONLY, GL_RGBA8);

        
        Shader impMetric({"impostor_dist.comp"},{});
            impMetric.use();
            impMetric.uniform("impostor_x", sizes.x);
            impMetric.uniform("impostor_y", sizes.y);
            impMetric.uniform("tasks_count",tasks_cnt);
            impMetric.uniform("impostor_x",sizes.x/sizes.z);
            impMetric.uniform("impostor_y",sizes.y/sizes.w);
            static int max_dispatches = settings.get_int("max_dispatches",16);
            static int threads = 8;
            int step = (max_dispatches*threads);
            int iters = ceil((float)tasks_cnt/step);

            int start_id = 0;
            for (int i=0;i<iters;i++)
            {
                impMetric.uniform("start_id",start_id);
                glDispatchCompute(max_dispatches, max_dispatches, 1);

                SDL_GL_SwapWindow(Tiny::view.gWindow);
                start_id += step;
            }
        
            glMemoryBarrier(GL_ALL_BARRIER_BITS);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
            GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
            memcpy(results,ptr,sizeof(float)*tasks_cnt);

            glDeleteBuffers(1, &tasks_buf);
            glDeleteBuffers(1, &results_buf);
            glDeleteBuffers(1, &raw_data_buf);

    /*for (int i=0;i<tasks_cnt;i++)
    {
        if (tasks[i].from.w > 0)
            logerr("res[%d] = %f",i, results[i]);
    }*/
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
                int start_pos = SQR(tasks_cnt_per_impostor)*(i*sz + j);
                float min_av_dist = 1;
                int best_rot = 0;
                
                for (int r=0;r<tasks_cnt_per_impostor;r++)
                {
                    float av_dst = 0;
                    for (int t=0;t<tasks_cnt_per_impostor;t++)
                    {
                        int pos = start_pos + t*tasks_cnt_per_impostor + (t + r)%tasks_cnt_per_impostor;
                        av_dst += results[pos];
                    }
                    av_dst /= tasks_cnt_per_impostor;
                    if (av_dst < min_av_dist)
                    {
                        min_av_dist = av_dst;
                        best_rot = r;
                    }
                }
                a = Answer(true,min_av_dist,min_av_dist);
                d = (2*PI*best_rot)/tasks_cnt_per_impostor;
            }
            data->ddt.set(i,j,a,d);
        }
    }
    delete tasks;
    delete results;
    delete ictx->self_impostors_raw_atlas;
    return data;
}