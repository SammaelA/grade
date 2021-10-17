#include "GPU_impostor_metric.h"
#include "dist_data_table.h"
#include "../billboard_cloud.h"
#include "../texture_manager.h"
#include "../tinyEngine/TinyEngine.h"

IntermediateClusteringData *GPUImpostorClusteringHelper::prepare_intermediate_data(Block &settings, 
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

    int tasks_cnt_per_impostor = MIN(4,real_branches[0]->self_impostor->slices.size());
    int bil_step = (real_branches[0]->self_impostor->slices.size())/tasks_cnt_per_impostor;
    int sz = real_branches.size();
    int tasks_cnt = tasks_cnt_per_impostor*tasks_cnt_per_impostor*sz*sz;
    Task *tasks = new Task[tasks_cnt];
    float *results = new float[tasks_cnt];
    int p = 0;
    auto &atlas = ictx->self_impostors_data->atlas;
    ProgressBar pb_prepare = ProgressBar("GPU impostor clustering prepare tasks", SQR(sz)*SQR(tasks_cnt_per_impostor), "tasks");

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
                        Billboard &b1 = real_branches[i]->self_impostor->slices[bil_step*k];
                        Billboard &b2 = real_branches[j]->self_impostor->slices[bil_step*l];
                        atlas.pixel_offsets(b1.id,tasks[p].from);
                        atlas.pixel_offsets(b2.id,tasks[p].to);
                    }
                    else
                    {
                        tasks[p].from = glm::uvec4(0,0,0,0);
                        tasks[p].to = glm::uvec4(0,0,0,0);
                    }
                    p++;
                    if (p % 100000 == 0)
                        pb_prepare.iter(p);
                }
            }
        }
    }
    pb_prepare.finish();
    ProgressBar pb_buffers = ProgressBar("GPU impostor clustering prepare buffers",3,"buffers");
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

    pb_buffers.finish();
        
        Shader impMetric({"impostor_dist.comp"},{});
            impMetric.use();
            impMetric.uniform("tex_sz_x", (float)sizes.x);
            impMetric.uniform("tex_sz_y", (float)sizes.y);
            impMetric.uniform("tasks_count",tasks_cnt);
            impMetric.uniform("impostor_x",sizes.x/sizes.z);
            impMetric.uniform("impostor_y",sizes.y/sizes.w);
            impMetric.texture("atlas", atlas.tex(0));
            static int max_dispatches = settings.get_int("max_dispatches",16);
            static int threads = 8;
            int step = (max_dispatches*threads);
            int iters = ceil((float)tasks_cnt/step);

            int start_id = 0;
            ProgressBar pb_dispatch = ProgressBar("GPU impostor clustering", iters, "iterations");
            for (int i=0;i<iters;i++)
            {
                impMetric.uniform("start_id",start_id);
                glDispatchCompute(max_dispatches, max_dispatches, 1);

                start_id += step;
                if (i % 100 == 0)
                {
                    pb_dispatch.iter(i);
                    SDL_GL_SwapWindow(Tiny::view.gWindow);
                }
            }
            pb_dispatch.finish();
            glMemoryBarrier(GL_ALL_BARRIER_BITS);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
            GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
            memcpy(results,ptr,sizeof(float)*tasks_cnt);

            glDeleteBuffers(1, &tasks_buf);
            glDeleteBuffers(1, &results_buf);
            glDeleteBuffers(1, &raw_data_buf);

    int cnt = 0;
    ProgressBar pb_finalize = ProgressBar("GPU impostor clustering finalize", SQR(real_branches.size()), "branches");
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
                #define SZ_DIFF(a,b) pow(MAX(1, MAX(a,b)/MIN(a,b) - isimParams.size_diff_tolerance), isimParams.size_diff_factor)
                glm::vec3 &s1 = real_branches[i]->min_bbox.sizes;
                glm::vec3 &s2 = real_branches[j]->min_bbox.sizes;
                float dist_discriminator = SZ_DIFF(s1.x, s2.x) *
                                        SZ_DIFF(sqrt(SQR(s1.y) + SQR(s1.z)), sqrt(SQR(s2.y) + SQR(s2.z)));
                min_av_dist *= dist_discriminator;
                a = Answer(true,min_av_dist,min_av_dist);
                d = (2*PI*best_rot)/tasks_cnt_per_impostor;
            }
            data->ddt.set(i,j,a,d);
            if (cnt % 10000 == 0)
            {
                pb_finalize.iter(cnt);
            }
            cnt++;
        }
    }
    pb_finalize.finish();
    
    delete tasks;
    delete results;
    return data;
}