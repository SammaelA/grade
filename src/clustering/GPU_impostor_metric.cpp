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
    float *results = new float[2*SQR(sz)];
    int slices_cnt = sz*tasks_cnt_per_impostor;
    glm::uvec4 *slices_info = new glm::uvec4[slices_cnt];
    int p = 0;
    auto &atlas = ictx->self_impostors_data->atlas;
    
    for (int i = 0; i < sz; i++)
    {
        for (int k = 0; k<tasks_cnt_per_impostor;k++)
        {
            Billboard &b1 = real_branches[i]->self_impostor->slices[bil_step*k];
            atlas.pixel_offsets(b1.id,slices_info[p]);
            p++;
        }
    }

    ProgressBar pb_buffers = ProgressBar("GPU impostor clustering prepare buffers",3,"buffers");
        glm::ivec4 sizes = atlas.get_sizes();
        int layers = atlas.layers_count();

        glGenBuffers(1, &results_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, results_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*2*SQR(sz), results, GL_STATIC_DRAW);

        glGenBuffers(1, &slices_info_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slices_info_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::uvec4)*slices_cnt, slices_info, GL_STATIC_DRAW);

    pb_buffers.finish();
        
        Shader impMetric({"impostor_dist.comp"},{});
            impMetric.use();
            impMetric.uniform("tex_sz_x", (float)sizes.x);
            impMetric.uniform("tex_sz_y", (float)sizes.y);
            impMetric.uniform("slice_count",tasks_cnt_per_impostor);
            impMetric.uniform("branch_count",sz);
            impMetric.uniform("impostor_x",sizes.x/sizes.z);
            impMetric.uniform("impostor_y",sizes.y/sizes.w);
            impMetric.texture("atlas", atlas.tex(0));
            static int max_dispatches = settings.get_int("max_dispatches",16);
            static int threads = 64;
            int step = (max_dispatches*threads);
            int iters_count = 0;
            for (int row = 0;row<sz;row++)
            {
                for (int column = row + 1; column < sz; column+=step)
                {
                    int dispatches_left = ceil(((float)(sz - column))/threads);
                    int dispatches = MIN(dispatches_left, max_dispatches);
                    iters_count += dispatches*threads;
                }
            }
            ProgressBar pb_dispatch = ProgressBar("GPU impostor clustering", iters_count, "iterations");
            p = 0;
            int prev_p = 0;
            int kk = 0;
            for (int row = 0;row<sz;row++)
            {
                impMetric.uniform("row", row);
                for (int column = row + 1; column < sz; column+=step)
                {
                    impMetric.uniform("start_column", column);
                    int dispatches_left = ceil(((float)(sz - column))/threads);
                    int dispatches = MIN(dispatches_left, max_dispatches);
                    glDispatchCompute(dispatches, 1, 1);
                    p += dispatches*threads;
                    if (p - prev_p > 5000)
                    {
                        if (kk % 5 == 0)
                            pb_dispatch.iter(p);
                        SDL_GL_SwapWindow(Tiny::view.gWindow);
                        prev_p = p;
                        kk++;
                    }
                }
            }

            pb_dispatch.finish();
            glMemoryBarrier(GL_ALL_BARRIER_BITS);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
            GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
            memcpy(results,ptr,sizeof(float)*2*SQR(sz));

            glDeleteBuffers(1, &slices_info_buf);
            glDeleteBuffers(1, &results_buf);

    int cnt = 0;
    ProgressBar pb_finalize = ProgressBar("GPU impostor clustering finalize", SQR(sz), "branches");
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
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
                float min_av_dist = results[2*(sz*i + j)];
                int best_rot = results[2*(sz*i + j) + 1];
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
            cnt++;
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
    pb_finalize.finish();
    
    delete slices_info;
    delete results;
    return data;
}