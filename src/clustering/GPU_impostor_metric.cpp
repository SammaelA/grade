#include "GPU_impostor_metric.h"
#include "dist_data_table.h"
#include "../billboard_cloud.h"
#include "graphics_utils/texture_manager.h"
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
    glm::vec4 *branches_sizes = new glm::vec4[sz];
    
    int gpu_impostor_start_mip = get_default_block().get_int("gpu_impostor_start_mip", 4);
    gpu_impostor_start_mip = settings.get_int("gpu_impostor_start_mip", gpu_impostor_start_mip);
    
    bool gpu_impostor_use_mips = get_default_block().get_bool("gpu_impostor_use_mips", true);
    gpu_impostor_use_mips = settings.get_bool("gpu_impostor_use_mips", gpu_impostor_use_mips);
    logerr("usin mip %d", gpu_impostor_use_mips);

    auto &atlas = ictx->self_impostors_data->atlas;
    if (gpu_impostor_use_mips)
    {
        atlas.gen_mipmaps("mipmap_render_average.fs");
        for (int i = 0; i < sz; i++)
        {
            branches_sizes[i] = glm::vec4(real_branches[i]->sizes,sqrt(SQR(real_branches[i]->sizes.y) + SQR(real_branches[i]->sizes.z)));
        }
    }
    
    int p = 0;
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
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*2*SQR(sz), nullptr, GL_STREAM_READ);

        glGenBuffers(1, &slices_info_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slices_info_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::uvec4)*slices_cnt, slices_info, GL_STATIC_DRAW);

        if (gpu_impostor_use_mips)
        {
            glGenBuffers(1, &branches_sizes_buf);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, branches_sizes_buf);
            glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec4)*sz, branches_sizes, GL_STATIC_DRAW);
        }

        delete slices_info;
        delete branches_sizes;

    pb_buffers.finish();
        
        std::string metric_shader_name = gpu_impostor_use_mips ? "impostor_dist_mips.comp": "impostor_dist.comp";
        Shader impMetric({metric_shader_name},{});
            impMetric.use();
            impMetric.uniform("tex_sz_x", (float)sizes.x);
            impMetric.uniform("tex_sz_y", (float)sizes.y);
            impMetric.uniform("slice_count",tasks_cnt_per_impostor);
            impMetric.uniform("branch_count",sz);
            impMetric.uniform("impostor_x",sizes.x/sizes.z);
            impMetric.uniform("impostor_y",sizes.y/sizes.w);
            impMetric.texture("atlas", atlas.tex(0));
            
            if (gpu_impostor_use_mips)
            {
                impMetric.uniform("start_mip", MIN((int)log2(sizes.x/sizes.z), gpu_impostor_start_mip));
                impMetric.uniform("max_dist", 1.0f);
                impMetric.uniform("size_diff_tolerance",isimParams.size_diff_tolerance);
                impMetric.uniform("size_diff_factor",isimParams.size_diff_factor);
            }

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
            glDeleteBuffers(1, &branches_sizes_buf);

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
                if (!gpu_impostor_use_mips)
                {
                    #define SZ_DIFF(a,b) pow(MAX(1, MAX(a,b)/MIN(a,b) - isimParams.size_diff_tolerance), isimParams.size_diff_factor)
                    glm::vec3 &s1 = real_branches[i]->sizes;
                    glm::vec3 &s2 = real_branches[j]->sizes;
                    float dist_discriminator = SZ_DIFF(s1.x, s2.x) *
                                            SZ_DIFF(sqrt(SQR(s1.y) + SQR(s1.z)), sqrt(SQR(s2.y) + SQR(s2.z)));
                    min_av_dist += dist_discriminator - 1;
                    if (min_av_dist > 1)
                        min_av_dist = 1e9;
                }
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
    
    delete results;
    return data;
}