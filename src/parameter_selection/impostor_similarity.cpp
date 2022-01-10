#include "impostor_similarity.h"
#include "tinyEngine/TinyEngine.h"

ImpostorSimilarityCalc::ImpostorSimilarityCalc(int _max_impostors, int _slices_per_impostor, bool _use_top_slice):
similarity_shader({"impostor_image_dist.comp"},{})
{
    use_top_slice = _use_top_slice;
    max_impostors = _max_impostors;
    slices_per_impostor = _slices_per_impostor;
    slices_stride = slices_per_impostor;
    results_data = new float[slices_per_impostor*max_impostors];
    slices_info_data = new glm::uvec4[slices_per_impostor*max_impostors];
    impostors_info_data = new TreeCompareInfo[max_impostors + 1];

    glGenBuffers(1, &results_buf);
    glGenBuffers(1, &slices_info_buf);
    glGenBuffers(1, &impostors_info_buf);
}

void ImpostorSimilarityCalc::get_tree_compare_info(Impostor &imp, TreeCompareInfo &info)
{
    info.BCyl_sizes = glm::vec2(imp.bcyl.r, 2*imp.bcyl.h_2);
    //logerr("imp bcyl size %f %f", imp.bcyl.r, 2*imp.bcyl.h_2);
}

void ImpostorSimilarityCalc::calc_similarity(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results)
{
    impostors_info_data[0] = reference.info;
    logerr("reference info %f %f", reference.info.BCyl_sizes.x, reference.info.BCyl_sizes.y);
    int impostors_cnt = grove.impostors[1].impostors.size();
    int cnt = 0;
    int imp_n = 1;
    for (auto &imp : grove.impostors[1].impostors)
    {
        get_tree_compare_info(imp, impostors_info_data[imp_n]);

        for (auto &sl : imp.slices)
        {
            if (cnt >= slices_per_impostor*max_impostors)
            {
                logerr("more impostor slices than expected at creation of ImpostorSimilarityCalc");
                cnt = slices_per_impostor*max_impostors-1;
            }
            grove.impostors[1].atlas.pixel_offsets(sl.id, slices_info_data[cnt]);
            cnt++;
        }
        if (use_top_slice)
        {
            if (cnt >= slices_per_impostor*max_impostors)
            {
                logerr("more impostor slices than expected at creation of ImpostorSimilarityCalc");
                cnt = slices_per_impostor*max_impostors-1;
            }
            grove.impostors[1].atlas.pixel_offsets(imp.top_slice.id, slices_info_data[cnt]);
            cnt++;
        }
        imp_n++;
    }
    if (cnt != impostors_cnt*slices_per_impostor)
    {
        logerr("wrong impostors count in grove %d, %d = %dx%d expected", cnt, impostors_cnt*slices_per_impostor,
               impostors_cnt, slices_per_impostor);
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, results_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*cnt, nullptr, GL_STREAM_READ);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slices_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::uvec4)*cnt, slices_info_data, GL_STATIC_DRAW);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, impostors_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(TreeCompareInfo)*(impostors_cnt+1), impostors_info_data, GL_STATIC_DRAW);

    glm::ivec2 slice_sizes = grove.impostors[1].atlas.get_slice_size();
    similarity_shader.use();
    similarity_shader.texture("atlas", grove.impostors[1].atlas.tex(0));
    similarity_shader.texture("reference_image", reference.tex);
    similarity_shader.uniform("impostor_x", slice_sizes.x);
    similarity_shader.uniform("impostor_y", slice_sizes.y);
    similarity_shader.uniform("impostors_count", impostors_cnt);
    similarity_shader.uniform("impostor_slice_count", slices_per_impostor);
    similarity_shader.uniform("slice_stride", slices_stride);
    similarity_shader.uniform("start_id", 0);
    
    glDispatchCompute(cnt, 1, 1);
    //SDL_GL_SwapWindow(Tiny::view.gWindow);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(results_data,ptr,sizeof(float)*cnt);
    //logerr("cnt = %d", cnt);
    sim_results = {};
    for (int i=0;i<impostors_cnt;i++)
    {   
        float dist = 0;
        for (int j =0;j<slices_per_impostor;j++)
        {
            dist += results_data[i*slices_per_impostor + j];
            //logerr("dist %d %d %f", i, j, results_data[i*slices_per_impostor + j]);
        }
        dist /= slices_per_impostor;
        glm::vec2 scale_fine = impostors_info_data[i+1].BCyl_sizes/impostors_info_data[0].BCyl_sizes; 
        if (scale_fine.x > 1)
            scale_fine.x = 1/scale_fine.x;
        if (scale_fine.y > 1)
            scale_fine.y = 1/scale_fine.y;
        float sf = (scale_fine.x)*(scale_fine.y);
        //logerr("scale fine %f %f %f", scale_fine.x, scale_fine.y, sf);
        dist = CLAMP(sf*(1 - dist), 0,1);
        sim_results.push_back(dist);
        //logerr("similarity data %f", sim_results.back());
    }
}

ImpostorSimilarityCalc::~ImpostorSimilarityCalc()
{
    if (results_data)
        delete[] results_data;
    if (slices_info_data)
        delete[] slices_info_data;
    if (impostors_info_data)
        delete[] impostors_info_data;
    glDeleteFramebuffers(1, &fbo);
    glDeleteBuffers(1, &results_buf);
    glDeleteBuffers(1, &slices_info_buf);
    glDeleteBuffers(1, &impostors_info_buf);
}