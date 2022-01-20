#include "impostor_similarity.h"
#include "tinyEngine/TinyEngine.h"
#include "graphics_utils/volumetric_occlusion.h"

ImpostorSimilarityCalc::ImpostorSimilarityCalc(int _max_impostors, int _slices_per_impostor, bool _use_top_slice):
similarity_shader({"impostor_atlas_dist.comp"},{})
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

void ImpostorSimilarityCalc::get_tree_compare_info(Impostor &imp, Tree &t, TreeCompareInfo &info)
{
    if (!t.root || !t.valid)
    {
        info.BCyl_sizes = glm::vec2(0,0);
        info.branches_curvature = 0;
        info.branches_density = 0;
        info.joints_cnt = 0;
        info.leaves_density = 0;
        info.trunk_thickness = 0;
        return;
    }
    info.BCyl_sizes = glm::vec2(imp.bcyl.r, 2*imp.bcyl.h_2);
    info.trunk_thickness = t.root->segments.front().rel_r_begin;
    
    glm::vec3 center = t.pos + imp.bcyl.h_2;
    glm::vec3 sz = 1.25f*glm::vec3(imp.bcyl.r, imp.bcyl.h_2, imp.bcyl.r);
    float v_max = MAX(sz.x, MAX(sz.y, sz.z))/8;
    float v_min = MIN(sz.x, MIN(sz.y, sz.z))/4;
    float v_sz = MIN(v_max, v_min);
    //logerr("creating voxels cube size %f %f %f sz %f", sz.x, sz.y, sz.z, v_sz);
    LightVoxelsCube branches_dens = LightVoxelsCube(center, sz, v_sz, 1.0f);
    LightVoxelsCube leaves_dens = LightVoxelsCube(center, sz, v_sz, 1.0f);
    double dots = 0;
    int seg_cnt = 1;

    for (auto &bh : t.branchHeaps)
    {
        for (auto &b : bh->branches)
        {
            if (b.dead)
                continue;
            for (auto &s : b.segments)
            {
                if (s.end.x == s.end.x)
                {
                    branches_dens.set_occluder_simple(s.end, glm::length(s.end - s.begin)*
                    (SQR(s.rel_r_begin) + s.rel_r_begin*s.rel_r_end + SQR(s.rel_r_end)));
                }
                //else 
                //    logerr("NAN detected");
            }
            if (b.segments.size() > 1)
            {
                auto sit = b.segments.begin();
                auto sit2 = b.segments.begin();
                sit2++;
                while (sit2 != b.segments.end())
                {
                    if (sit->end.x == sit->end.x && sit2->end.x == sit2->end.x)
                        dots += glm::dot(normalize(sit->end - sit->begin), normalize(sit2->end - sit2->begin));
                    sit++;
                    sit2++;
                }
                seg_cnt += b.segments.size()-1;
            }
        }
    }

    float V = v_sz*v_sz*v_sz;
    int b_vox_cnt = 1;
    double res = 0;
    std::function<void (float)> f = [&](float val)
    {
        res += val;
        b_vox_cnt += (val >= 1);
    };
    if (t.leaves)
    {
        for (auto &l : t.leaves->leaves)
        {
            float sz = MAX(glm::length(l.edges[0] - l.edges[1]), glm::length(l.edges[0] - l.edges[2]));
            if (l.pos.x == l.pos.x)
                leaves_dens.set_occluder_simple(l.pos, sz*sz);
        }
    }
    leaves_dens.read_func_simple(f);
    info.leaves_density = res/(V*MAX(b_vox_cnt, 1));
    //logerr("leaves dens %f %d %f", (float)res, b_vox_cnt, info.leaves_density);

    b_vox_cnt = 1;
    res = 0;
    branches_dens.read_func_simple(f);

    info.branches_density = res/(V*MAX(b_vox_cnt,1));
    //logerr("branch dens %f %d %f", (float)res, b_vox_cnt, info.branches_density);

    info.branches_curvature = dots/MAX(seg_cnt,1);
    //logerr("branch curvative %f %d %f", (float)dots, seg_cnt, info.branches_curvature);
    //logerr("imp bcyl size %f %f", imp.bcyl.r, 2*imp.bcyl.h_2);
    info.joints_cnt = 0;
    for (auto &bh : t.branchHeaps)
    {
        for (auto &b : bh->branches)
        {
            info.joints_cnt += b.joints.size();
        }
    }
    //logerr("%d joints", info.joints_cnt);
}

void ImpostorSimilarityCalc::calc_similarity(GrovePacked &grove, ReferenceTree &reference, 
                                             std::vector<float> &sim_results, Tree *original_trees)
{
    //grove.impostors[1].atlas.gen_mipmaps("mipmap_render_average.fs");
    
    impostors_info_data[0] = reference.info;
    //logerr("reference info %f %f", reference.info.BCyl_sizes.x, reference.info.BCyl_sizes.y);
    int impostors_cnt = grove.impostors[1].impostors.size();
    int cnt = 0;
    int imp_n = 0;

    for (auto &imp : grove.impostors[1].impostors)
    {
        //logerr("tree %d/%d", imp_n, impostors_cnt);
        get_tree_compare_info(imp, original_trees[imp_n], impostors_info_data[imp_n+1]);

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
    glm::ivec4 atlas_sizes = reference.atlas.get_sizes();
    similarity_shader.use();
    similarity_shader.texture("atlas", grove.impostors[1].atlas.tex(0));
    similarity_shader.texture("reference_image", reference.atlas.tex(0));

    similarity_shader.uniform("reference_images_cnt", atlas_sizes.z);
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
        float d_ld = abs(impostors_info_data[i+1].leaves_density - impostors_info_data[0].leaves_density); 
        d_ld /= MAX(1e-4, (impostors_info_data[i+1].leaves_density + impostors_info_data[0].leaves_density));
        float d_bd = abs(impostors_info_data[i+1].branches_density - impostors_info_data[0].branches_density); 
        d_bd /= MAX(1e-4, (impostors_info_data[i+1].branches_density + impostors_info_data[0].branches_density));
        float d_bc = abs(impostors_info_data[i+1].branches_curvature - impostors_info_data[0].branches_curvature);
        d_bc /= MAX(1e-4, (impostors_info_data[i+1].branches_curvature + impostors_info_data[0].branches_curvature)); 
        float d_jcnt = abs(impostors_info_data[i+1].joints_cnt - impostors_info_data[0].joints_cnt);
        d_jcnt /= MAX(1, (impostors_info_data[i+1].joints_cnt + impostors_info_data[0].joints_cnt)); 
        if (scale_fine.x > 1)
            scale_fine.x = 1/scale_fine.x;
        if (scale_fine.y > 1)
            scale_fine.y = 1/scale_fine.y;
        float sf = sqrt((scale_fine.x)*(scale_fine.y));
        //dist = dist;
        //logerr("dist %f %f %f %f %f %f", 1- sf, d_ld, d_bd, d_bc, d_jcnt, dist);
        dist = CLAMP(sf*(1 - dist)*(1 - d_ld)*(1 - d_bd)*(1 - d_bc)*(1-d_jcnt), 0,1);
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