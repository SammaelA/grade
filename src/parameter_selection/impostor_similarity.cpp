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
    glGenBuffers(1, &dbg_buf);
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
    //logerr("imp size %f %f", info.BCyl_sizes.x, info.BCyl_sizes.y);
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

    info.branches_density = 100*res/(V*MAX(b_vox_cnt,1));
    //logerr("branch dens %f %d %f", (float)res, b_vox_cnt, info.branches_density);

    info.branches_curvature = 1 - dots/MAX(seg_cnt,1);
    //logerr("branch curvative %f %d %f", (float)dots, seg_cnt, info.branches_curvature);
    //logerr("imp bcyl size %f %f", info.BCyl_sizes.x, info.BCyl_sizes.y);
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

void ImpostorSimilarityCalc::calc_similarity(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results,
                                             Tree *original_trees, bool debug_print, bool image_debug)
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
    
    if (image_debug)
    { 
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, dbg_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 4*slice_sizes.x*slice_sizes.y*sizeof(glm::uvec4), nullptr, GL_STREAM_READ);
    }

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
    if (image_debug)
       similarity_shader.uniform("image_debug", (int)image_debug); 
    int relative_scale = (reference.width_status == TCIFeatureStatus::FROM_IMAGE && reference.height_status == TCIFeatureStatus::FROM_IMAGE);
    similarity_shader.uniform("relative_scale", relative_scale); 

    glDispatchCompute(cnt, 1, 1);
    //SDL_GL_SwapWindow(Tiny::view.gWindow);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(results_data,ptr,sizeof(float)*cnt);
    
    if (image_debug)
    { 
        glm::uvec4 *data = new glm::uvec4 [4*slice_sizes.x*slice_sizes.y];
        unsigned char *data_ch = new unsigned char[4*4*slice_sizes.x*slice_sizes.y];
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, dbg_buf);
        ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
        memcpy(data,ptr,4*slice_sizes.x*slice_sizes.y*sizeof(glm::uvec4));
        for (int i=0;i<4*slice_sizes.x*slice_sizes.y;i++)
        {
            data_ch[4*i] = data[i].x;
            data_ch[4*i+1] = data[i].y;
            data_ch[4*i+2] = data[i].z;
            data_ch[4*i+3] = data[i].w;
            //debug("%d %d %d %d  ",data[i].x, data[i].y, data[i].z, data[i].w);
        }
        textureManager.save_png_raw(data_ch, slice_sizes.x, 4*slice_sizes.y, 4, "atlas_comp");
        delete[] data;
        delete[] data_ch;
    }

    sim_results = {};
    for (int i=0;i<impostors_cnt;i++)
    {   
        float dist = 0;
        for (int j =0;j<slices_per_impostor;j++)
            dist += results_data[i*slices_per_impostor + j];
        
        dist /= slices_per_impostor;
        dist = smoothstep5(dist);
        float d_ld = abs(impostors_info_data[i+1].leaves_density - impostors_info_data[0].leaves_density); 
        d_ld /= MAX(1e-4, (impostors_info_data[i+1].leaves_density + impostors_info_data[0].leaves_density));
        float d_bd = abs(impostors_info_data[i+1].branches_density - impostors_info_data[0].branches_density); 
        d_bd /= MAX(1e-4, (impostors_info_data[i+1].branches_density + impostors_info_data[0].branches_density));
        float d_bc = abs(impostors_info_data[i+1].branches_curvature - impostors_info_data[0].branches_curvature);
        d_bc /= MAX(1e-4, (impostors_info_data[i+1].branches_curvature + impostors_info_data[0].branches_curvature)); 
        float d_jcnt = abs(impostors_info_data[i+1].joints_cnt - impostors_info_data[0].joints_cnt);
        d_jcnt /= MAX(1, (impostors_info_data[i+1].joints_cnt + impostors_info_data[0].joints_cnt)); 
        
        float d_th;
        if (reference.width_status == TCIFeatureStatus::FROM_IMAGE && 
            reference.height_status == TCIFeatureStatus::FROM_IMAGE &&
            reference.trunk_thickness_status == TCIFeatureStatus::FROM_IMAGE)
        {
            d_th = (impostors_info_data[i+1].trunk_thickness/MAX(impostors_info_data[i+1].BCyl_sizes.y, 1e-6))/
                   (impostors_info_data[0].trunk_thickness/impostors_info_data[0].BCyl_sizes.y);
            if (d_th > 1)
                d_th = 1/d_th;
            d_th = 1 - d_th;
        }
        else
        {
            d_th = abs(impostors_info_data[i+1].trunk_thickness - impostors_info_data[0].trunk_thickness); 
            d_th /= MAX(1e-4, (impostors_info_data[i+1].trunk_thickness + impostors_info_data[0].trunk_thickness)); 
        }

        glm::vec2 scale_fine;
        if (reference.width_status == TCIFeatureStatus::FROM_IMAGE && reference.height_status == TCIFeatureStatus::FROM_IMAGE)
        {
            scale_fine.x = (impostors_info_data[i+1].BCyl_sizes.x/MAX(impostors_info_data[i+1].BCyl_sizes.y, 1e-6))/
                           (impostors_info_data[0].BCyl_sizes.x/impostors_info_data[0].BCyl_sizes.y);
            if (scale_fine.x > 1)
                scale_fine.x = 1/scale_fine.x;
            scale_fine.y = 1;
        }
        else
        {
            scale_fine = impostors_info_data[i+1].BCyl_sizes/impostors_info_data[0].BCyl_sizes; 
            //scale_fine.x *= 2;
            //logerr("bcyl %f %f %f %f", impostors_info_data[i+1].BCyl_sizes.x, impostors_info_data[i+1].BCyl_sizes.y, 
            //       impostors_info_data[0].BCyl_sizes.x, impostors_info_data[0].BCyl_sizes.y);
            if (scale_fine.x > 1)
                scale_fine.x = 1/scale_fine.x;
            if (scale_fine.y > 1)
                scale_fine.y = 1/scale_fine.y;

                    
            if (reference.width_status == TCIFeatureStatus::DONT_CARE)
                scale_fine.x = 1;
            if (reference.height_status == TCIFeatureStatus::DONT_CARE)
                scale_fine.y = 1;
            scale_fine = sqrt(scale_fine);
        }
        if (reference.branches_density_status == TCIFeatureStatus::DONT_CARE)
            d_bd = -1e-5;
        if (reference.leaves_density_status == TCIFeatureStatus::DONT_CARE)
            d_ld = -1e-5;
        if (reference.branches_curvature_status == TCIFeatureStatus::DONT_CARE)
            d_bc = -1e-5;
        if (reference.trunk_thickness_status == TCIFeatureStatus::DONT_CARE)
            d_th = -1e-5;
        if (reference.joints_cnt_status == TCIFeatureStatus::DONT_CARE)
            d_jcnt = -1e-5;
        if (reference.reference_image_status == TCIFeatureStatus::DONT_CARE)
            dist = -1e-5;    
        float d_sd = 1 - (scale_fine.x)*(scale_fine.y);

        if (debug_print && i == 0)
            logerr("dist %f %f %f %f %f %f %f", d_sd, d_ld, d_bd, d_bc, d_jcnt, d_th, dist);
        //dist = CLAMP((1 - d_sd)*(1 - dist)*(1 - d_ld)*(1 - d_bd)*(1 - d_bc)*(1-d_jcnt)*(1-d_th), 0,1);
        dist = CLAMP(((1 - d_sd) + (1 - dist) + (1 - d_ld) + (1 - d_bd) + (1 - d_bc) + (1-d_jcnt) + (1-d_th))/7, 0,1);
        dist = dist*dist*dist;
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
    glDeleteBuffers(1, &dbg_buf);
}