#include "impostor_similarity.h"
#include "tinyEngine/TinyEngine.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "tinyEngine/postfx.h"
#include "generation/grove_packer.h"

glm::vec4 tc_tr_mult2(glm::vec4 tc_tr_1, glm::vec4 tc_tr_2)
{
    return glm::vec4(tc_tr_2.x + tc_tr_1.x*tc_tr_2.z,
                    tc_tr_2.y + tc_tr_1.y*tc_tr_2.w,
                    tc_tr_1.z*tc_tr_2.z,
                    tc_tr_1.w*tc_tr_2.w);
}

std::vector<std::vector<StripeInfo>> 
ImpostorSimilarityCalc::get_alternative_tree_image_info(TextureAtlas &images_atl, const std::vector<std::pair<int, int>> &slice_id_impostor_n,
                                                        TreeCompareInfo *impostors_info, ReferenceTree *reference, glm::vec4 *tc_transform)
{
    //set slices data
    int slices_cnt = slice_id_impostor_n.size();
    if (slices_cnt > (slices_per_impostor + 1)*max_impostors)
    {
        logerr("atlas for tree image info has too many slices");
        slices_cnt = (slices_per_impostor + 1)*max_impostors;
    }
{
    int cnt = 0;
    for (auto &id : slice_id_impostor_n)
    {
        images_atl.pixel_offsets(id.first, slices_info_data[cnt]);
        cnt++;
        if (cnt > slices_cnt)
            break;
    }
}
    //first step - calculate borders for every slice

    glm::ivec2 slice_size = images_atl.get_slice_size();
    int stripe_height = 16;
    int stripes_cnt = slice_size.y > stripe_height ? slice_size.y/stripe_height : 1;
    if (stripes_cnt*stripe_height != slice_size.y)
    {
        logerr("wrong impostors size %d, it should divide to %d", slice_size.y, stripe_height);
    }
    //bind buffers
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, slices_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::uvec4)*slices_cnt, slices_info_data, GL_STATIC_DRAW);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, stripes_results_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec4)*(2*slices_cnt + stripes_cnt*slices_cnt), nullptr, GL_STREAM_READ);

    //shader to get borders
    constexpr int THREADS = 32;
    Shader get_tree_borders({"get_tree_info_2.comp"},{});
    get_tree_borders.use();
    get_tree_borders.texture("atlas", images_atl.tex(0));

    get_tree_borders.uniform("impostor_x", slice_size.x);
    get_tree_borders.uniform("impostor_y", slice_size.y);
    get_tree_borders.uniform("slices_count", slices_cnt);
    get_tree_borders.uniform("start_id", 0);
    get_tree_borders.uniform("image_debug", 0); 
    get_tree_borders.uniform("stripes_count", stripes_cnt); 
    
    get_tree_borders.uniform("pass", 0); 
    glDispatchCompute(ceil((float)slices_cnt*stripes_cnt/THREADS), 1, 1);
    //SDL_GL_SwapWindow(Tiny::viezw.gWindow);
    glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

    get_tree_borders.uniform("pass", 1); 
    glDispatchCompute(ceil((float)slices_cnt/THREADS), 1, 1);  

    glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

    glm::vec4 *borders_data = new glm::vec4[2*slices_cnt + stripes_cnt*slices_cnt];
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, stripes_results_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(borders_data,ptr,sizeof(glm::vec4)*(2*slices_cnt + stripes_cnt*slices_cnt));

/*
    for (int i=0;i<slices_cnt;i++)
    {
        logerr("borders [%f %f %f %f] center %f",borders_data[2*i].x, borders_data[2*i].y,borders_data[2*i].z,borders_data[2*i].w,borders_data[2*i+1].x);
    }
*/

    if (tc_transform)
    {
        *tc_transform = glm::vec4(borders_data[0].x, borders_data[0].z, borders_data[0].y - borders_data[0].x, borders_data[0].w - borders_data[0].z);
    }
    //transform slices in atlas
    glm::ivec4 sizes = images_atl.get_sizes();
    slice_size = images_atl.get_slice_size();
    TextureAtlas atl_tmp = TextureAtlas(slice_size.x*slice_id_impostor_n.size(), slice_size.y, images_atl.layers_count(), 1);
    atl_tmp.set_grid(slice_size.x, slice_size.y, false);
{
    int slice_n = 0;
    for (auto &id : slice_id_impostor_n)
    {
        glm::vec3 tc = glm::vec3(0,0,0);
        images_atl.process_tc(id.first, tc);
        glm::vec4 border = borders_data[2*slice_n];
        glm::vec4 border_in_atl = glm::vec4(tc.x + border.x/sizes.z, tc.x + border.y/sizes.z, 
                                            tc.y + border.z/sizes.w, tc.y + border.w/sizes.w);//(min_x,max_x,min_y,max_y)
        
        glm::vec2 sz = glm::vec2(0.5,0.5);
        if (reference)
        {
            //scale to the size of reference tree
            float width = impostors_info[id.second].BCyl_sizes.x * (border.y - border.x);
            float height = impostors_info[id.second].BCyl_sizes.y * (border.w - border.z);
            impostors_info[id.second].real_sizes = glm::vec2(width, height);
            sz.x *= width/reference->info.BCyl_sizes.x; 
            sz.y *= height/reference->info.BCyl_sizes.y; 
        }
        else
        {
            impostors_info[id.second].real_sizes = impostors_info[id.second].BCyl_sizes;
        }
        glm::vec2 scale = 1.0f/sz;
        float root_x = borders_data[2*slice_n+1].x;
        float root_x_dst = 0.5;
        float x_shift =  root_x_dst - scale.x*root_x;
        float y_shift = border.z;
        glm::vec4 in_slice_tr = glm::vec4(x_shift,y_shift,scale.x,scale.y);
        glm::vec4 slice_tr = glm::vec4(tc.x, tc.y, 1.0/sizes.z,1.0/sizes.w);
        //logerr("in slice tr %f %f %f %f",in_slice_tr.x, in_slice_tr.y, in_slice_tr.z, in_slice_tr.w);
        //logerr("slice tr %f %f %f %f",slice_tr.x, slice_tr.y, slice_tr.z, slice_tr.w);
        glm::vec4 tex_transform = tc_tr_mult2(in_slice_tr, slice_tr);
        PostFx copy_arr("copy_arr_slice.fs");
        int tex_id = atl_tmp.add_tex();
        atl_tmp.target_slice(tex_id, 0);
        copy_arr.use();
        copy_arr.get_shader().texture("tex", images_atl.tex(0));
        copy_arr.get_shader().uniform("in_slice_tr", in_slice_tr);
        copy_arr.get_shader().uniform("slice_tr", slice_tr);
        copy_arr.get_shader().uniform("layer", tc.z);
        copy_arr.render();


        atl_tmp.pixel_offsets(tex_id, slices_info_data[slice_n]);
        slice_n++; 
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    textureManager.save_png(atl_tmp.tex(0), "alternative_atlas");
}

    //get stripes data
    //bind buffers
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, slices_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::uvec4)*slices_cnt, slices_info_data, GL_STATIC_DRAW);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, stripes_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(StripeInfo)*(stripes_cnt*slices_cnt), nullptr, GL_STREAM_READ);

    //shader to stripes info
    Shader get_tree_stripes({"get_tree_info_stripes.comp"},{});
    get_tree_stripes.use();
    get_tree_stripes.texture("atlas", atl_tmp.tex(0));

    get_tree_stripes.uniform("impostor_x", slice_size.x);
    get_tree_stripes.uniform("impostor_y", slice_size.y);
    get_tree_stripes.uniform("slices_count", slices_cnt);
    get_tree_stripes.uniform("start_id", 0);
    get_tree_stripes.uniform("stripes_count", stripes_cnt); 
    
    get_tree_stripes.uniform("pass", 0); 
    glDispatchCompute(ceil((float)slices_cnt*stripes_cnt/THREADS), 1, 1);
    //SDL_GL_SwapWindow(Tiny::viezw.gWindow);
    glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

    get_tree_stripes.uniform("pass", 1); 
    glDispatchCompute(ceil((float)slices_cnt*stripes_cnt/THREADS), 1, 1);  

    glMemoryBarrier(GL_COMPUTE_SHADER_BIT);

    StripeInfo *stripes_data = new StripeInfo[stripes_cnt*slices_cnt];
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, stripes_info_buf);
    ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(stripes_data,ptr,sizeof(StripeInfo)*(stripes_cnt*slices_cnt));

    /*
    for (int i=0;i<stripes_cnt*slices_cnt;i++)
    {
        logerr("stripe %f %f %f %f -- %f %f -- %f",stripes_data[i].crown_bord.x, stripes_data[i].crown_bord.y, stripes_data[i].crown_bord.z, stripes_data[i].crown_bord.w,
                stripes_data[i].crown_ymin_ymax.x, stripes_data[i].crown_ymin_ymax.y, stripes_data[i].branches_share);
    }
    */
    delete[] borders_data;


    std::vector<std::vector<StripeInfo>> result;
    for (int i=0;i<slices_cnt;i++)
    {
        result.emplace_back();
        for (int j=0;j<stripes_cnt;j++)
        {
            result.back().push_back(stripes_data[i*stripes_cnt + j]);
        }
    }

    delete[] stripes_data;

    return result;
}

void ImpostorSimilarityCalc::get_reference_tree_image_info_alt(ReferenceTree &reference, float clsh_mult )
{
    std::vector<std::pair<int,int>> sl_imp;
    sl_imp.push_back(std::pair<int,int>(0,0));
    auto res = get_alternative_tree_image_info(reference.atlas, sl_imp, &reference.info, nullptr, &reference.image_info.tc_transform);
    reference.alternative_image_info = res[0];
}

void ImpostorSimilarityCalc::calc_similarity_alt(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results,
                                                 Tree *original_trees, int original_trees_cnt, bool debug_print, bool image_debug)
{
    int impostors_cnt = grove.impostors[1].impostors.size();
    int slices_cnt = grove.impostors[1].impostors.size()*(slices_per_impostor);

    int imp_n = 0;
    int t_n = 0;
    impostors_info_data[0] = reference.info;
    tree_image_info_data[0] = reference.image_info;


    std::vector<std::pair<int,int>> sl_imp;

    for (auto &imp : grove.impostors[1].impostors)
    {
        while (t_n < original_trees_cnt && !GrovePacker::is_valid_tree(original_trees[t_n]))
        {
            t_n++;
        }
        if (t_n < original_trees_cnt)
        {
            get_tree_compare_info(imp, original_trees[t_n], impostors_info_data[imp_n+1]);
        }
        else
        {
            auto &info = impostors_info_data[imp_n+1];
            info.BCyl_sizes = glm::vec2(0,0);
            info.branches_curvature = 0;
            info.branches_density = 0;
            info.joints_cnt = 0;
            info.leaves_density = 0;
            info.trunk_thickness = 0; 
        }
        int slices_cnt = 0;
        for (auto &slice : imp.slices)
        {
            sl_imp.push_back(std::pair<int,int>(slice.id,imp_n));
            slices_cnt++;
        }
        imp_n++;
        t_n++;
    }

    auto alt_tree_infos = get_alternative_tree_image_info(grove.impostors[1].atlas, sl_imp, impostors_info_data+1, &reference);


    sim_results = {};
    for (int i=0;i<impostors_cnt;i++)
    {           
        float d_ld = abs(impostors_info_data[i+1].leaves_density - impostors_info_data[0].leaves_density); 
        d_ld /= MAX(1e-4, (impostors_info_data[i+1].leaves_density + impostors_info_data[0].leaves_density));
        float d_bd = abs(impostors_info_data[i+1].branches_density - impostors_info_data[0].branches_density); 
        d_bd /= MAX(1e-4, (impostors_info_data[i+1].branches_density + impostors_info_data[0].branches_density));
        float d_bc = CLAMP(10*abs(impostors_info_data[i+1].branches_curvature - impostors_info_data[0].branches_curvature),0,1);
        float d_jcnt = abs(impostors_info_data[i+1].joints_cnt - impostors_info_data[0].joints_cnt);
        d_jcnt /= MAX(1, (impostors_info_data[i+1].joints_cnt + impostors_info_data[0].joints_cnt)); 
        float d_trop = CLAMP(abs(impostors_info_data[i+1].tropism - impostors_info_data[0].tropism)/2,-1,1);
        
        glm::vec2 scale_fine;
        glm::vec2 imp_real_size = max(glm::vec2(1e-6,1e-6),impostors_info_data[i+1].real_sizes);
        glm::vec2 ref_real_size = max(glm::vec2(1e-6,1e-6),impostors_info_data[0].BCyl_sizes);

        if (reference.width_status == TCIFeatureStatus::FROM_IMAGE && reference.height_status == TCIFeatureStatus::FROM_IMAGE)
        {
            scale_fine.x = (imp_real_size.x/MAX(imp_real_size.y, 1e-6))/
                           (ref_real_size.x/ref_real_size.y);
            if (scale_fine.x > 1)
                scale_fine.x = 1/scale_fine.x;
            scale_fine.y = 1;
        }
        else
        {
            scale_fine = imp_real_size/ref_real_size; 
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
        if (reference.tropism_status == TCIFeatureStatus::DONT_CARE)
            d_trop = -1e-5;
        if (reference.branches_curvature_status == TCIFeatureStatus::DONT_CARE)
            d_bc = -1e-5;
        if (reference.joints_cnt_status == TCIFeatureStatus::DONT_CARE)
            d_jcnt = -1e-5;   
        float d_sd = 1 - (scale_fine.x)*(scale_fine.y);
        float res_dist = 0;
        float alt_imp_dist = 0;

        for (int sl = 0; sl < slices_per_impostor; sl++)
        {
            auto &stripes = alt_tree_infos[i*slices_per_impostor + sl];//TODO: fix for multisliced impostors
            float trunk_sum = 0;
            float trunk_val = 0;
            if (stripes.size() == reference.alternative_image_info.size())
            {
                float max_points = 0;
                float points = 0;
                for (int s = 0; s< stripes.size(); s++)
                {
                    glm::vec4 &b1 = reference.alternative_image_info[s].crown_bord;
                    glm::vec4 &b2 = stripes[s].crown_bord;

                    float max_p_add = MAX((b1.y - b1.x) + 3*(b1.z - b1.y) + (b1.w - b1.z), (b2.y - b2.x) + 3*(b2.z - b2.y) + (b2.w - b2.z));
                    max_points += max_p_add;
                    //logerr("stripe %d max points %f", s, (b1.y - b1.x) + 3*(b1.z - b1.y) + (b1.w - b1.z));
                    float p = 0;

                    p += MAX(0, MIN(b1.y, b2.y) - MAX(b1.x, b2.x)) + MAX(0, MIN(b1.y, b2.z) - MAX(b1.x, b2.y)) + MAX(0, MIN(b1.y, b2.w) - MAX(b1.x, b2.z));
                    p += MAX(0, MIN(b1.z, b2.y) - MAX(b1.y, b2.x)) + 3*MAX(0, MIN(b1.z, b2.z) - MAX(b1.y, b2.y)) + MAX(0, MIN(b1.z, b2.w) - MAX(b1.x, b2.z));
                    p += MAX(0, MIN(b1.w, b2.y) - MAX(b1.z, b2.x)) + MAX(0, MIN(b1.w, b2.z) - MAX(b1.z, b2.y)) + MAX(0, MIN(b1.w, b2.w) - MAX(b1.x, b2.z));

                    glm::vec2 y1 = reference.alternative_image_info[s].crown_ymin_ymax;
                    glm::vec2 y2 = stripes[s].crown_ymin_ymax;
                    float b_share_fine = 3*abs(reference.alternative_image_info[s].branches_share - stripes[s].branches_share);
                    float p_add = p*MAX(0, (MIN(y1.y, y2.y) - MAX(y1.x, y2.x)))/MAX(1, y1.y - y1.x)*MAX(0.01,1 - b_share_fine);
                    points += p_add;
                    //logerr("reference stripe %f %f %f %f -- %f %f -- %f",b1.x,b1.y,b1.z,b1.w,y1.x,y1.y,reference.alternative_image_info[s].branches_share);
                    //logerr("impostor  stripe %f %f %f %f -- %f %f -- %f",b2.x,b2.y,b2.z,b2.w,y2.x,y2.y,stripes[s].branches_share);
                    //logerr("stripe %d max points %f", s, max_p_add);
                    //logerr("stripe %d points %f", s, p_add);

                    if (b1.z - b1.y + 2 > 0.67*(b1.w - b1.x) && reference.alternative_image_info[s].branches_share > 0.9)
                    {
                        //it is probably trunk stripe
                        trunk_sum += max_p_add;
                        trunk_val += p_add;
                    }
                }
                float d = 1 - (0.1 + 0.9*(points/max_points))*(0.1 + 0.9*sqrt(trunk_val/(trunk_sum+0.1)));
                //logerr("altermative image dist %f", d);
                alt_imp_dist += d;
            }
            else
            {
                logerr("reference and impostor have different number of stripes %d %d",stripes.size(), reference.alternative_image_info.size());
                alt_imp_dist += 1;
            }
        }
        alt_imp_dist /= slices_per_impostor;
        alt_imp_dist = 1 - (0.5 - sinf(asinf(1-2*sqrtf(1 - alt_imp_dist))/3));
        res_dist = (1 - d_sd)*(1 - alt_imp_dist)*(1 - d_ld)*(1 - d_bd)*(1 - d_bc)*(1-d_jcnt)*(1-d_trop);
        //logerr("%f %f %f %f %f %f %f", (1 - d_sd),(1 - alt_imp_dist),(1 - d_ld),(1 - d_bd),(1 - d_bc),(1-d_jcnt),(1-d_trop));
        //res_dist = 0.5 - sinf(asinf(1-2*sqrtf(res_dist))/3);
        sim_results.push_back(res_dist);
    }
}