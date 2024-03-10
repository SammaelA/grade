#include "impostor_similarity.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "tinyEngine/postfx.h"
#include "generation/grove_packer.h"

ImpostorSimilarityCalc::ImpostorSimilarityCalc(int _max_impostors, int _slices_per_impostor, bool _use_top_slice):
similarity_shader({"impostor_atlas_dist.comp"},{}),
tree_info_shader({"get_tree_info.comp"},{})
{
    use_top_slice = _use_top_slice;
    max_impostors = _max_impostors;
    slices_per_impostor = _slices_per_impostor;
    slices_stride = slices_per_impostor;
    results_data = new float[slices_per_impostor*max_impostors];
    slices_info_data = new uint4[(slices_per_impostor  + 1)*max_impostors];
    impostors_info_data = new TreeCompareInfo[max_impostors + 1];
    shader_imp_data = new float4[max_impostors + 1];
    tree_image_info_data = new TreeImageInfo[(slices_per_impostor + 1)*max_impostors];

    results_buf = create_buffer();
    slices_info_buf = create_buffer();
    impostors_info_buf = create_buffer();
    dbg_buf = create_buffer();
    tree_image_info_buf = create_buffer();
    stripes_results_buf = create_buffer();
    stripes_info_buf = create_buffer();
}

void ImpostorSimilarityCalc::get_tree_compare_info(Impostor &imp, Tree &t, TreeCompareInfo &info)
{
    if (!t.root || !t.valid)
    {
        info.BCyl_sizes = float2(0,0);
        info.branches_curvature = 0;
        info.branches_density = 0;
        info.joints_cnt = 0;
        info.leaves_density = 0;
        info.trunk_thickness = 0;
        return;
    }
    info.BCyl_sizes = float2(imp.bcyl.r, 2*imp.bcyl.h_2);
    //logerr("imp size %f %f", info.BCyl_sizes.x, info.BCyl_sizes.y);
    info.trunk_thickness = t.root->segments.front().rel_r_begin;
    
    float3 center = t.pos + imp.bcyl.h_2;
    float3 sz = 1.25f*float3(imp.bcyl.r, imp.bcyl.h_2, imp.bcyl.r);
    float v_sz = MAX(sz.x, MAX(sz.y, sz.z))/32;
    //logerr("creating voxels cube size %f %f %f sz %f", sz.x, sz.y, sz.z, v_sz);
    LightVoxelsCube branches_dens = LightVoxelsCube(center, sz, v_sz);
    LightVoxelsCube leaves_dens = LightVoxelsCube(center, sz, v_sz);
    double dots = 0;
    int seg_cnt = 1;

    double trop_sum = 0;
    int seg_for_trop_cnt = 1;
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
                    branches_dens.set_occluder_simple(s.end, length(s.end - s.begin)*
                    (SQR(s.rel_r_begin) + s.rel_r_begin*s.rel_r_end + SQR(s.rel_r_end)));
                    leaves_dens.set_occluder_simple(s.end, 0.001);
                }
                //else 
                //    logerr("NAN detected");
            }
            if (b.segments.size() > 1)
            {
                auto sit = b.segments.begin();
                auto sit2 = b.segments.begin();
                auto jit = b.joints.begin();
                jit++;
                sit2++;
                float b_trop = 0;
                int b_chb = 0;
                while (sit2 != b.segments.end())
                {
                    if (sit->end.x == sit->end.x && sit2->end.x == sit2->end.x)
                    {
                        float l1 = MAX(length(sit->end - sit->begin),1e-6);
                        float l2 = MAX(length(sit2->end - sit2->begin),1e-6);
                        dots += dot((sit->end - sit->begin)/l1, (sit2->end - sit2->begin)/l2);
                        b_trop += dot((sit2->end - sit2->begin)/l2, float3(0,1,0)) - dot((sit->end - sit->begin)/l1, float3(0,1,0));
                    }
                    b_chb += jit->childBranches.size();
                    sit++;
                    sit2++;
                    jit++;
                }
                seg_cnt += b.segments.size()-1;
                if (b_chb > 0)
                {
                    //we do not use small branches(with no child branches) for tropism, as they are usually have random direction
                    trop_sum += b_trop;
                    seg_for_trop_cnt += b.segments.size()-1;
                }
            }
        }
    }

    float V = v_sz*v_sz*v_sz;
    int b_vox_cnt = 1;
    double res = 0;
    std::function<void (float)> f = [&](float val)
    {
        res += val;
        b_vox_cnt += (val >= 0.001);
    };
    if (t.leaves)
    {
        for (auto &l : t.leaves->leaves)
        {
            int cnt = l.edges.size()/4;
            //logerr("l.edges %d",cnt);
            for (int i=0;i<cnt;i++)
            {
                float occ = length(l.edges[4*i] - l.edges[4*i+1])*length(l.edges[4*i] - l.edges[4*i+2]);
                if (l.edges[4*i].x == l.edges[4*i].x)
                {
                    leaves_dens.set_occluder_simple(l.edges[4*i], occ);
                    branches_dens.set_occluder_simple(l.edges[4*i], 0.001);
                }
            }
        }
    }
    info.id = t.id;
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
    info.tropism = 10 * trop_sum/MAX(seg_for_trop_cnt,1);
    //logerr("tree tropism %f %d %f", trop_sum, seg_for_trop_cnt, info.tropism);
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

void print_tree_image_info(TreeImageInfo &info)
{
         debug("tr thick %f %f\n",info.trunk_thickness,info.trunk_info.average_thickness);
        debug("tr end sym_line len %f %f %f\n",info.trunk_info.trunk_split_level,
                                               info.trunk_info.symmetry_x,
                                               info.trunk_info.trunk_len);
        debug("thickness [ ");
        for (int j=0;j<16;j++)
        {
            debug("%f ",info.trunk_info.thickness[j]);
        }
        debug("]\n");   
}

void ImpostorSimilarityCalc::get_tree_image_info(TextureAtlas &images_atl, std::map<int, TreeImageInfo> &results, 
                                                 bool image_debug)
{
    //image_debug = true;
    int4 sizes = images_atl.get_sizes();
    int2 slice_size = images_atl.get_slice_size();
    TextureAtlas atl_tmp = TextureAtlas(sizes.x, sizes.y, images_atl.layers_count(), 1);
    atl_tmp.set_grid(slice_size.x, slice_size.y, false);
    PostFx copy = PostFx("copy_arr2.fs");
    for (int i = 0; i < images_atl.layers_count(); i++)
    {
        int tex_id = atl_tmp.add_tex();
        atl_tmp.target(tex_id, 0);
        copy.use();
        copy.get_shader().texture("tex", images_atl.tex(0));
        copy.get_shader().uniform("tex_transform", float4(0,0,1,1));
        copy.get_shader().uniform("layer", (float)i);
        copy.render();
    }
    ref_atlas_transform(images_atl);

    std::vector<int> ids = images_atl.get_all_valid_slices_ids();
    int cnt = 0;
    for (auto &id : ids)
    {
        images_atl.pixel_offsets(id, slices_info_data[cnt]);
        cnt++;
    }

    int slices_cnt = ids.size();
    if (slices_cnt > (slices_per_impostor + 1)*max_impostors)
    {
        logerr("atlas for tree image info has too many slices");
        slices_cnt = (slices_per_impostor + 1)*max_impostors;
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, tree_image_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(TreeImageInfo)*slices_cnt, nullptr, GL_STREAM_READ);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, slices_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint4)*slices_cnt, slices_info_data, GL_STATIC_DRAW);

    int2 slice_sizes = images_atl.get_slice_size();
    int4 atlas_sizes = images_atl.get_sizes();
    
    if (image_debug)
    { 
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, dbg_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 4*slice_sizes.x*slice_sizes.y*sizeof(uint4), nullptr, GL_STREAM_READ);
    }

    tree_info_shader.use();
    tree_info_shader.texture("initial_atlas", atl_tmp.tex(0));
    tree_info_shader.texture("transformed_atlas", images_atl.tex(0));

    tree_info_shader.uniform("impostor_x", slice_sizes.x);
    tree_info_shader.uniform("impostor_y", slice_sizes.y);
    tree_info_shader.uniform("slices_count", slices_cnt);
    tree_info_shader.uniform("start_id", 0);
    tree_info_shader.uniform("image_debug", (int)image_debug); 

    glDispatchCompute(ceil(slices_cnt), 1, 1);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, tree_image_info_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(tree_image_info_data,ptr,sizeof(TreeImageInfo)*slices_cnt);
        
    if (image_debug)
    { 
        for (int i=0;i<slices_cnt;i++)
        {
            print_tree_image_info(tree_image_info_data[i]);
            logerr("line %f %f %f %f", tree_image_info_data[i].crown_start_level, tree_image_info_data[i].crown_branches_share, 
                                       tree_image_info_data[i].crown_leaves_share, tree_image_info_data[i].trunk_thickness);
            float4 tc = tree_image_info_data[i].tc_transform;
            logerr("tc transform %f %f %f %f", tc.x, tc.y, tc.z, tc.w);
        }
        uint4 *data = new uint4[4*slice_sizes.x*slice_sizes.y];
        unsigned char *data_ch = new unsigned char[4*4*slice_sizes.x*slice_sizes.y];
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, dbg_buf);
        ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
        memcpy(data,ptr,4*slice_sizes.x*slice_sizes.y*sizeof(uint4));
        for (int i=0;i<4*slice_sizes.x*slice_sizes.y;i++)
        {
            data_ch[4*i] = data[i].x;
            data_ch[4*i+1] = data[i].y;
            data_ch[4*i+2] = data[i].z;
            data_ch[4*i+3] = data[i].w;
            //debug("%d %d %d %d  ",data[i].x, data[i].y, data[i].z, data[i].w);
        }
        engine::textureManager->save_png_raw(data_ch, slice_sizes.x, 4*slice_sizes.y, 4, "get_tree_info_debug");
        delete[] data;
        delete[] data_ch;
    }

    for (int i=0;i<slices_cnt;i++)
    {
        results.emplace(ids[i], tree_image_info_data[i]);
    }
}

void ImpostorSimilarityCalc::set_slices_data(GrovePacked &grove)
{
    int impostors_cnt = grove.impostors[1].impostors.size();
    int slices_cnt = grove.impostors[1].impostors.size()*(slices_per_impostor);
    int cnt = 0;
    int imp_n = 0;

    for (auto &imp : grove.impostors[1].impostors)
    {
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
            logerr("use_top_slice is not supported");
            /*
            if (cnt >= slices_per_impostor*max_impostors)
            {
                logerr("more impostor slices than expected at creation of ImpostorSimilarityCalc");
                cnt = slices_per_impostor*max_impostors-1;
            }
            grove.impostors[1].atlas.pixel_offsets(imp.top_slice.id, slices_info_data[cnt]);
            cnt++;
            */
        }
        imp_n++;
    }
    if (cnt != impostors_cnt*slices_per_impostor)
    {
        logerr("wrong impostors count in grove %d, %d = %dx%d expected", cnt, impostors_cnt*slices_per_impostor,
               impostors_cnt, slices_per_impostor);
    }
}

void ImpostorSimilarityCalc::get_reference_tree_image_info(ReferenceTree &reference, float clsh_mult)
{
    get_reference_tree_image_info_alt(reference, clsh_mult);
    return;

    std::map<int, TreeImageInfo> results;
    get_tree_image_info(reference.atlas, results, true);
    if (results.empty())
    {
        logerr("refrence trees atlas is empty!!!");
        return;
    }
    TreeImageInfo av_info;
    av_info.tc_transform = float4(0,0,0,0);
    for (int i=0;i<16;i++)
        av_info.trunk_info.thickness[i] = 0;
    for (auto &p : results)
    {
        av_info.trunk_thickness += p.second.trunk_thickness;
        av_info.crown_start_level += p.second.crown_start_level;
        av_info.crown_leaves_share += p.second.crown_leaves_share;
        av_info.crown_branches_share += p.second.crown_branches_share;
        av_info.tc_transform += p.second.tc_transform;

        av_info.trunk_info.average_thickness += p.second.trunk_info.average_thickness;
        av_info.trunk_info.symmetry_x += p.second.trunk_info.symmetry_x;
        av_info.trunk_info.trunk_len += p.second.trunk_info.trunk_len;
        av_info.trunk_info.trunk_split_level += p.second.trunk_info.trunk_split_level;
        for (int i=0;i<16;i++)
            av_info.trunk_info.thickness[i] += p.second.trunk_info.thickness[i];
    }
    av_info.trunk_thickness /= results.size();
    av_info.crown_start_level /= results.size();
    av_info.crown_leaves_share /= results.size();
    av_info.crown_branches_share /= results.size();
    av_info.tc_transform /= results.size();

    av_info.trunk_info.average_thickness /= results.size();
    av_info.trunk_info.symmetry_x /= results.size();
    av_info.trunk_info.trunk_len /= results.size();
    av_info.trunk_info.trunk_split_level /= results.size();
    for (int i=0;i<16;i++)
        av_info.trunk_info.thickness[i] /= results.size();

    av_info.crown_leaves_share *= clsh_mult;//usually we don't see all the gaps on tree's crone on image, so it has larger leaves share
                                            //it is an euristic to reduce this error
    reference.image_info = av_info;
}

float sigmoid(float x)
{
    #define SIGM(x,k) (1/(1 + exp(-k*(x))));
    constexpr int k = 5;
    float s_1 = SIGM(-1,k);
    return s_1 + SIGM(2*x-1,k);
}

void ImpostorSimilarityCalc::calc_similarity(GrovePacked &grove, ReferenceTree &reference, std::vector<float> &sim_results,
                                             Tree *original_trees, int original_trees_cnt, bool debug_print, bool image_debug)
{
    calc_similarity_alt(grove, reference, sim_results, original_trees, original_trees_cnt, debug_print, image_debug);
    return;
    
    int impostors_cnt = grove.impostors[1].impostors.size();
    int slices_cnt = grove.impostors[1].impostors.size()*(slices_per_impostor);

    std::map<int, TreeImageInfo> images_info;
    //engine::textureManager->save_png(grove.impostors[1].atlas.tex(0),"atlas_original");
    get_tree_image_info(grove.impostors[1].atlas, images_info, false);

    int imp_n = 0;
    int t_n = 0;
    impostors_info_data[0] = reference.info;
    tree_image_info_data[0] = reference.image_info;
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
            info.BCyl_sizes = float2(0,0);
            info.branches_curvature = 0;
            info.branches_density = 0;
            info.joints_cnt = 0;
            info.leaves_density = 0;
            info.trunk_thickness = 0; 
        }
        int slices_cnt = 0;
        TreeImageInfo av_info;
        av_info.tc_transform = float4(0,0,0,0);
        for (int i=0;i<16;i++)
            av_info.trunk_info.thickness[i] = 0;
        for (auto &slice : imp.slices)
        {
            auto it = images_info.find(slice.id);
            if (it == images_info.end())
            {
                logerr("cannot find image info for slice %d", slice.id);
            }
            else
            {
                av_info.trunk_thickness += it->second.trunk_thickness;
                av_info.crown_start_level += it->second.crown_start_level;
                av_info.crown_leaves_share += it->second.crown_leaves_share;
                av_info.crown_branches_share += it->second.crown_branches_share;
                av_info.tc_transform += it->second.tc_transform;

                av_info.trunk_info.average_thickness += it->second.trunk_info.average_thickness;
                av_info.trunk_info.symmetry_x += it->second.trunk_info.symmetry_x;
                av_info.trunk_info.trunk_len += it->second.trunk_info.trunk_len;
                av_info.trunk_info.trunk_split_level += it->second.trunk_info.trunk_split_level;
                for (int i=0;i<16;i++)
                    av_info.trunk_info.thickness[i] += it->second.trunk_info.thickness[i];
                
                slices_cnt++;
            }
        }
        av_info.trunk_thickness /= slices_cnt;
        av_info.crown_start_level /= slices_cnt;
        av_info.crown_leaves_share /= slices_cnt;
        av_info.crown_branches_share /= slices_cnt;
        av_info.tc_transform /= slices_cnt;

        av_info.trunk_info.average_thickness /= slices_cnt;
        av_info.trunk_info.symmetry_x /= slices_cnt;
        av_info.trunk_info.trunk_len /= slices_cnt;
        av_info.trunk_info.trunk_split_level /= slices_cnt;
        for (int i=0;i<16;i++)
            av_info.trunk_info.thickness[i] /= slices_cnt;
        tree_image_info_data[imp_n+1] = av_info;
        imp_n++;
        t_n++;
    }
    for (int i=0;i<imp_n+1;i++)
    {
        if (i!=0)
        {
        impostors_info_data[i].BCyl_sizes = float2(tree_image_info_data[i].tc_transform.z*impostors_info_data[i].BCyl_sizes.x,
                                                      tree_image_info_data[i].tc_transform.w*impostors_info_data[i].BCyl_sizes.y);
        }
        shader_imp_data[i].x = impostors_info_data[i].BCyl_sizes.x;
        shader_imp_data[i].y = impostors_info_data[i].BCyl_sizes.y;
        shader_imp_data[i].z = tree_image_info_data[i].crown_start_level;
        shader_imp_data[i].w = tree_image_info_data[i].trunk_info.symmetry_x;
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    set_slices_data(grove);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, results_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*slices_cnt, nullptr, GL_STREAM_READ);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slices_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint4)*slices_cnt, slices_info_data, GL_STATIC_DRAW);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, impostors_info_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float4)*(impostors_cnt+1), shader_imp_data, GL_STATIC_DRAW);

    int2 slice_sizes = grove.impostors[1].atlas.get_slice_size();
    int4 atlas_sizes = reference.atlas.get_sizes();
    
    if (image_debug)
    { 
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, dbg_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 4*slice_sizes.x*slice_sizes.y*sizeof(uint4), nullptr, GL_STREAM_READ);
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
    similarity_shader.uniform("image_debug", (int)image_debug); 
    int relative_scale = (reference.width_status == TCIFeatureStatus::FROM_IMAGE && reference.height_status == TCIFeatureStatus::FROM_IMAGE);
    similarity_shader.uniform("relative_scale", relative_scale); 

    glDispatchCompute(slices_cnt, 1, 1);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_buf);
    GLvoid* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
    memcpy(results_data,ptr,sizeof(float)*slices_cnt);
    
    if (image_debug)
    { 
        uint4 *data = new uint4 [4*slice_sizes.x*slice_sizes.y];
        unsigned char *data_ch = new unsigned char[4*4*slice_sizes.x*slice_sizes.y];
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, dbg_buf);
        ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
        memcpy(data,ptr,4*slice_sizes.x*slice_sizes.y*sizeof(uint4));
        for (int i=0;i<4*slice_sizes.x*slice_sizes.y;i++)
        {
            data_ch[4*i] = data[i].x;
            data_ch[4*i+1] = data[i].y;
            data_ch[4*i+2] = data[i].z;
            data_ch[4*i+3] = data[i].w;
            //debug("%d %d %d %d  ",data[i].x, data[i].y, data[i].z, data[i].w);
        }
        engine::textureManager->save_png_raw(data_ch, slice_sizes.x, 4*slice_sizes.y, 4, "atlas_comp");
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
        //dist = sigmoid(dist);
        //dist = 1 - SQR(1-dist);
        
        float d_ld = abs(impostors_info_data[i+1].leaves_density - impostors_info_data[0].leaves_density); 
        d_ld /= MAX(1e-4, (impostors_info_data[i+1].leaves_density + impostors_info_data[0].leaves_density));
        float d_bd = abs(impostors_info_data[i+1].branches_density - impostors_info_data[0].branches_density); 
        d_bd /= MAX(1e-4, (impostors_info_data[i+1].branches_density + impostors_info_data[0].branches_density));
        float d_bc = abs(impostors_info_data[i+1].branches_curvature - impostors_info_data[0].branches_curvature);
        d_bc /= MAX(1e-4, (impostors_info_data[i+1].branches_curvature + impostors_info_data[0].branches_curvature)); 
        float d_jcnt = abs(impostors_info_data[i+1].joints_cnt - impostors_info_data[0].joints_cnt);
        d_jcnt /= MAX(1, (impostors_info_data[i+1].joints_cnt + impostors_info_data[0].joints_cnt)); 
        float d_trop = CLAMP(abs(impostors_info_data[i+1].tropism - impostors_info_data[0].tropism)/2,-1,1);
        
        float d_b_crone = abs(tree_image_info_data[i+1].crown_branches_share - tree_image_info_data[0].crown_branches_share);
        float d_b_leaves = abs(tree_image_info_data[i+1].crown_leaves_share - tree_image_info_data[0].crown_leaves_share);
        float d_cs = CLAMP(4*abs(tree_image_info_data[i+1].crown_start_level - tree_image_info_data[0].crown_start_level),0,1);
        float d_tl = CLAMP(abs(tree_image_info_data[i+1].trunk_info.trunk_len - tree_image_info_data[0].trunk_info.trunk_len) /
                           (1e-4 + tree_image_info_data[0].trunk_info.trunk_len + tree_image_info_data[i+1].trunk_info.trunk_len),
                           0,1);
        float d_th = 0;
        float v1=1e-6,v2=1e-6;
        for (int j=0;j<16;j++)
        {
           float r1 = tree_image_info_data[0].trunk_info.thickness[j];
           float r2 = tree_image_info_data[i+1].trunk_info.thickness[j]; 
           v1 += r1*r1;
           v2 += r2*r2;
           d_th += abs(r1*r1 - r2*r2);
        }
        d_th /= MAX(v1,v2);
        //logerr("th %f %f", tree_image_info_data[i+1].trunk_thickness, tree_image_info_data[0].trunk_thickness);
        float2 scale_fine;
        float2 ref_real_size = (1-2*ReferenceTree::border_size)*
                                   impostors_info_data[0].BCyl_sizes;//reference image has empty borders 
        float2 imp_real_size = float2(tree_image_info_data[i+1].tc_transform.z*impostors_info_data[i+1].BCyl_sizes.x,
                                            tree_image_info_data[i+1].tc_transform.w*impostors_info_data[i+1].BCyl_sizes.y);
        imp_real_size = impostors_info_data[i+1].BCyl_sizes;
        ref_real_size = impostors_info_data[0].BCyl_sizes;
        //logerr("%f %f -- %f %f", ref_real_size.x, ref_real_size.y, imp_real_size.x, imp_real_size.y);
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
        if (reference.tropism_status == TCIFeatureStatus::DONT_CARE)
            d_trop = -1e-5;
        if (reference.branches_curvature_status == TCIFeatureStatus::DONT_CARE)
            d_bc = -1e-5;
        //if (reference.trunk_thickness_status == TCIFeatureStatus::DONT_CARE)
        //    d_th = -1e-5;
        if (reference.joints_cnt_status == TCIFeatureStatus::DONT_CARE)
            d_jcnt = -1e-5;
        if (reference.reference_image_status == TCIFeatureStatus::DONT_CARE)
            dist = -1e-5;    
        float d_sd = 1 - (scale_fine.x)*(scale_fine.y);

        if (debug_print && i == 0)
        {
            logerr("size %.3f l_dens %.3f b_dens %.3f b_curv %.3f b_cnt %.3f tr_len %.3f tr_th %.3f trop %.3f crone_st %.3f l_share %.3f b_share %.3f imp_sim %.3f", 
                   d_sd, d_ld, d_bd, d_bc, d_jcnt, d_tl, d_th, d_trop, d_b_crone, 
                   d_b_leaves, d_cs, dist);
        }
        float res_dist = CLAMP(((1 - d_sd) + (1 - d_ld) + (1 - d_bd) + (1 - d_bc) + (1-d_jcnt) + (1-d_tl) + (1-d_th) + 
                                (1-d_trop) + (1 - d_b_crone) + (1 - d_b_leaves) + (1 - d_cs) + 5*(1 - dist))/16, 0,1);
        res_dist = res_dist*res_dist*res_dist;
        
        bool valid = d_sd < 0.8 && d_jcnt < 0.9 && dist < 0.9;
        if (!valid)
            res_dist = 0;

        sim_results.push_back(res_dist);
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
    if (shader_imp_data)
        delete[] shader_imp_data;
    if (tree_image_info_data)
        delete[] tree_image_info_data;
    delete_framebuffer(fbo);
    delete_buffer(results_buf);
    delete_buffer(slices_info_buf);
    delete_buffer(impostors_info_buf);
    delete_buffer(dbg_buf);
    delete_buffer(tree_image_info_buf);
    delete_buffer(stripes_results_buf);
    delete_buffer(stripes_info_buf);
}

void ImpostorSimilarityCalc::ref_atlas_transform(TextureAtlas &atl)
{
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    int4 sizes = atl.get_sizes();
    int2 slice_size = atl.get_slice_size();
    TextureAtlas atl_tmp = TextureAtlas(sizes.x, sizes.y, atl.layers_count(), 1);
    atl_tmp.set_grid(slice_size.x, slice_size.y, false);

    PostFx gauss = PostFx("gaussian_blur_atlas.fs");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl_tmp.target(l, 0);
        gauss.use();
        gauss.get_shader().texture("tex", atl.tex(0));
        gauss.get_shader().uniform("tex_transform", float4(0, 0, 1, 1));
        gauss.get_shader().uniform("layer", (float)l);
        gauss.get_shader().uniform("pass", 0);
        gauss.get_shader().uniform("tex_size_inv", float2(1.0f / sizes.x, 1.0f / sizes.y));
        gauss.get_shader().uniform("slice_size", to_float2(slice_size));
        gauss.get_shader().uniform("blur_step", 0.33f);
        gauss.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //engine::textureManager->save_png(atl_tmp.tex(0), "atlass_gauss_0");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl.target(l, 0);
        gauss.use();
        gauss.get_shader().texture("tex", atl_tmp.tex(0));
        gauss.get_shader().uniform("tex_transform", float4(0, 0, 1, 1));
        gauss.get_shader().uniform("layer", l);
        gauss.get_shader().uniform("pass", 1);
        gauss.get_shader().uniform("tex_size_inv", float2(1.0f / sizes.x, 1.0f / sizes.y));
        gauss.get_shader().uniform("slice_size", to_float2(slice_size));
        gauss.get_shader().uniform("blur_step", 0.33f);
        gauss.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //engine::textureManager->save_png(atl.tex(0), "atlass_gauss_1");

    PostFx sil_fill = PostFx("silhouette_fill.fs");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl_tmp.target(l, 0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl.tex(0));
        sil_fill.get_shader().uniform("tex_transform", float4(0, 0, 1, 1));
        sil_fill.get_shader().uniform("layer", (float)l);
        sil_fill.get_shader().uniform("radius", 4);
        sil_fill.get_shader().uniform("dir_threshold", 4);
        sil_fill.get_shader().uniform("tex_size_inv", float2(1.0f / sizes.x, 1.0f / sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", to_float2(slice_size));
        sil_fill.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl.target(l, 0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl_tmp.tex(0));
        sil_fill.get_shader().uniform("tex_transform", float4(0, 0, 1, 1));
        sil_fill.get_shader().uniform("layer", (float)l);
        sil_fill.get_shader().uniform("radius", 4);
        sil_fill.get_shader().uniform("dir_threshold", 4);
        sil_fill.get_shader().uniform("tex_size_inv", float2(1.0f / sizes.x, 1.0f / sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", to_float2(slice_size));
        sil_fill.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl_tmp.target(l, 0);
        sil_fill.use();
        sil_fill.get_shader().texture("tex", atl.tex(0));
        sil_fill.get_shader().uniform("tex_transform", float4(0, 0, 1, 1));
        sil_fill.get_shader().uniform("layer", (float)l);
        sil_fill.get_shader().uniform("radius", 8);
        sil_fill.get_shader().uniform("dir_threshold", 6);
        sil_fill.get_shader().uniform("tex_size_inv", float2(1.0f / sizes.x, 1.0f / sizes.y));
        sil_fill.get_shader().uniform("threshold", 0.05f);
        sil_fill.get_shader().uniform("slice_size", to_float2(slice_size));
        sil_fill.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //engine::textureManager->save_png(atl_tmp.tex(0), "atlass_gauss_2");

    PostFx sil_sharp = PostFx("silhouette_sharpen.fs");
    for (int l = 0; l < atl.layers_count(); l++)
    {
        atl.target(l, 0);
        sil_sharp.use();
        sil_sharp.get_shader().texture("tex", atl_tmp.tex(0));
        sil_sharp.get_shader().uniform("tex_transform", float4(0, 0, 1, 1));
        sil_sharp.get_shader().uniform("layer", (float)l);
        sil_sharp.get_shader().uniform("threshold", 0.05f);
        sil_sharp.render();
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    //engine::textureManager->save_png(atl.tex(0), "atlass_gauss_3");
}
