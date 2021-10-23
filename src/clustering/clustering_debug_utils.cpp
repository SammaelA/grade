#include "clustering_debug_utils.h"
#include "../tinyEngine/utility/postfx.h"
#include "../texture_manager.h"
#include "../grove_packer.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include "../tinyEngine/save_utils/csv.h"
#include "clustering_debug_status.h"
#include "hasing.h"

ClusteringDebugInfo clusteringDebugInfo;
using namespace glm;

void visualize_clusters(Block &settings, std::vector<BranchClusteringData *> branches, 
                        std::vector<ClusteringBase::ClusterStruct> &clusters,
                        ClusteringContext *ctx, std::string file_name, int w, int h)
{
    ctx->self_impostors_raw_atlas = new TextureAtlasRawData(ctx->self_impostors_data->atlas);
    int max_size = 0;
    std::vector<std::pair<glm::ivec4, Billboard *>> icons;
    std::vector<ivec4> cluster_frames;
    std::vector<int> cluster_nums;
    int tex_w = 0, tex_h = 0;
    for (auto &cs : clusters)
    {
        max_size = MAX(max_size, cs.members.size());
    }
    tex_w = w * ceil(sqrt(max_size));
    int cur_h = 0;
    int cur_w = 0;
    int layer_h = 0;
    int n = 0;
    for (auto &cs : clusters)
    {
        std::vector<Billboard *> child_clusters;
        for (auto &p : cs.members)
        {
            BranchClusteringDataImpostor *imp = dynamic_cast<BranchClusteringDataImpostor *>(branches[p.first]);
            
            int pos = round(imp->self_impostor->slices.size()*p.second.rot/(2*PI));
            Billboard *icon = &(imp->self_impostor->slices[pos]);
            child_clusters.push_back(icon);
        }

        int sz = child_clusters.size();
        int icons_w = ceil(sqrt(sz));
        int cluster_w = w * icons_w;
        int icons_h = ceil(sz / icons_w);
        int cluster_h = h * icons_h;
        if (cur_w + cluster_w > tex_w)
        {
            cur_h += layer_h + 0.4 * w;
            layer_h = 0;
            cur_w = 0;
        }
        cluster_frames.push_back(glm::ivec4(cur_w, cur_h, cluster_w, cluster_h));
        cluster_nums.push_back(n);
        for (int i = 0; i < icons_w; i++)
        {
            for (int j = 0; j < icons_h; j++)
            {
                int n = i * icons_h + j;
                if (n < child_clusters.size())
                {
                    glm::ivec4 bounds = glm::ivec4(cur_w + w * i, cur_h + h * j, w, h);
                    icons.push_back(std::pair<glm::ivec4, Billboard *>(bounds, child_clusters[n]));
                }
            }
        }
        layer_h = MAX(layer_h, cluster_h);
        cur_w += cluster_w;
        n++;
    }
    tex_h = cur_h + layer_h;

    Texture res(textureManager.create_unnamed(tex_w, tex_h));

    if (false)
    {
        PostFx copy = PostFx("copy_arr2.fs");
        PostFx frame = PostFx("thick_frame.fs");
        GLuint fbo;
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, res.texture, 0);
        glViewport(0, 0, tex_w, tex_h);
        glClearColor(0.8, 0.8, 1, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        copy.use();
        copy.get_shader().texture("tex", ctx->self_impostors_data->atlas.tex(0));

        for (auto &p : icons)
        {

            auto &bill = *(p.second);
            glm::vec3 tc_from = glm::vec3(0, 0, 0);
            glm::vec3 tc_to = glm::vec3(1, 1, 0);
            ctx->self_impostors_data->atlas.process_tc(bill.id, tc_from);
            ctx->self_impostors_data->atlas.process_tc(bill.id, tc_to);
            glm::vec4 transform = glm::vec4(tc_from.x, tc_from.y, tc_to.x - tc_from.x, tc_to.y - tc_from.y);
            copy.get_shader().uniform("tex_transform", transform);
            copy.get_shader().uniform("layer", tc_from.z);
            logerr("%f %f %f %f %f transform", transform.x, transform.y, transform.z, transform.w, tc_from.z);
            glViewport(p.first.x, p.first.y, p.first.z, p.first.w);
            copy.get_shader().texture("tex", ctx->self_impostors_data->atlas.tex(0));

            copy.render();
        }

        frame.use();
        frame.get_shader().uniform("thickness", 0.06f);
        frame.get_shader().uniform("color", glm::vec4(0, 0, 0, 1));
        for (auto &p : icons)
        {
            glViewport(p.first.x, p.first.y, p.first.z, p.first.w);
            logerr("%d %d %d %d", p.first.x, p.first.y, p.first.z, p.first.w);
            frame.render();
        }

        frame.get_shader().uniform("thickness", 0.06f);
        frame.get_shader().uniform("color", glm::vec4(1, 0, 0, 1));
        for (auto &cl : cluster_frames)
        {
            glViewport(cl.x, cl.y, cl.z, cl.w);
            frame.render();
        }
    }
    unsigned char *data = new unsigned char[4 * tex_w * tex_h];
    for (int i = 0; i < 4 * tex_w * tex_h; i += 4)
    {
        data[i] = 0;
        data[i + 1] = 0;
        data[i + 2] = 0;
        data[i + 3] = 255;
    }

    for (auto &p : icons)
    {
        if (!p.second)
        {
            continue;
        }
        auto &bill = *(p.second);
        glm::ivec4 sizes = ctx->self_impostors_data->atlas.get_sizes();
        sizes.x /= sizes.z;
        sizes.y /= sizes.w;
        for (int i = 0; i < sizes.x; i++)
        {
            for (int j = 0; j < sizes.y; j++)
            {
                if (((p.first.x + i) * sizes.y + (p.first.y + j)) > tex_w * tex_h)
                    logerr("(%d %d) %d %d %d %d %d %d", sizes.x, sizes.y, tex_w, tex_h, p.first.y, j, p.first.x, i);
                if (ctx->self_impostors_raw_atlas->get_pixel_uc(i, j, Channel::A, bill.id) > 10)
                {
                    data[4 * ((p.first.y + j) * tex_w + (p.first.x + i))] =
                        ctx->self_impostors_raw_atlas->get_pixel_uc(i, j, Channel::R, bill.id);
                    data[4 * ((p.first.y + j) * tex_w + (p.first.x + i)) + 1] =
                        ctx->self_impostors_raw_atlas->get_pixel_uc(i, j, Channel::G, bill.id);
                    data[4 * ((p.first.y + j) * tex_w + (p.first.x + i)) + 2] =
                        ctx->self_impostors_raw_atlas->get_pixel_uc(i, j, Channel::B, bill.id);
                    data[4 * ((p.first.y + j) * tex_w + (p.first.x + i)) + 3] = 255;
                }
            }
        }
    }
    int thickness = 4;
    glm::ivec4 base_color = glm::ivec4(255,150,150,255);
    glm::ivec4 color = base_color;
    std::vector<glm::vec3> colors =
    {
        glm::vec3(0,0,0), glm::vec3(0,0,1), glm::vec3(0,1,0), glm::vec3(0,1,1), 
        glm::vec3(1,0,0), glm::vec3(1,0,1), glm::vec3(1,1,0), glm::vec3(0,0.5,0), 
        glm::vec3(0,0,0.5)
    };
    for (int i=0;i<cluster_frames.size();i++)
    {
        auto &cl = cluster_frames[i];
        color = i >= colors.size() ? base_color : glm::ivec4(255*colors[i].x,255*colors[i].y,255*colors[i].z,255);
        #define FRAME(i,j) \
        data[4 * ((cl.y + j) * tex_w + (cl.x + i))] = color.x;\
        data[4 * ((cl.y + j) * tex_w + (cl.x + i))+1] = color.y;\
        data[4 * ((cl.y + j) * tex_w + (cl.x + i))+2] = color.z;\
        data[4 * ((cl.y + j) * tex_w + (cl.x + i))+3] = color.w;

        for (int i=0;i<cl.z;i++)
        {
            for (int j=0;j<thickness;j++)
            {
                FRAME(i,j)
            }
        }
        for (int i=0;i<cl.z;i++)
        {
            for (int j=cl.w-thickness;j<cl.w;j++)
            {
                FRAME(i,j)
            }
        }
        for (int i=0;i<thickness;i++)
        {
            for (int j=0;j<cl.w;j++)
            {
                FRAME(i,j)
            }
        }
        for (int i=cl.z-thickness;i<cl.z;i++)
        {
            for (int j=0;j<cl.w;j++)
            {
                FRAME(i,j)
            }
        }
    }

    textureManager.save_bmp_raw(data, tex_w, tex_h, 4, file_name);
    textureManager.delete_tex(res);
    delete data;

    if (ctx->self_impostors_raw_atlas)
    {
        delete ctx->self_impostors_raw_atlas;
    }
}

void save_csv_impostor(std::string &save_path, ClusteringContext *ctx, std::vector<ClusterPackingLayer> &packingLayers)
{

    int cnt = 0;
    TextureAtlasRawData raw_atlas = TextureAtlasRawData(ctx->self_impostors_data->atlas);
    unsigned char *sl_data = safe_new<unsigned char>(2 * raw_atlas.get_slice_size(0), "sl_data");
    memset(sl_data, 0, 2 * raw_atlas.get_slice_size(0));
    unsigned char *tmp_data = safe_new<unsigned char>(raw_atlas.get_slice_size(0), "tmp_data");
    int ww = 0, hh = 0;
    int sl = 0;

    //find largest branch to rescale dataset images properly
    glm::vec3 max_sizes = glm::vec3(0, 0, 0);
    for (int i = 0; i < packingLayers.size(); i++)
    {
        for (auto &c : packingLayers[i].clusters)
        {
            for (auto *cd : c.ACDA.clustering_data)
            {
                auto *imp_cd = dynamic_cast<BranchClusteringDataImpostor *>(cd);
                if (imp_cd)
                {
                    max_sizes = max(max_sizes, imp_cd->min_bbox.sizes);
                }
            }
        }
    }
    float q = 0.85;
    float max_size = q * MAX(max_sizes.x, MAX(max_sizes.y, max_sizes.z));
    bool need_rescale = true;
    try
    {
        raw_atlas.get_slice(0, sl_data, &ww, &hh);
        auto *base_imp = dynamic_cast<BranchClusteringDataImpostor *>(packingLayers[0].clusters[0].ACDA.clustering_data[0]);
        std::vector<std::string> columns;
        for (int i=0;i<base_imp->self_impostor->slices.size();i++)
        {
            for (int y = 0; y < hh; y++)
            {
                for (int x = 0; x < ww; x++)
                {
                    columns.push_back(std::to_string(i) + "_pixel_" + std::to_string(y) + "_" + std::to_string(x) + "_red");
                    columns.push_back(std::to_string(i) + "_pixel_" + std::to_string(y) + "_" + std::to_string(x) + "_green");
                }
            }
        }
        columns.push_back("target");
        CSVData table = CSVData(columns);
        for (int i = 0; i < packingLayers.size(); i++)
        {
            for (auto &c : packingLayers[i].clusters)
            {
                for (auto *cd : c.ACDA.clustering_data)
                {
                    auto *imp_cd = dynamic_cast<BranchClusteringDataImpostor *>(cd);
                    if (imp_cd)
                    {
                        std::vector<int> row;
                        for (auto &bill : imp_cd->self_impostor->slices)
                        {
                            raw_atlas.get_slice(bill.id, sl_data, &ww, &hh);

                            float sz = MAX(imp_cd->min_bbox.sizes.x, MAX(imp_cd->min_bbox.sizes.y, imp_cd->min_bbox.sizes.z));
                            float scale = max_size / sz;

                            for (int y = 0; y < hh; y++)
                            {
                                float y_src = scale * y;
                                int y0 = y_src;
                                float qy = y_src - y0;
                                for (int x = 0; x < ww; x++)
                                {
                                    float x_src = scale * x;
                                    if (x_src > ww + 1 || y_src > hh + 1)
                                    {
                                        for (int ch = 0; ch < 3; ch++)
                                            tmp_data[4 * (y * ww + x) + ch] = 0;
                                    }
                                    else
                                    {
                                        int x0 = x_src;
                                        float qx = x_src - x0;
                                        for (int ch = 0; ch < 3; ch++)
                                        {
                                            tmp_data[4 * (y * ww + x) + ch] =
                                                (1 - qy) * ((1 - qx) * sl_data[4 * (y0 * ww + x0) + ch] +
                                                            qx * sl_data[4 * (y0 * ww + x0 + 1) + ch]) +
                                                qy * ((1 - qx) * sl_data[4 * ((y0 + 1) * ww + x0) + ch] +
                                                      qx * sl_data[4 * ((y0 + 1) * ww + x0 + 1) + ch]);
                                            tmp_data[4 * (y * ww + x) + ch] = sl_data[4 * (y0 * ww + x0) + ch];
                                        }
                                    }
                                    tmp_data[4 * (y * ww + x) + 3] = 255;
                                }
                            }

                            for (int y = 0; y < hh; y++)
                            {
                                for (int x = 0; x < ww; x++)
                                {
                                    row.push_back(tmp_data[4 * (y * ww + x)]);
                                    row.push_back(tmp_data[4 * (y * ww + x) + 1]);
                                }
                            }

                            sl++;
                        }
                        row.push_back(cnt);
                        logerr("%d cnt", cnt);
                        table.add_row(row, 0);
                    }
                    else
                    {
                        logerr("error - trying to save impostors csv without impostor. Check clustering settings");
                    }
                }
                cnt++;
            }
        }
        CSVSaver::save_csv_in_file(table, clusteringDebugInfo.csv_file_name);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }

    safe_delete<unsigned char>(sl_data, "sl_data");
    safe_delete<unsigned char>(tmp_data, "tmp_data");
    raw_atlas.clear();
}

void save_csv_hash(std::string &save_path, ClusteringContext *ctx, std::vector<ClusterPackingLayer> &packingLayers)
{
    logerr("adsdsdsds");
    int cnt = 0;
    auto *h0 = dynamic_cast<BranchHash *>(packingLayers[0].clusters[0].ACDA.clustering_data[0]);
    std::vector<std::string> columns;
    for (int y = 0; y < h0->hashes[0].data.size(); y++)
    {
        columns.push_back("hash_" + std::to_string(y));
    }
    columns.push_back("target");
    CSVData table = CSVData(columns);
    for (int i = 0; i < packingLayers.size(); i++)
    {
        for (auto &c : packingLayers[i].clusters)
        {
            for (auto *cd : c.ACDA.clustering_data)
            {
                auto *hash_cd = dynamic_cast<BranchHash *>(cd);
                if (hash_cd)
                {
                    for (int hash_n = 0; hash_n < hash_cd->hashes.size(); hash_n++)
                    {

                        std::vector<float> row = hash_cd->hashes[hash_n].data;

                        row.push_back(cnt);
                        logerr("%d cnt %f %f %f %f %f", cnt, row[0],row[1],row[2],row[3],row[4]);
                        if (hash_n == 0)
                            table.add_row(row, 0);
                    }
                }
                else
                {
                    logerr("error - trying to save hashes csv without hash. Check clustering settings");
                }
            }
            cnt++;
        }
    }
    CSVSaver::save_csv_in_file(table, clusteringDebugInfo.csv_file_name);
}

void save_csv(std::string &save_path, ClusteringContext *ctx, std::vector<ClusterPackingLayer> &packingLayers)
{
    if (dynamic_cast<BranchClusteringDataImpostor *>(packingLayers[0].clusters[0].ACDA.clustering_data[0]))
        save_csv_impostor(save_path, ctx, packingLayers); 
    else if (dynamic_cast<BranchHash *>(packingLayers[0].clusters[0].ACDA.clustering_data[0]))
        save_csv_hash(save_path, ctx, packingLayers);
    else
        logerr("cannot save branch clustering data in csv file - it is not impostors or hashes");
}

void prepare_dataset(std::string &save_path, ClusteringContext *ctx, std::vector<ClusterPackingLayer> &packingLayers)
{

    int cnt = 0;
    TextureAtlasRawData raw_atlas = TextureAtlasRawData(ctx->self_impostors_data->atlas);
    unsigned char *sl_data = safe_new<unsigned char>(2*raw_atlas.get_slice_size(0), "sl_data");
    memset(sl_data,0,2*raw_atlas.get_slice_size(0));
    unsigned char *tmp_data = safe_new<unsigned char>(raw_atlas.get_slice_size(0), "tmp_data");
    int ww = 0, hh = 0;
    bool info_files = true;
    
    std::string database_file;
    std::string test_file;
    std::string train_file;
    
    float train_part = 0.9;
    float test_part = 0.05;
    
    std::string database_name = "database.txt";
    std::string test_name = "test.txt";
    std::string train_name = "train.txt";
    
    std::string folder_name = "images";

    int sl = 0;
    int train_elems = 0;
    int test_elems = 0;

    std::string dir_path;

    //find largest branch to rescale dataset images properly
    glm::vec3 max_sizes = glm::vec3(0,0,0);
    for (int i = 0; i< packingLayers.size();i++)
    {
        for (auto &c : packingLayers[i].clusters)
        {
            for (auto *cd : c.ACDA.clustering_data)
            {
                auto *imp_cd = dynamic_cast<BranchClusteringDataImpostor *>(cd);
                if (imp_cd)
                {
                    max_sizes = max(max_sizes, imp_cd->min_bbox.sizes);
                }
            }
        }
    }
    float q = 0.85;
    float max_size = q*MAX(max_sizes.x,MAX(max_sizes.y,max_sizes.z));
    bool need_rescale = true;
    try
    {
        int clusters_count = 0;
        if (info_files)
        {
            dir_path = save_path + "/"+folder_name;
            boost::filesystem::create_directory(dir_path);
            boost::filesystem::permissions(dir_path, boost::filesystem::perms::all_all); 

            for (int i = 0; i< packingLayers.size();i++)
            {
                clusters_count += packingLayers[i].clusters.size();
            }         
        }
        raw_atlas.get_slice(0, sl_data, &ww, &hh);
        std::vector<std::string> columns;
        for (int y = 0; y < hh; y++)
        {
            for (int x = 0; x < ww; x++)
            {
                columns.push_back("pixel_" + std::to_string(y) + "_" + std::to_string(x) + "_red");
                columns.push_back("pixel_" + std::to_string(y) + "_" + std::to_string(x) + "_green");
            }
        }
        columns.push_back("target");
        CSVData table = CSVData(columns);
        for (int i = 0; i< packingLayers.size();i++)
        {
            for (auto &c : packingLayers[i].clusters)
            {
                std::string cluster_labels;
                if (info_files)
                {
                    for (int j=0;j<clusters_count;j++)
                    {
                        cluster_labels += ((j == cnt) ? " 1" : " 0");
                    }
                }
                else
                {
                    dir_path = save_path + "/" + std::to_string(cnt);
                    boost::filesystem::create_directory(dir_path);
                    boost::filesystem::permissions(dir_path, boost::filesystem::perms::all_all);
                    sl = 0;
                }
                for (auto *cd : c.ACDA.clustering_data)
                {
                    auto *imp_cd = dynamic_cast<BranchClusteringDataImpostor *>(cd);
                    if (imp_cd)
                    {
                        for (auto &bill : imp_cd->self_impostor->slices)
                        {
                            std::string file_path = dir_path + "/" + std::to_string(sl)+".bmp";
                            raw_atlas.get_slice(bill.id, sl_data, &ww, &hh);

                            if (need_rescale)
                            {
                                float sz = MAX(imp_cd->min_bbox.sizes.x, MAX(imp_cd->min_bbox.sizes.y, imp_cd->min_bbox.sizes.z));
                                float scale = max_size / sz;

                                for (int y = 0; y < hh; y++)
                                {
                                    float y_src = scale * y;
                                    int y0 = y_src;
                                    float qy = y_src - y0;
                                    for (int x = 0; x < ww; x++)
                                    {
                                        float x_src = scale * x;
                                        if (x_src > ww + 1 || y_src > hh + 1)
                                        {
                                            for (int ch = 0; ch < 3; ch++)
                                                tmp_data[4 * (y * ww + x) + ch] = 0;
                                        }
                                        else
                                        {
                                            int x0 = x_src;
                                            float qx = x_src - x0;
                                            for (int ch = 0; ch < 3; ch++)
                                            {
                                                tmp_data[4 * (y * ww + x) + ch] =
                                                    (1 - qy) * ((1 - qx) * sl_data[4 * (y0 * ww + x0) + ch] +
                                                                qx * sl_data[4 * (y0 * ww + x0 + 1) + ch]) +
                                                    qy * ((1 - qx) * sl_data[4 * ((y0 + 1) * ww + x0) + ch] +
                                                          qx * sl_data[4 * ((y0 + 1) * ww + x0 + 1) + ch]);
                                                tmp_data[4 * (y * ww + x) + ch] = sl_data[4 * (y0 * ww + x0) + ch];
                                            }
                                        }
                                        tmp_data[4 * (y * ww + x) + 3] = 255;
                                    }
                                }
                                textureManager.save_bmp_raw_directly(tmp_data, ww, hh, 4, file_path);
                            }
                            else
                                textureManager.save_bmp_raw_directly(sl_data, ww, hh, 4, file_path);
                            
                            std::vector<int> row;
                            for (int y = 0; y < hh; y++)
                            {
                                for (int x = 0; x < ww; x++)
                                {
                                    row.push_back(tmp_data[4*(y*ww + x)]);
                                    row.push_back(tmp_data[4*(y*ww + x) + 1]);
                                }
                            }
                            row.push_back(cnt);
                            logerr("%d cnt",cnt);
                            if (bill.id == imp_cd->self_impostor->slices.front().id)
                                table.add_row(row,0);
                            if (info_files)
                            {
                                //add a record about images to database and (maybe) test or train lists;
                                std::string record = folder_name + "/" + std::to_string(sl)+".bmp" + cluster_labels + "\n";
                                database_file += record;
                                if (urand() < train_part)
                                {
                                    train_file += record;
                                    train_elems++;
                                }
                                if (urand() < test_part)
                                {
                                    test_file += record;
                                    test_elems++;
                                }
                            }
                            sl++;
                        }
                    }
                    else
                    {
                        logerr("error - trying to save to dataset branch without impostor. Check clustering settings");
                    }
                }
                cnt++;
                if (info_files)
                {
                    std::ofstream database_ofs;
                    database_ofs.open(save_path+"/database.txt");
                    database_ofs << database_file;
                    database_ofs.close();

                    std::ofstream test_ofs;
                    test_ofs.open(save_path+"/test.txt");
                    test_ofs << test_file;
                    test_ofs.close();

                    std::ofstream train_ofs;
                    train_ofs.open(save_path+"/train.txt");
                    train_ofs << train_file;
                    train_ofs.close();

                    std::ofstream info_ofs;
                    info_ofs.open(save_path+"/info.txt");
                    info_ofs << std::string(std::to_string(clusters_count)+" "+std::to_string(sl)+" "+
                                            std::to_string(train_elems)+" "+std::to_string(test_elems));
                    info_ofs.close();

                }
            }
        }
        CSVSaver::save_csv_in_file(table, "scripts/test.csv");
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    safe_delete<unsigned char>(sl_data, "sl_data");
    safe_delete<unsigned char>(tmp_data, "tmp_data");
    raw_atlas.clear();
}