#include "visualize_clusters.h"
#include "../tinyEngine/utility/postfx.h"
#include "../texture_manager.h"

using namespace glm;
void visualize_clusters(Block &settings, std::vector<BranchClusteringData *> branches, 
                        std::vector<ClusteringBase::ClusterStruct> &clusters,
                        ClusteringContext *ctx, std::string file_name, int w, int h)
{
    ctx->self_impostors_raw_atlas = new TextureAtlasRawData(ctx->self_impostors_data->atlas);
    int max_size = 0;
    std::vector<std::pair<glm::ivec4, Billboard *>> icons;
    std::vector<ivec4> cluster_frames;
    int tex_w = 0, tex_h = 0;
    for (auto &cs : clusters)
    {
        max_size = MAX(max_size, cs.members.size());
    }
    tex_w = w * ceil(sqrt(max_size));
    int cur_h = 0;
    int cur_w = 0;
    int layer_h = 0;
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
    glm::ivec4 color = glm::ivec4(255,150,150,255);
    for (auto &cl : cluster_frames)
    {
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