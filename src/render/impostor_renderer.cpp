#include "graphics_utils/impostor.h"
#include "tinyEngine/engine.h"
#include "impostor_renderer.h"
#include "rendering_SSBO_structs.h"
#include "grove_renderer.h"

using glm::vec3;
using glm::vec4;
using glm::mat4;

ImpostorRenderer::~ImpostorRenderer()
{
    #define DELBUF(a) if (a) { delete_buffer((a)); a = 0;}

    DELBUF(slicesBuffer);
    DELBUF(impostorsDataBuffer);

    for (int i=0;i<models.size();i++)
    {
        if (models[i])
        {
            delete models[i];
            models[i] = nullptr;
        }
    }
}
ImpostorRenderer::ImpostorRenderer(ImpostorsData *data):
impostorRenderer({"impostor_render.vs", "impostor_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
impostorRendererInstancing({"impostor_render_instancing.vs", "impostor_render_instancing.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{
    if (!data || !data->valid)
    {
        logerr("empty or corrupted impostor data");
        return;
    }
    this->data = data;
    if (data->impostors.size() == 1)
    {
        hth = 0.1;
        vth = 0.5;
    }
    else
    {
        hth = 0;
        vth = 0.1;
    }
    std::vector<float> s_verts;
    for (Impostor &imp : data->impostors)
    {
        offsets.push_back(s_verts.size()/8);
        std::vector<glm::vec3> tcs = imp.top_slice.get_tc(data->atlas);
        for (int i = 0;i<MIN(imp.top_slice.positions.size(),tcs.size());i++)
        {
            vec3 &pos = imp.top_slice.positions[i];
            
            s_verts.push_back(pos.x);
            s_verts.push_back(pos.y);
            s_verts.push_back(pos.z);
            s_verts.push_back(1);

            s_verts.push_back(tcs[i].x);
            s_verts.push_back(tcs[i].y);
            s_verts.push_back(tcs[i].z);
            s_verts.push_back(1);
        }
        for (Billboard &b : imp.slices)
        {
            tcs = b.get_tc(data->atlas);
            for (int i = 0;i<MIN(b.positions.size(),tcs.size());i++)
            {
                vec3 &pos = b.positions[i];

                s_verts.push_back(pos.x);
                s_verts.push_back(pos.y);
                s_verts.push_back(pos.z);
                s_verts.push_back(1);
                            
                s_verts.push_back(tcs[i].x);
                s_verts.push_back(tcs[i].y);
                s_verts.push_back(tcs[i].z);
                s_verts.push_back(1);
            }
        }
        Model *bm = new Model();
        std::vector<float> vertexes = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
        std::vector<float> tc = {0,0,0,0, 1,0,0,0, 0,1,0,0, 1,1,0,0};
        std::vector<GLuint> indices = {0, 1, 3, 2, 0, 3};

        std::function<void(Model *)> _c_mip = [&](Model *h) 
        {
            bm->positions = vertexes;
            bm->colors = tc;
            bm->indices = indices;
        };
        bm->construct(_c_mip);
        models.push_back(bm);
    }
    std::vector<ImpostorData> imp_data_buffer;
    auto it = data->impostors.begin();
    for (int i=0;i<models.size();i++)
    {
        ImpostorData dat;
        dat.slice_count = (it->slices.size());
        dat.slice_offset = offsets[i];
        dat.slice_verts = it->slices[0].positions.size();
        dat.imp_center = glm::vec4(it->bcyl.center,1);
        imp_data_buffer.push_back(dat);
        it++;
    }
    impostorsDataBuffer = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, impostorsDataBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(ImpostorData)*imp_data_buffer.size(), imp_data_buffer.data(), GL_STATIC_DRAW);
    slicesBuffer = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slicesBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*s_verts.size(), s_verts.data(), GL_STATIC_DRAW);


}
void ImpostorRenderer::render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light,
                              glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                              glm::vec4 screen_size, bool to_shadow, GroveRendererDebugParams dbgpar)
{
    if (to_shadow)//we do not render impostors to shadow
        return;
    if (!data || !data->valid)
        return;
    
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, slicesBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, impostorsDataBuffer);

    impostorRendererInstancing.use();
    impostorRendererInstancing.uniform("projection", projection);
    impostorRendererInstancing.uniform("view", view);
    impostorRendererInstancing.texture("color_tex", data->atlas.tex(0));
    impostorRendererInstancing.texture("normal_tex", data->atlas.tex(1));
    impostorRendererInstancing.uniform("camera_pos", camera_pos);
    impostorRendererInstancing.uniform("screen_size", screen_size);
    impostorRendererInstancing.texture("noise", engine::textureManager->get("noise"));
    impostorRendererInstancing.uniform("hor_vert_transition_thr", glm::vec2(hth,vth));
    impostorRendererInstancing.uniform("delta",0.5f);
    impostorRendererInstancing.uniform("type_id",(uint)mdrd.type_id);
    impostorRendererInstancing.uniform("vertex_id_offset",mdrd.base_vertex_id);
    impostorRendererInstancing.uniform("debug_model_id",dbgpar.need_focus_model ? dbgpar.model_focused : -1);
    glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                        mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
}