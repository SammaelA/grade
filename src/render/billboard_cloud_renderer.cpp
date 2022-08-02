#include "billboard_cloud_renderer.h"
#include "graphics_utils/texture_manager.h"

using namespace glm;

BillboardCloudRenderer::BillboardCloudRenderer(const BillboardCloudData *_data):
Countable(3),
rendererToTexture({"render_to_billboard.vs", "render_to_billboard.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRenderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
billboardRendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                            {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"})
{
    this->data = (BillboardCloudData *)_data;
    if (!data || !data->valid)
    {
        return;
    }
    for (BillboardData &bill : data->billboards)
    {
        Model *m = new Model();
        for (Billboard &b : bill.billboards)
        {
            b.instancing = true;
            b.to_model(m, data->atlas);
        }
        for (int i=0;i<bill.IDA.centers_par.size();i++)
        {
            if (data->level > GroveRenderer::base_level)
                bill.IDA.centers_par[i] = bill.IDA.centers_self[i];
        }

        instances.push_back(m);
    }
}
BillboardCloudRenderer::~BillboardCloudRenderer()
{
    debugl(11,"deleting BCR\n");
    if (cloud)
        delete (cloud);
    debugl(11,"cloud deleted\n");
    for (int i=0;i<instances.size();i++)
    {
        delete instances[i];
       debugl(11,"instance %d deleted\n",i);
    }
}
void BillboardCloudRenderer::render(MultiDrawRendDesc &mdrd, glm::mat4 &projection, glm::mat4 &view, DirectedLight &light, 
                                    glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                                    glm::vec4 screen_size, bool to_shadow, GroveRendererDebugParams dbgpar)
{
    if (to_shadow)
        return;
    if (!data || !data->valid)
        return;
    
    if (renderMode == ONLY_SINGLE || renderMode == BOTH)
    {
        std::function<void(Model *)> _ce = [&](Model *h) {
            for (BillboardData &bill : data->billboards)
            {
                for (Billboard &b : bill.billboards)
                    b.to_model(h, data->atlas);
            }
        };
        cloud->construct(_ce);
        billboardRenderer.use();
        billboardRenderer.texture("tex", data->atlas.tex(0));
        billboardRenderer.uniform("model", cloud->model);
        billboardRenderer.uniform("projectionCamera", projection *view);
     
        cloud->render(GL_TRIANGLES);
    }
    if (renderMode == ONLY_INSTANCES || renderMode == BOTH)
    {
        billboardRendererInstancing.use();
        billboardRendererInstancing.uniform("camera_pos", camera_pos);
        billboardRendererInstancing.uniform("screen_size",screen_size);
        billboardRendererInstancing.texture("color_tex", data->atlas.tex(0));
        billboardRendererInstancing.texture("normal_tex", data->atlas.tex(1));
        billboardRendererInstancing.texture("noise",textureManager.get("noise"));
        billboardRendererInstancing.uniform("projection", projection);
        billboardRendererInstancing.uniform("view", view);
        billboardRendererInstancing.uniform("type_id", (uint)mdrd.type_id);
        billboardRendererInstancing.uniform("debug_model_id",dbgpar.need_focus_model ? dbgpar.model_focused : -1);
        
        glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                            mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
    }
}