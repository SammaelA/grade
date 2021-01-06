#include "grove.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "billboard_cloud.h"

GroveRenderer::GroveRenderer(GrovePacked *_source, int LODs_count) : renderer({"billboard_render.vs", "billboard_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
                                                                     rendererInstancing({"billboard_render_instancing.vs", "billboard_render_instancing.fs"},
                                                                                        {"in_Position", "in_Normal", "in_Tex", "in_Model"}),
                                                                     leaf(textureManager.get("leaf")),
                                                                     wood(textureManager.get("wood"))
{
    Visualizer v = Visualizer();
    source = _source;
    debug("creating grove renderer with %d LODs\n", _source->clouds.size());
    for (int i = 0; i < _source->clouds.size(); i++)
    {
        LODs.emplace_back();
        LODs.back().cloud = _source->clouds[i];
        LODs.back().cloud->set_render_mode(BillboardCloud::ONLY_INSTANCES);
        Model *m = new Model();
        for (int j = 0; j < i; j++)
        {
            auto packed_branches = _source->uniqueCatalogue.get_level(j);
            for (PackedBranch &pb : packed_branches)
            {
                v.packed_branch_to_model(pb, m, false);
            }
        }
        LODs.back().m = m;
    }

    LODs.emplace_back();
    LODs.back().cloud = nullptr;
    Model *m = new Model();
    for (int j = 0; j < source->uniqueCatalogue.levels(); j++)
    {
        auto packed_branches = source->uniqueCatalogue.get_level(j);
        for (PackedBranch &pb : packed_branches)
        {
            v.packed_branch_to_model(pb, m, false);
        }
    }
    LODs.back().m = m;
    for (InstancedBranch &b : source->instancedBranches)
    {
        add_instance_model(LODs.back(), source, b);
    }
}
void GroveRenderer::render(int lod, glm::mat4 prc)
{
    if (lod < 0 || lod >= LODs.size())
    {
        //logerr("trying to render grove with wrong LOD number %d. Grove has %d LODs",lod, LODs.size());
        //return;
        lod = LODs.size() - 1;
    }
    Model *m = LODs[lod].m;
    renderer.use();
    renderer.texture("tex", wood);
    renderer.uniform("projectionCamera", prc);
    renderer.uniform("model", m->model);
    m->update();
    m->render(GL_TRIANGLES);
    if (LODs[lod].cloud)
        LODs[lod].cloud->render(prc);

    rendererInstancing.use();
    rendererInstancing.texture("tex", wood);
    rendererInstancing.uniform("projectionCamera", prc);
    for (Instance *in : LODs[lod].instances)
    {
        Model *m = (Model *)(in->m);
        m->update();
        in->render(GL_TRIANGLES);
    }
}
void GroveRenderer::add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch)
{
    Visualizer v = Visualizer();
    Model *m = new Model();
    for (int id : branch.branches)
    {
        v.packed_branch_to_model(source->instancedCatalogue.get(id), m, false);
    }
    m->update();
    Instance *in = new Instance(m);
    debugl(4, "transforms size %d\n", branch.transforms.size());
    in->addBufferCopy(branch.transforms);
    lod.instances.push_back(in);
}