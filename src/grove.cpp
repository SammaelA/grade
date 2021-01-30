#include "grove.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "billboard_cloud.h"

GroveRenderer::GroveRenderer(): 
renderer({"simple_render.vs", "simple_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
rendererInstancing({"simple_render_instancing.vs", "simple_render_instancing.fs"},
                   {"in_Position", "in_Normal", "in_Tex", "in_Model"}),
leaf(textureManager.get("leaf")),
wood(textureManager.get("wood"))
{

}
GroveRenderer::GroveRenderer(GrovePacked *_source, int LODs_count) :
GroveRenderer()
{
    Visualizer v = Visualizer();
    source = _source;
    debug("creating grove renderer with %d LODs\n", _source->clouds.size());
    for (int i = 0; i < _source->clouds.size(); i++)
    {
        LODs.emplace_back();
        LODs.back().cloud = new BillboardCloudRenderer(&_source->clouds[i]);
        LODs.back().cloud->set_render_mode(BillboardCloudRenderer::ONLY_INSTANCES);
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
        for (InstancedBranch &b : source->instancedBranches)
        {
            add_instance_model(LODs.back(), source, b, i-1);
        }
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
        add_instance_model(LODs.back(), source, b,1000);
    }
}
GroveRenderer::~GroveRenderer()
{
    for (int i=0;i<LODs.size();i++)
    {
        if (LODs[i].m)
            delete LODs[i].m;
        if (LODs[i].cloud)
            delete LODs[i].cloud;
        for (int j=0;j<LODs[i].instances.size();j++)
        {
            delete LODs[i].instances[j]->m;
            delete LODs[i].instances[j];
        }
    }
    LODs.clear();
    source = nullptr;
}
void GroveRenderer::render(int lod, glm::mat4 prc)
{
    if (LODs.size()==0)
        return;
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
void GroveRenderer::add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch, int up_to_level)
{
    Visualizer v = Visualizer();
    Model *m = new Model();
    for (int id : branch.branches)
    {
        PackedBranch &b = source->instancedCatalogue.get(id);
        if (b.level <= up_to_level)
            v.packed_branch_to_model(b, m, false);
    }
    m->update();
    Instance *in = new Instance(m);
    in->addBufferCopy(branch.transforms);
    lod.instances.push_back(in);
}