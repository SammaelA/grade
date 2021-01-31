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
GroveRenderer::GroveRenderer(GrovePacked *_source, int LODs_count, std::vector<float> &max_distances) :
GroveRenderer()
{
    if (LODs_count != _source->clouds.size() + 1|| max_distances.size() != _source->clouds.size() + 1)
    {
        logerr("Can not calculate LODs count for GroveRenderer. Given LODs_count doesn't match given data");
        LODs_count = _source->clouds.size();
        for (int i=max_distances.size(); i<_source->clouds.size();i++)
        {
            max_distances.push_back(2*max_distances.back());
        }
    }
    Visualizer v = Visualizer();
    source = _source;
    debug("creating grove renderer with %d LODs\n", _source->clouds.size());
    for (int i = 0; i < _source->clouds.size(); i++)
    {
        LODs.emplace_back();
        LODs.back().max_dist = max_distances[i];
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
    LODs.back().max_dist = max_distances[LODs_count - 1];
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
        add_instance_model(LODs.back(), source, b,1000,true);
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
void GroveRenderer::render_auto_LOD(glm::mat4 prc, glm::vec3 camera_pos)
{
    int lod = -1;
    float len = glm::length(source->center - camera_pos);
    for (int i=1;i<LODs.size();i++)
    {
        if (LODs[i].max_dist < len && LODs[i-1].max_dist >= len)
        {
            lod = i - 1;
            break;
        }
    }
    if (LODs.back().max_dist > len)
        lod = LODs.size() - 1;
    if (lod == -1)
        return;
    else
        render(lod,prc,camera_pos);
    
}
void GroveRenderer::render(int lod, glm::mat4 prc, glm::vec3 camera_pos)
{
    if (LODs.size()==0)
        return;
    if (lod == -1)
    {
        render_auto_LOD(prc,camera_pos);
        return;
    }
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

    rendererInstancing.texture("tex", leaf);
    for (Instance *in : LODs[lod].leaves_instances)
    {
        Model *m = (Model *)(in->m);
        m->update();
        in->render(GL_TRIANGLES);
    }
}
void GroveRenderer::add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch, int up_to_level, bool need_leaves)
{
    Visualizer v = Visualizer();
    Model *m = new Model();
    Model *lm = need_leaves ? new Model() : nullptr;
    for (int id : branch.branches)
    {
        PackedBranch &b = source->instancedCatalogue.get(id);
        if (b.level <= up_to_level)
            v.packed_branch_to_model(b, m, false);
        if (need_leaves)
            v.packed_branch_to_model(b, lm, true);
    }
    m->update();
    Instance *in = new Instance(m);
    in->addBufferCopy(branch.transforms);
    lod.instances.push_back(in);

    if (need_leaves)
    {
        lm->update();
        Instance *lin = new Instance(lm);
        lin->addBufferCopy(branch.transforms);
        lod.leaves_instances.push_back(lin);   
    }
}