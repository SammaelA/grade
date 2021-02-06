#include "grove.h"
#include "tree.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "billboard_cloud.h"

GroveRenderer::GroveRenderer(): 
renderer({"simple_render.vs", "simple_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
rendererInstancing({"simple_render_instancing.vs", "simple_render_instancing.fs"},
                   {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"})
{

}
GroveRenderer::GroveRenderer(GrovePacked *_source, GroveGenerationData *_ggd, int LODs_count, std::vector<float> &max_distances) :
GroveRenderer()
{
    if (LODs_count != _source->clouds.size() + 1 || max_distances.size() != _source->clouds.size() + 1)
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
    ggd = _ggd;
    debug("creating grove renderer with %d LODs\n", _source->clouds.size());
    for (int i = 0; i < _source->clouds.size(); i++)
    {
        LODs.emplace_back();
        LODs.back().max_dist = max_distances[i];
        LODs.back().cloud = new BillboardCloudRenderer(&_source->clouds[i]);
        LODs.back().cloud->set_render_mode(BillboardCloudRenderer::ONLY_INSTANCES);
        for (int j = 0; j < i; j++)
        {
            auto packed_branches = _source->uniqueCatalogue.get_level(j);
            for (PackedBranch &pb : packed_branches)
            {
                Model *m = new Model();
                v.packed_branch_to_model(pb, m, false);
                LODs.back().models.push_back(std::pair<uint,Model *>(pb.type_id,m));
            }
        }
        for (InstancedBranch &b : source->instancedBranches)
        {
            add_instance_model(LODs.back(), source, b, i-1);
        }
    }

    LODs.emplace_back();
    LODs.back().cloud = nullptr;
    LODs.back().max_dist = max_distances[LODs_count - 1];

    for (int j = 0; j < source->uniqueCatalogue.levels(); j++)
    {
        auto packed_branches = source->uniqueCatalogue.get_level(j);
        for (PackedBranch &pb : packed_branches)
        {
            Model *m = new Model();
            v.packed_branch_to_model(pb, m, false);
            LODs.back().models.push_back(std::pair<uint,Model *>(pb.type_id,m));
        }
    }

    for (InstancedBranch &b : source->instancedBranches)
    {
        add_instance_model(LODs.back(), source, b,1000,true);
    }

    LODs.emplace_back();
    LODs.back().cloud = nullptr;
    LODs.back().max_dist = 0*max_distances[LODs_count - 1];

    for (int j = 0; j < source->uniqueCatalogue.levels(); j++)
    {
        auto packed_branches = source->uniqueCatalogue.get_level(j);
        for (PackedBranch &pb : packed_branches)
        {
            Model *m = new Model();
            v.packed_branch_to_model(pb, m, false);
            LODs.back().models.push_back(std::pair<uint,Model *>(pb.type_id,m));
        }
    }

    for (InstancedBranch &b : source->instancedBranches)
    {
        add_instance_model(LODs.back(), source, b,1000,false);
    }
}
GroveRenderer::~GroveRenderer()
{
    for (int i=0;i<LODs.size();i++)
    {
        for (int j=0;j < LODs[j].models.size();j++)
            delete LODs[i].models[j].second;
        if (LODs[i].cloud)
            delete LODs[i].cloud;
        for (int j=0;j<LODs[i].instances.size();j++)
        {
            delete LODs[i].instances[j].second->m;
            delete LODs[i].instances[j].second;
        }
    }
    LODs.clear();
    source = nullptr;
}
void GroveRenderer::render_auto_LOD(glm::mat4 prc, glm::vec3 camera_pos, glm::vec2 screen_size)
{
    float len = glm::length(source->center - camera_pos);
    float sz = sqrt(ggd->size.x*ggd->size.x + ggd->size.z*ggd->size.z);
    for (int i=0;i<LODs.size() - 1;i++)
    {
        if (LODs[i + 1].max_dist < len + sz && LODs[i].max_dist >= len - sz)
            render(i,prc,camera_pos, screen_size);
    }
    /*
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
        render(lod,prc,camera_pos, screen_size);
    */
}
void GroveRenderer::render(int lod, glm::mat4 prc, glm::vec3 camera_pos, glm::vec2 screen_size)
{
    if (LODs.size()==0)
        return;
    if (lod == -1)
    {
        render_auto_LOD(prc,camera_pos, screen_size);
        return;
    }
    if (lod < 0 || lod >= LODs.size())
    {
        logerr("trying to render grove with wrong LOD number %d. Grove has %d LODs",lod, LODs.size());
        //return;
        lod = LODs.size() - 1;
    }
    renderer.use();
    renderer.uniform("projectionCamera", prc);
    for (int j=0;j < LODs[lod].models.size();j++)
    {
        Model *m = LODs[lod].models[j].second;
        renderer.texture("tex", ggd->types[LODs[lod].models[j].first].wood);
        renderer.uniform("model", m->model);
        m->update();
        m->render(GL_TRIANGLES);
    }
    float mx = LODs[lod].max_dist;
    float mn = lod + 1 == LODs.size() ? 0 : LODs[lod + 1].max_dist;
    glm::vec2 mn_mx = glm::vec2(mn,mx);
    Texture noise = textureManager.get("noise");
    glm::vec4 ss = glm::vec4(screen_size.x,screen_size.y,1/screen_size.x,1/screen_size.y);
    if (LODs[lod].cloud)
        LODs[lod].cloud->render(prc,camera_pos,mn_mx,ss);
    
    rendererInstancing.use();
    rendererInstancing.uniform("projectionCamera", prc);
    rendererInstancing.uniform("screen_size",ss);
    rendererInstancing.texture("noise",noise);
    rendererInstancing.uniform("LOD_dist_min_max",mn_mx);
    rendererInstancing.uniform("camera_pos",camera_pos);
    for (auto &in : LODs[lod].instances)
    {
        rendererInstancing.texture("tex", ggd->types[in.first].wood);
        Model *m = (Model *)(in.second->m);
        m->update();
        in.second->render(GL_TRIANGLES);
    }

    for (auto &in : LODs[lod].leaves_instances)
    {
        rendererInstancing.texture("tex", ggd->types[in.first].leaf);
        Model *m = (Model *)(in.second->m);
        m->update();
        in.second->render(GL_TRIANGLES);
    }
}
void GroveRenderer::add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch, int up_to_level, bool need_leaves)
{
    if (branch.branches.empty())
        return;
    //clusterization process guarantees that type of all branches in instance
    //will be the same
    uint type = source->instancedCatalogue.get(branch.branches.front()).type_id;
    glm::vec3 pos = source->instancedCatalogue.get(branch.branches.front()).joints.front().pos;
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
    std::vector<glm::vec4> pt1,pt2;
    if (branch.IDA.transforms.size() != branch.IDA.centers_par.size())
    {
        logerr("sizes do not match %d %d",branch.IDA.transforms.size(),branch.IDA.centers_par.size());
    }
    for (int i=0;i<branch.IDA.transforms.size();i++)
    {
        glm::vec4 center = glm::vec4(branch.IDA.centers_par[i],1);
        if (branch.branches.size() == 1)
        {
            int level = source->instancedCatalogue.get(branch.branches.front()).level;
            if (level <= up_to_level)
                center = glm::vec4(branch.IDA.centers_self[i],1);
        }
        if (up_to_level >= base_level)
            center = glm::vec4(branch.IDA.centers_self[i],1);
        pt1.push_back(center);
        pt2.push_back(glm::vec4(branch.IDA.centers_self[i],1));
    }
    in->addBufferCopy(pt1);
    in->addBufferCopy(pt2);
    in->addBufferCopy(branch.IDA.transforms);
    if (!m->positions.empty() && !branch.IDA.transforms.empty())
        lod.instances.push_back(std::pair<uint, Instance *>(type,in));
    else
        delete in;
    
    if (need_leaves)
    {
        lm->update();
        Instance *lin = new Instance(lm);
        lin->addBufferCopy(pt1);
        lin->addBufferCopy(pt2);
        lin->addBufferCopy(branch.IDA.transforms);
        if (!lm->positions.empty() && !branch.IDA.transforms.empty())
            lod.leaves_instances.push_back(std::pair<uint, Instance *>(type,lin));   
        else
            delete lin;
    }
}