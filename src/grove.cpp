#include "grove.h"
#include "tree.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "billboard_cloud.h"
#include "impostor.h"

GroveRenderer::GroveRenderer(): 
renderer({"simple_render.vs", "simple_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
rendererInstancing({"simple_render_instancing.vs", "simple_render_instancing.fs"},
                   {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"}),
lodCompute({"lod_compute.hlsl"},{}),
clearCompute({"clear_compute.hlsl"},{})
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
        for (int j = 0; j < i - 1; j++)
        {
            auto packed_branches = _source->uniqueCatalogue.get_level(j);
            for (PackedBranch &pb : packed_branches)
            {
                Model *m = new Model();
                v.packed_branch_to_model(pb, m, false, i-2);
                LODs.back().models.push_back(std::pair<uint,Model *>(pb.type_id,m));
            }
        }
        for (InstancedBranch &b : source->instancedBranches)
        {
            add_instance_model(LODs.back(), source, b, i-2);
        }
    }

    LODs.emplace_back();
    LODs.back().cloud = nullptr;
    LODs.back().max_dist = max_distances[LODs_count - 1];

    int max_level = 3;

    for (int j = 0; j < source->uniqueCatalogue.levels(); j++)
    {
        auto packed_branches = source->uniqueCatalogue.get_level(j);
        for (PackedBranch &pb : packed_branches)
        {
            Model *m = new Model();
            v.packed_branch_to_model(pb, m, false,max_level);
            LODs.back().models.push_back(std::pair<uint,Model *>(pb.type_id,m));
        }
    }

    for (InstancedBranch &b : source->instancedBranches)
    {
        add_instance_model(LODs.back(), source, b,max_level,true);
    }

    LODs.emplace_back();
    LODs.back().cloud = nullptr;
    LODs.back().max_dist = -10;

    for (int j = 0; j < source->uniqueCatalogue.levels(); j++)
    {
        auto packed_branches = source->uniqueCatalogue.get_level(j);
        for (PackedBranch &pb : packed_branches)
        {
            Model *m = new Model();
            v.packed_branch_to_model(pb, m, false,max_level);
            LODs.back().models.push_back(std::pair<uint,Model *>(pb.type_id,m));
        }
    }

    for (InstancedBranch &b : source->instancedBranches)
    {
        add_instance_model(LODs.back(), source, b,max_level,false);
    }

    if (source->impostors.size() == 2 && LODs.size() >= 3)
    {
        LODs[0].imp_rend = new ImpostorRenderer(&(source->impostors[0]));
        LODs[1].imp_rend = new ImpostorRenderer(&(source->impostors[1]));
    }
    std::vector<LOD_info> lod_infos;
    std::vector<InstanceData> instances;
    std::vector<glm::uvec2> models_intervals;

    for (LOD &lod : LODs)
    {
        LOD_info Li;
        Li.min_max.y = lod.max_dist;
        Li.offset = instances.size();
        if (!lod_infos.empty())
        {
            lod_infos.back().min_max.x = Li.min_max.y;
        }
        lod_infos.push_back(Li);
        
        if (lod.imp_rend)
        {
            for (Impostor &imp :lod.imp_rend->data->impostors)
            {
                imp.id = models_intervals.size()/2;
                IDA_to_bufer(imp.IDA, lod_infos, instances, models_intervals);
            }
        }

        if (lod.cloud)
        {
            for (BillboardData &bill :lod.cloud->data->billboards)
            {
                bill.id = models_intervals.size()/2;
                IDA_to_bufer(bill.IDA, lod_infos, instances, models_intervals);
            }
        }


            for (int i=0;i<lod.instances.size();i++)
            {
                lod.instances[i].id = models_intervals.size()/2;
                IDA_to_bufer(lod.instances[i].ida, lod_infos, instances, models_intervals);
            }
        
        for (int i=0;i<lod.leaves_instances.size();i++)
        {
            lod.leaves_instances[i].id = lod.instances[lod.leaves_instances[i].id].id;
        }
        
    }

    //TODO - add billboards and models
    instance_models_count = models_intervals.size()/2;
    
    glGenBuffers(1, &lods_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(LOD_info)*lod_infos.size(), lod_infos.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &inst_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, inst_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(InstanceData)*instances.size(), instances.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &indexes_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, indexes_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, 16*sizeof(uint)*instances.size(), NULL, GL_DYNAMIC_COPY);

    glGenBuffers(1, &intervals_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, intervals_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, 2*sizeof(uint)*models_intervals.size(), models_intervals.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &counts_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, counts_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, 2*sizeof(uint)*models_intervals.size(), NULL, GL_DYNAMIC_READ);

    lod_infos.clear();
    instances.clear();
    models_intervals.clear();
}
void GroveRenderer::IDA_to_bufer(InstanceDataArrays &ida, std::vector<LOD_info> &lod_infos, 
                                 std::vector<InstanceData> &instances, std::vector<glm::uvec2> &models_intervals)
{
    uint st = instances.size();
    if (ida.centers_par.size() != ida.centers_self.size() || ida.centers_par.size() != ida.transforms.size())
    {
        logerr("invalid Instance Data Arrays");
        return;
    }
    for (int i=0; i<ida.centers_self.size(); i++)
    {
        instances.push_back(InstanceData());
        instances.back().projection_camera = ida.transforms[i];
        instances.back().center_par = glm::vec4(ida.centers_par[i],1);
        instances.back().center_self = glm::vec4(ida.centers_self[i],1);
    }
    models_intervals.push_back(glm::uvec2(st,instances.size()));
    models_intervals.push_back(glm::uvec2(100,100));
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
            //if (LODs[i].instances[j].m)
            //   delete LODs[i].instances[j].m;
        }
    }
    LODs.clear();
    source = nullptr;
}
void GroveRenderer::render(int explicit_lod, glm::mat4 prc, glm::vec3 camera_pos, glm::vec2 screen_size)
{
    if (LODs.size() == 0)
        return;
    std::vector<int> lods_to_render;
    std::vector<uint> counts;

    if (explicit_lod == -1)
    {
        float len = glm::length(source->center - camera_pos);
        float sz = sqrt(ggd->size.x * ggd->size.x + ggd->size.z * ggd->size.z);
        for (int i = 0; i < LODs.size() - 1; i++)
        {
            if (LODs[i + 1].max_dist < len + sz && LODs[i].max_dist >= len - sz)
                lods_to_render.push_back(i);
        }
    }
    else
    {
        if (explicit_lod >= LODs.size())
        {
            logerr("trying to render grove with wrong explicit LOD number %d. Grove has %d LODs", explicit_lod, LODs.size());
            explicit_lod = LODs.size() - 1;
        }
        else if (explicit_lod < 0)
        {
            logerr("trying to render grove with wrong explicit LOD number %d. Grove has %d LODs", explicit_lod, LODs.size());
            explicit_lod = 0;
        }
        lods_to_render.push_back(explicit_lod);
    }

    ts.start("clearCompute");
    clearCompute.use();
    clearCompute.uniform("count", (uint)instance_models_count);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, counts_buf);  
    glDispatchCompute(1, 1, 1);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    ts.end("clearCompute");
    ts.start("lodCompute");
    lodCompute.use();
    lodCompute.uniform("lods_count", (uint)LODs.size());
    lodCompute.uniform("camera_pos", camera_pos);
    lodCompute.uniform("trans", 20.0f);
    lodCompute.uniform("objects_count", (uint)instance_models_count);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, inst_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, indexes_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, intervals_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, counts_buf);

    glDispatchCompute((GLuint)ceil(instance_models_count/128.0), 1, 1);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    ts.end("lodCompute");

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, 0);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, counts_buf);
    GLuint *ptr = nullptr;
    ptr = (GLuint *)glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);

    for (int i = 0; i < 4 * instance_models_count; i += 4)
    {
        counts.push_back(ptr[i]);
        //logerr("%d) %d",i/4,counts.back());
    }
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    for (int lod : lods_to_render)
    {
        ts.start("lod"+std::to_string(lod)+"_models");
        renderer.use();
        renderer.uniform("projectionCamera", prc);
        for (int j = 0; j < LODs[lod].models.size(); j++)
        {
            Model *m = LODs[lod].models[j].second;
            renderer.texture("tex", ggd->types[LODs[lod].models[j].first].wood);
            renderer.uniform("model", m->model);
            m->update();
            m->render(GL_TRIANGLES);
        }
        ts.end("lod"+std::to_string(lod)+"_models");
        float mx = LODs[lod].max_dist == -10 ? 1000 : LODs[lod].max_dist;
        float mn = lod + 1 == LODs.size() ? 0 : LODs[lod + 1].max_dist;
        glm::vec2 mn_mx = glm::vec2(mn, mx);
        Texture noise = textureManager.get("noise");
        glm::vec4 ss = glm::vec4(screen_size.x, screen_size.y, 1 / screen_size.x, 1 / screen_size.y);

        ts.start("lod"+std::to_string(lod)+"_billboards");
        if (LODs[lod].cloud)
            LODs[lod].cloud->render(counts, prc, camera_pos, mn_mx, ss);
        ts.end("lod"+std::to_string(lod)+"_billboards");

        ts.start("lod"+std::to_string(lod)+"_imposters");
        if (LODs[lod].imp_rend)
            LODs[lod].imp_rend->render(counts, prc, camera_pos, mn_mx, ss);
        ts.end("lod"+std::to_string(lod)+"_imposters");

        ts.start("lod"+std::to_string(lod)+"_instances");
        rendererInstancing.use();
        rendererInstancing.uniform("projectionCamera", prc);
        rendererInstancing.uniform("screen_size", ss);
        rendererInstancing.texture("noise", noise);
        rendererInstancing.uniform("LOD_dist_min_max", mn_mx);
        rendererInstancing.uniform("camera_pos", camera_pos);
        for (auto &in : LODs[lod].instances)
        {
            rendererInstancing.texture("tex", ggd->types[in.type].wood);
            rendererInstancing.uniform("id",(uint)in.id);
            if (in.id >= 0 && in.id < counts.size() && counts[in.id] > 0)
            {
                Model *m = in.m;
                if(m && m->indexed)
                {
                    glBindVertexArray(m->vao);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m->ibo);
                    glDrawElementsInstanced(GL_TRIANGLES, m->SIZE, GL_UNSIGNED_INT, 0, counts[in.id]);
                }
            }
        }

        for (auto &in : LODs[lod].leaves_instances)
        {
            rendererInstancing.texture("tex", ggd->types[in.type].leaf);
            rendererInstancing.uniform("id",(uint)in.id);
            if (in.id >= 0 && in.id < counts.size() && counts[in.id] > 0)
            {
                Model *m = in.m;
                if(m && m->indexed)
                {
                    glBindVertexArray(m->vao);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m->ibo);
                    glDrawElementsInstanced(GL_TRIANGLES, m->SIZE, GL_UNSIGNED_INT, 0, counts[in.id]);
                }
            }
        }
        ts.end("lod"+std::to_string(lod)+"_instances");
    }
    ts.resolve();
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
            v.packed_branch_to_model(b, m, false, up_to_level);
        if (need_leaves)
            v.packed_branch_to_model(b, lm, true, up_to_level);
    }
    m->update();

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
        
        branch.IDA.centers_par[i] = center;
    }

    if (!m->positions.empty() && !branch.IDA.transforms.empty())
    {
        lod.instances.push_back(Instance2(m,type,branch.IDA));
    }
    else
    {   
        delete m;
    }
    
    if (need_leaves)
    {
        lm->update();
        if (!lm->positions.empty() && !branch.IDA.transforms.empty())
        {
            lod.leaves_instances.push_back(Instance2(lm,type,branch.IDA));   
            lod.leaves_instances.back().id = lod.instances.size() - 1;
        }
        else
            delete lm;
    }
}
GroveRenderer::Instance2::Instance2(Model *_m, uint _type, InstanceDataArrays &_ida):
ida(_ida)
{
    m = _m;
    type = _type;
}
GroveRenderer::Instance2::~Instance2()
{
    //if (m)
    //    delete m;
}