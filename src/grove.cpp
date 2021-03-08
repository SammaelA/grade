#include "grove.h"
#include "tree.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "billboard_cloud.h"
#include "impostor.h"
#include <chrono>

GroveRenderer::GroveRenderer(): 
renderer({"simple_render.vs", "simple_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
rendererInstancing({"simple_render_instancing.vs", "simple_render_instancing.fs"},
                   {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"}),
lodCompute({"lod_compute.comp"},{}),
clearCompute({"clear_compute.comp"},{})
{

}
GroveRenderer::GroveRenderer(GrovePacked *_source, GroveGenerationData *_ggd, int LODs_count, std::vector<float> &max_distances,
                             bool print_perf) :
GroveRenderer()
{
    ts = Timestamp(print_perf);
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
    base_container = new Model();
    prepare_wood_types_atlas();
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
    std::vector<LodData> lods;
    std::vector<InstanceData> instances;
    std::vector<ModelData> models;
    std::vector<TypeData> types;

    //fill LODs data
    for (LOD &lod : LODs)
    {
        LodData Li;
        Li.min_max.y = lod.max_dist;

        if (!lods.empty())
        {
            lods.back().min_max.x = Li.min_max.y;
        }
        lods.push_back(Li);
    }
    //types - impostors;
    int l;
    l = 0;
    for (LOD &lod : LODs)
    {
        if (!lod.imp_rend)
            continue;
        if (lod.imp_rend->data->impostors.size() != lod.imp_rend->models.size())
        {
            logerr("broken impostor data");
            continue;
        }
        
        types_descs.push_back(TypeDescriptionForRender());
        types_descs.back().imp = lod.imp_rend;
        types_descs.back().rendDesc.cmd_size = sizeof(DrawElementsIndirectCommand);
        types_descs.back().rendDesc.current_types_offset = sizeof(currentTypesData)*types.size();
        types_descs.back().rendDesc.cmd_buffer_offset = sizeof(DrawElementsIndirectCommand)*models.size();
        types_descs.back().rendDesc.type_id = types.size();
        types_descs.back().rendDesc.max_models = lod.imp_rend->models.size();
        types_descs.back().rendDesc.base_vertex_id = base_container->indices.size();
        TypeData t;
        t.offset = models.size();
        
        int i = 0;
        for (Impostor &imp :lod.imp_rend->data->impostors)
        {
            Model *m = lod.imp_rend->models[i];
            BBox bb;
            auto cmd = model_to_base(m,bb);
            imp.id = types.size();
            IDA_to_bufer(imp.IDA, lods, instances, models, types);

            models.back().LOD = l;
            models.back().type = types.size();
            models.back().vertexes = cmd.count;
            models.back().first_index = cmd.firstIndex;
            i++;
        }
        types.push_back(t);
        l++;
    }

    l = 0;
    for (LOD &lod : LODs)
    {
        if (!lod.cloud)
            continue;

        if (lod.cloud->data->billboards.size() != lod.cloud->instances.size())
        {
            logerr("broken impostor data");
            continue;
        }

        types_descs.push_back(TypeDescriptionForRender());
        types_descs.back().bill = lod.cloud;
        types_descs.back().rendDesc.cmd_size = sizeof(DrawElementsIndirectCommand);
        types_descs.back().rendDesc.current_types_offset = sizeof(currentTypesData)*types.size();
        types_descs.back().rendDesc.cmd_buffer_offset = sizeof(DrawElementsIndirectCommand)*models.size();
        types_descs.back().rendDesc.type_id = types.size();
        types_descs.back().rendDesc.max_models = lod.cloud->instances.size();
        TypeData t;
        t.offset = models.size();

        int i = 0;
        for (BillboardData &bill :lod.cloud->data->billboards)
        {
            Model *m = lod.cloud->instances[i];
            BBox bb;
            auto cmd = model_to_base(m,bb);
            bill.id = types.size();
            IDA_to_bufer(bill.IDA, lods, instances, models, types);
            models.back().LOD = l;
            models.back().type = types.size();
            models.back().vertexes = cmd.count;
            models.back().first_index = cmd.firstIndex;
            pack_bb_to_model_data(models.back(), bb);

            i++;
        }
        types.push_back(t);
        l++;
    }

    l = 0;
    int mods = 0;
    TypeData t;
    t.offset = models.size();
    types_descs.push_back(TypeDescriptionForRender());
    types_descs.back().rendDesc.cmd_size = sizeof(DrawElementsIndirectCommand);
    types_descs.back().rendDesc.current_types_offset = sizeof(currentTypesData)*types.size();
    types_descs.back().rendDesc.cmd_buffer_offset = sizeof(DrawElementsIndirectCommand)*models.size();
    types_descs.back().rendDesc.type_id = types.size();

    for (LOD &lod : LODs)
    {
        for (Instance2 in : lod.instances)
        {
            in.id = types.size();

            IDA_to_bufer(in.ida, lods, instances, models, types);
            models.back().LOD = l;
            models.back().type = types.size();
            models.back().vertexes = in.cmd.count;
            models.back().first_index = in.cmd.firstIndex;

            pack_bb_to_model_data(models.back(),in.bbox);
            mods++;
        }
        for (Instance2 in : lod.leaves_instances)
        {
            in.id = types.size();

            IDA_to_bufer(in.ida, lods, instances, models, types);
            models.back().LOD = l;
            models.back().type = types.size();
            models.back().vertexes = in.cmd.count;
            models.back().first_index = in.cmd.firstIndex;
            
            pack_bb_to_model_data(models.back(),in.bbox);
            mods++;
        }
        l++;
    }
    types.push_back(t);
    types_descs.back().rendDesc.max_models = mods;
    for (int i=0;i<types.size();i++)
    {
        logerr("%d offset %d",i,types[i].offset);
    }
    glGenBuffers(1, &lods_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(LodData)*lods.size(), lods.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(TypeData)*types.size(), types.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &instances_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, instances_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(InstanceData)*instances.size(), instances.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &models_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, models_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(ModelData)*models.size(), models.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &cur_insts_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, cur_insts_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(currentInstancesData)*instances.size(), NULL, GL_STREAM_DRAW);

    glGenBuffers(1, &cur_models_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, cur_models_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(currentModelsData)*models.size(), NULL, GL_STREAM_DRAW);

    glGenBuffers(1, &draw_indirect_buffer);
	glBindBuffer(GL_DRAW_INDIRECT_BUFFER, draw_indirect_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, draw_indirect_buffer);
	glBufferData(GL_DRAW_INDIRECT_BUFFER, sizeof(DrawElementsIndirectCommand)*models.size(), NULL, GL_STREAM_DRAW);
	glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0); // unbind

    glGenBuffers(1, &cur_types_buf);
	glBindBuffer(GL_PARAMETER_BUFFER_ARB, cur_types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, cur_types_buf);
	glBufferData(GL_PARAMETER_BUFFER_ARB, sizeof(currentTypesData)*types.size(), NULL, GL_STREAM_DRAW);
	glBindBuffer(GL_PARAMETER_BUFFER_ARB, 0); // unbind

    total_models_count = models.size();

    lods.clear();
    instances.clear();
    models.clear();
    types.clear();

    base_container->update();
}
void GroveRenderer::pack_bb_to_model_data(ModelData &md, BBox &bb)
{
    md.x_s = glm::vec4(bb.position.x, bb.sizes.x * bb.a.x, bb.sizes.x * bb.b.x, bb.sizes.x * bb.c.x);
    md.y_s = glm::vec4(bb.position.y, bb.sizes.y * bb.a.y, bb.sizes.y * bb.b.y, bb.sizes.y * bb.c.y);
    md.z_s = glm::vec4(bb.position.z, bb.sizes.z * bb.a.z, bb.sizes.z * bb.b.z, bb.sizes.z * bb.c.z);
    md.culling = 1;
}
DrawElementsIndirectCommand GroveRenderer::model_to_base(Model *m, BBox &bb)
{
    DrawElementsIndirectCommand cmd;
    int firstVertex = base_container->positions.size()/3;
    cmd.firstIndex = base_container->indices.size();
    cmd.count = m->indices.size();
    for (GLuint ind : m->indices)
    {
        base_container->indices.push_back(firstVertex + ind);
    }
    for (GLfloat pos : m->positions)
    {
        base_container->positions.push_back(pos);
    }
    for (GLfloat norm : m->normals)
    {
        base_container->normals.push_back(norm);
    }
    for (GLfloat tex : m->colors)
    {
        base_container->colors.push_back(tex);
    }
    //logerr("%d ind", cmd.count );
    glm::vec3 mx, mn;
    mx = glm::vec3(-1e9);
    mn = glm::vec3(1e9);
    for (int i = 0; i<m->positions.size(); i+=3)
    {
        mx = glm::max(mx,glm::vec3(m->positions[i],m->positions[i+1],m->positions[i+2]));
        mn = glm::min(mn,glm::vec3(m->positions[i],m->positions[i+1],m->positions[i+2]));
    }
    BBox br_bbox;
    br_bbox.position = mn;
    br_bbox.sizes = mx - mn;
    br_bbox.a = glm::vec3(1,0,0);
    br_bbox.b = glm::vec3(0,1,0);
    br_bbox.c = glm::vec3(0,0,1);
    bb = br_bbox;
    return cmd;
}
void GroveRenderer::IDA_to_bufer(InstanceDataArrays &ida, std::vector<LodData> &lods, std::vector<InstanceData> &instances,
                                 std::vector<ModelData> &models, std::vector<TypeData> &types)
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
    models.push_back(ModelData());
    models.back().interval = glm::uvec2(st,instances.size());

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
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
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
    clearCompute.uniform("count", (uint)types_descs.size());
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, cur_types_buf);  
    glDispatchCompute(1, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
    
    ts.end("clearCompute");
    
    ts.start("lodCompute");
    lodCompute.use();
    lodCompute.uniform("lods_count", (uint)LODs.size());
    lodCompute.uniform("camera_pos", camera_pos);
    lodCompute.uniform("trans", 20.0f);
    lodCompute.uniform("projectionCamera", prc);
    lodCompute.uniform("objects_count", (uint)total_models_count);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, instances_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, models_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, cur_insts_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, cur_models_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, draw_indirect_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, cur_types_buf);

    glDispatchCompute((GLuint)ceil(total_models_count/128.0), 1, 1);
    
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
    
    ts.end("lodCompute");

    std::chrono::steady_clock::time_point tt1 = std::chrono::steady_clock::now();
    ts.start("bind");
    glBindVertexArray(base_container->vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, base_container->ibo);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, draw_indirect_buffer);
    glBindBuffer(GL_PARAMETER_BUFFER_ARB, cur_types_buf);
    ts.end("bind");
    std::chrono::steady_clock::time_point tt2 = std::chrono::steady_clock::now();
    glm::vec4 ss = glm::vec4(screen_size.x, screen_size.y, 1 / screen_size.x, 1 / screen_size.y);

    for (auto &type : types_descs)
    {
        if (type.imp)
        {
            type.imp->render(type.rendDesc, prc, camera_pos, ss);
        }
        else if (type.bill)
        {
            type.bill->render(type.rendDesc, prc, camera_pos, ss);
        }
        else
        {
            auto &mdrd = type.rendDesc;
            Texture noise = textureManager.get("noise");
            rendererInstancing.use();
            rendererInstancing.uniform("projectionCamera", prc);
            rendererInstancing.uniform("screen_size", ss);
            rendererInstancing.texture("noise", noise);
            rendererInstancing.texture("tex", atlas->tex());
            rendererInstancing.uniform("type_id", (uint)mdrd.type_id);
            rendererInstancing.uniform("camera_pos", camera_pos);
            glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                                mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
        }
    }
    ts.resolve();
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    //std::cerr << "grove rendered " << std::chrono::duration_cast<std::chrono::microseconds>(tt2 - tt1).count() << "[us]" << std::endl;
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
    //Model *m = new Model();
    //Model *lm = need_leaves ? new Model() : nullptr;
    uint ind_offset = base_container->indices.size();
    uint verts = base_container->positions.size();
    for (int id : branch.branches)
    {
        PackedBranch &b = source->instancedCatalogue.get(id);
        if (b.level <= up_to_level)
            v.packed_branch_to_model(b, base_container, false, up_to_level, glm::vec2(2*type, 0));
    }

    uint l_ind_offset = base_container->indices.size();
    uint l_verts = base_container->positions.size();
    glm::vec3 mx, mn;
    mx = glm::vec3(-1e9);
    mn = glm::vec3(1e9);
    for (int i = verts; i<l_verts; i+=3)
    {
        mx = glm::max(mx,glm::vec3(base_container->positions[i],base_container->positions[i+1],base_container->positions[i+2]));
        mn = glm::min(mn,glm::vec3(base_container->positions[i],base_container->positions[i+1],base_container->positions[i+2]));
    }
    BBox br_bbox;
    br_bbox.position = mn;
    br_bbox.sizes = mx - mn;
    br_bbox.a = glm::vec3(1,0,0);
    br_bbox.b = glm::vec3(0,1,0);
    br_bbox.c = glm::vec3(0,0,1);
    verts = l_verts;
    if (need_leaves)
    {
        for (int id : branch.branches)
        {
            PackedBranch &b = source->instancedCatalogue.get(id);
            //if (b.level <= up_to_level)
            v.packed_branch_to_model(b, base_container, true, up_to_level, glm::vec2(2*type + 1, 0));
        }
    }
    uint l_end = base_container->indices.size();
    l_verts = base_container->positions.size();
    int count = l_ind_offset - ind_offset;
    int l_count = l_end - l_ind_offset;
    
    BBox l_bbox;
    if (need_leaves)
    {
        mx = glm::vec3(-1e9);
        mn = glm::vec3(1e9);
        for (int i = verts; i<l_verts; i+=3)
        {
            //logerr("eee");
            mx = glm::max(mx,glm::vec3(base_container->positions[i],base_container->positions[i+1],base_container->positions[i+2]));
            mn = glm::min(mn,glm::vec3(base_container->positions[i],base_container->positions[i+1],base_container->positions[i+2]));
        }

        l_bbox.position = mn;
        l_bbox.sizes = mx - mn;
        l_bbox.a = glm::vec3(1,0,0);
        l_bbox.b = glm::vec3(0,1,0);
        l_bbox.c = glm::vec3(0,0,1);
        //logerr("sizes %f %f %f", l_bbox.sizes.x,l_bbox.sizes.y,l_bbox.sizes.z);
    }
    //m->update();

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

    if (count > 0 && !branch.IDA.transforms.empty())
    {
        lod.instances.push_back(Instance2(nullptr,type,branch.IDA));
        lod.instances.back().cmd.firstIndex = ind_offset;
        lod.instances.back().cmd.count = count;
        lod.instances.back().cmd.instanceCount = lod.instances.back().ida.transforms.size();
        lod.instances.back().bbox = br_bbox;
    }
    else
    {   
        //delete m;
    }
    
    if (need_leaves)
    {
        //lm->update();
        if (l_count > 0 && !branch.IDA.transforms.empty())
        {
            logerr("llllllaaaaaaa %d",l_count);
            lod.leaves_instances.push_back(Instance2(nullptr,type,branch.IDA));   
            lod.leaves_instances.back().id = lod.instances.size() - 1;
            lod.leaves_instances.back().cmd.firstIndex = l_ind_offset;
            lod.leaves_instances.back().cmd.count = l_count;
            lod.leaves_instances.back().cmd.instanceCount = lod.leaves_instances.back().ida.transforms.size();
            lod.leaves_instances.back().bbox = l_bbox;
        }
        //else
        //    delete lm;
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
void GroveRenderer::prepare_wood_types_atlas()
{
    if (!ggd)
        return;

    const int tex_size = 512;
    int num_texs = 2*ggd->types.size();
    if (num_texs <= 0)
        return;
    atlas = new TextureAtlas(tex_size,tex_size,num_texs);
    atlas->set_grid(tex_size,tex_size);
    int k = 0;

    Shader copy({"copy.vs", "copy.fs"}, {"in_Position", "in_Tex"});

    Model bm;
    std::vector<float> vertexes = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
    std::vector<float> tc = {0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0};
    std::vector<GLuint> indices = {0, 1, 2, 2, 1, 3};

    std::function<void(Model *)> _c_mip = [&](Model *h) {
        bm.positions = vertexes;
        bm.colors = tc;
        bm.indices = indices;
    };
    for (auto &type : ggd->types)
    {
        k = atlas->add_tex();
        atlas->target(k);
        bm.construct(_c_mip);
        copy.use();
        copy.texture("tex", type.wood);
        bm.render(GL_TRIANGLES);

        k = atlas->add_tex();
        atlas->target(k);
        bm.construct(_c_mip);
        copy.use();
        copy.texture("tex", type.leaf);
        bm.render(GL_TRIANGLES);
    }
    atlas->gen_mipmaps();
}