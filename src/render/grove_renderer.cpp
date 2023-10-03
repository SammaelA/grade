#include "grove_renderer.h"
#include "core/tree.h"
#include "tree_utils/tree_modeling.h"
#include "tinyEngine/engine.h"
#include "tree_utils/billboard_cloud.h"
#include "tree_utils/impostor.h"
#include "tinyEngine/camera.h"
#include "render/billboard_cloud_renderer.h"
#include "render/impostor_renderer.h"
#include "generation/generation_task.h"
#include <chrono>

GroveRenderer::GroveRenderer(const GrovePacked *_source, AABB2D _scene_bbox, const std::vector<TreeTypeData> &_types, 
                             int LODs_count, std::vector<float> &max_distances, bool print_perf, Precision precision) :
renderer({"simple_render.vs", "simple_render.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
rendererInstancing({"simple_render_instancing.vs", "simple_render_instancing.fs"},
                   {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"}),
shadowRendererInstancing({"simple_render_instancing.vs", "depth_billboard_array.fs"},
                   {"in_Position", "in_Normal", "in_Tex", "in_Center_par", "in_Center_self", "in_Model"}),
lodCompute({"lod_compute.comp"},{}),
cellsCompute({"cells_compute.comp"},{}),
types(_types)
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
    
    source = (GrovePacked *)_source;
    scene_bbox = _scene_bbox;
    debugl(10,"creating grove renderer with %d LODs\n", _source->clouds.size());
    base_container = new Model();
    prepare_wood_types_atlas();
    for (int i = 0; i < _source->clouds.size(); i++)
    {
        LODs.emplace_back();
        LODs.back().max_dist = max_distances[i];
        LODs.back().cloud = new BillboardCloudRenderer(&_source->clouds[i]);
        LODs.back().cloud->set_render_mode(BillboardCloudRenderer::ONLY_INSTANCES);

        for (InstancedBranch &b : source->instancedBranches)
        {
            add_instance_model(LODs.back(), source, b, i-2);
        }
    }
    if (precision >= Precision::MEDIUM)
    {
        LODs.emplace_back();
        LODs.back().cloud = nullptr;
        LODs.back().max_dist = max_distances.empty() ? 0 : max_distances.back();

        int max_level = 12;

        for (InstancedBranch &b : source->instancedBranches)
        {
            add_instance_model(LODs.back(), source, b,max_level,true);
        }
    }
    if (precision >= Precision::DEBUG)
    {
        LODs.emplace_back();
        LODs.back().cloud = nullptr;
        LODs.back().max_dist = -10;
        int max_level = 3;

        for (InstancedBranch &b : source->instancedBranches)
        {
            add_instance_model(LODs.back(), source, b,max_level,false);
        }
    }
    if (source->impostors.size() == 2 && LODs.size() >= 3)
    {
        LODs[0].imp_rend = new ImpostorRenderer(&(source->impostors[1]));
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

            IDA_to_bufer(in.ida, lods, instances, models, types, false);
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

            IDA_to_bufer(in.ida, lods, instances, models, types, true);
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

    //save cells_info;
    cellsInfo.x_cells = 32;
    cellsInfo.y_cells = 32;
    cellsInfo.x_size = (scene_bbox.max_pos.x - scene_bbox.min_pos.x)/cellsInfo.x_cells;
    cellsInfo.y_size = (scene_bbox.max_pos.y - scene_bbox.min_pos.y)/cellsInfo.y_cells;
    cellsInfo.start_pos = glm::vec3(scene_bbox.min_pos.x, 0, scene_bbox.min_pos.y);
    for (int i=0;i<instances.size();i++)
    {
        int par_cell_id, self_cell_id;
        {
        glm::vec3 ps = instances[i].center_self;
        int xc = (ps.x - cellsInfo.start_pos.x)/cellsInfo.x_size;
        int yc = (ps.y - cellsInfo.start_pos.y)/cellsInfo.y_size;
        self_cell_id = xc*cellsInfo.y_cells + yc;
        }
        {
        glm::vec3 ps = instances[i].center_par;
        int xc = (ps.x - cellsInfo.start_pos.x)/cellsInfo.x_size;
        int yc = (ps.y - cellsInfo.start_pos.y)/cellsInfo.y_size;
        par_cell_id = xc*cellsInfo.y_cells + yc;
        }
        if (par_cell_id == self_cell_id && par_cell_id >= 0 && par_cell_id < cellsInfo.x_cells*cellsInfo.y_cells)
        {
            instances[i].center_self.w = par_cell_id;
        }
        else
        {
            instances[i].center_self.w = -1;
        }
    }
    lods_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(LodData)*lods.size(), lods.data(), GL_STATIC_DRAW);

    types_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(TypeData)*types.size(), types.data(), GL_STATIC_DRAW);

    instances_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, instances_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(InstanceData)*instances.size(), instances.data(), GL_STATIC_DRAW);

    models_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, models_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(ModelData)*models.size(), models.data(), GL_STATIC_DRAW);

    cur_insts_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, cur_insts_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(currentInstancesData)*instances.size(), NULL, GL_STREAM_DRAW);

    cur_models_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, cur_models_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(currentModelsData)*models.size(), NULL, GL_STREAM_DRAW);

    draw_indirect_buffer = create_buffer();
	glBindBuffer(GL_DRAW_INDIRECT_BUFFER, draw_indirect_buffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, draw_indirect_buffer);
	glBufferData(GL_DRAW_INDIRECT_BUFFER, sizeof(DrawElementsIndirectCommand)*models.size(), NULL, GL_STREAM_DRAW);
	glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0); // unbind

    cur_types_buf = create_buffer();
	glBindBuffer(GL_PARAMETER_BUFFER_ARB, cur_types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, cur_types_buf);
	glBufferData(GL_PARAMETER_BUFFER_ARB, sizeof(currentTypesData)*types.size(), NULL, GL_STREAM_DRAW);
	glBindBuffer(GL_PARAMETER_BUFFER_ARB, 0); // unbind

    cells_buf = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, cells_buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(CellInfo)*(cellsInfo.x_cells*cellsInfo.y_cells), NULL, GL_STREAM_DRAW);

    total_models_count = models.size();

    lods.clear();
    instances.clear();
    models.clear();
    types.clear();

    base_container->update();
    long size = (base_container->SIZE)*8*sizeof(float) + sizeof(LodData)*lods.size() + sizeof(TypeData)*types.size() +
                sizeof(InstanceData)*instances.size() + sizeof(ModelData)*models.size() +
                sizeof(currentInstancesData)*instances.size() + sizeof(currentModelsData)*models.size() +
                sizeof(DrawElementsIndirectCommand)*models.size() + sizeof(currentTypesData)*types.size();
    debugl(1,"total size %f Mbytes\n",1e-6*size);
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
                                 std::vector<ModelData> &models, std::vector<TypeData> &, bool is_leaf)
{
    uint st = instances.size();
    if (ida.centers_par.size() != ida.centers_self.size() || ida.centers_par.size() != ida.transforms.size())
    {
        logerr("invalid Instance Data Arrays");
        return;
    }
    for (int i=0; i<ida.centers_self.size(); i++)
    {
        float type_slice = 0;
        if (ida.type_ids.size() > i)
            type_slice = is_leaf ? types[ida.type_ids[i]].leaf_id : types[ida.type_ids[i]].wood_id;
        instances.push_back(InstanceData());
        instances.back().projection_camera = ida.transforms[i];
        instances.back().center_par = glm::vec4(ida.centers_par[i], type_slice);
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
        if (LODs[i].imp_rend)
            delete LODs[i].imp_rend;
    }

    #define DELBUF(a) if (a) { delete_buffer((a)); a = 0;}

    DELBUF(lods_buf);
    DELBUF(instances_buf);
    DELBUF(models_buf);
    DELBUF(types_buf);
    DELBUF(cur_insts_buf);
    DELBUF(cur_models_buf);
    DELBUF(cur_types_buf);
    DELBUF(draw_indirect_buffer);
    DELBUF(cells_buf);

    LODs.clear();
    delete base_container;
    delete atlas;
    source = nullptr;
}

void GroveRenderer::render(int explicit_lod, glm::mat4 &projection, glm::mat4 &view, Camera &camera, glm::vec2 screen_size, 
                           DirectedLight &light, GroveRendererDebugParams dbgpar, glm::mat4 &shadow_tr,
                           GLuint shadow_tex, bool to_shadow)
{
    glm::vec3 camera_pos = camera.pos;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    if (LODs.size() == 0)
        return;
    std::vector<int> lods_to_render;
    std::vector<uint> counts;

    if (explicit_lod == -1)
    {
        glm::vec2 center = 0.5f*(scene_bbox.max_pos + scene_bbox.min_pos);
        float len = glm::length(center - glm::vec2(camera_pos.x, camera_pos.z));
        glm::vec2 sz_2 = 0.5f*(scene_bbox.max_pos - scene_bbox.min_pos);
        float sz = sqrt(sz_2.x*sz_2.x + sz_2.y*sz_2.y);
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
            static int last_lod_err = -1;
            if (explicit_lod != last_lod_err)
            {
                last_lod_err = explicit_lod;
                logerr("trying to render grove with wrong explicit LOD number %d. Grove has %d LODs", explicit_lod, LODs.size());
            }
            explicit_lod = LODs.size() - 1;
        }
        else if (explicit_lod < 0)
        {
            logerr("trying to render grove with wrong explicit LOD number %d. Grove has %d LODs", explicit_lod, LODs.size());
            explicit_lod = 0;
        }
        lods_to_render.push_back(explicit_lod);
    }
    bool camera_changed = (length(camera.pos - prev_camera.pos) > 1e-4) || 
                          (length(camera.front - prev_camera.front) > 1e-4) || 
                          (length(camera.up - prev_camera.up) > 1e-4);
    if (true || frames<10 || (frames % 2 == 0 && camera_changed))
    {
        //ts.start("cellsCompute");
        cellsCompute.use();
        cellsCompute.uniform("count", (uint)types_descs.size());
        cellsCompute.uniform("grid_sizes", glm::vec2(cellsInfo.x_size,cellsInfo.y_size));
        cellsCompute.uniform("grid_counts_x", (uint)cellsInfo.x_cells);
        cellsCompute.uniform("grid_counts_y", (uint)cellsInfo.y_cells);
        cellsCompute.uniform("pos", cellsInfo.start_pos);
        cellsCompute.uniform("camera_pos", camera_pos);
        cellsCompute.uniform("trans", 20.0f);
        cellsCompute.uniform("lods_cnt", 5);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, cur_types_buf);  
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, cells_buf);


        glDispatchCompute(1, 1, 1);
        //ts.end("cellsCompute");
        
        //ts.start("flushCompute");
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
        
        //ts.end("flushCompute");
        
        //ts.start("lodCompute");
        lodCompute.use();
        lodCompute.uniform("lods_count", (uint)LODs.size());
        lodCompute.uniform("camera_pos", camera_pos);
        lodCompute.uniform("trans", 20.0f);
        lodCompute.uniform("projectionCamera", projection * view);
        lodCompute.uniform("objects_count", (uint)total_models_count);
        lodCompute.uniform("forced_lod",explicit_lod);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, lods_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, instances_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, models_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, cur_insts_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, cur_models_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, draw_indirect_buffer);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, cur_types_buf);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, cells_buf);

        glDispatchCompute((GLuint)ceil(total_models_count/128.0), 1, 1);
        
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
        
        //ts.end("lodCompute");
        prev_camera = camera;
    }
    std::chrono::steady_clock::time_point tt1 = std::chrono::steady_clock::now();
    //ts.start("bind");
    glBindVertexArray(base_container->vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, base_container->ibo);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, draw_indirect_buffer);
    glBindBuffer(GL_PARAMETER_BUFFER_ARB, cur_types_buf);
    ////ts.end("bind");
    std::chrono::steady_clock::time_point tt2 = std::chrono::steady_clock::now();
    glm::vec4 ss = glm::vec4(screen_size.x, screen_size.y, 1 / screen_size.x, 1 / screen_size.y);


    auto rend_lod = [&](TypeDescriptionForRender &type)
    {
        if (type.imp)
        {
            type.imp->render(type.rendDesc, projection, view, light, shadow_tr, shadow_tex, camera_pos, ss, to_shadow);
        }
        else if (type.bill)
        {
            type.bill->render(type.rendDesc, projection, view, light, shadow_tr, shadow_tex, camera_pos, ss, to_shadow);
        }
        else
        {
            auto &mdrd = type.rendDesc;
            Texture noise = engine::textureManager->get("noise");
            Shader &shader = to_shadow ? shadowRendererInstancing : rendererInstancing;
            shader.use();

            shader.uniform("projection", projection);
            shader.uniform("view", view);
            shader.uniform("screen_size", ss);
            shader.texture("noise", noise);
            shader.texture("tex", atlas->tex(0));
            shader.uniform("type_id", (uint)mdrd.type_id);
            shader.uniform("camera_pos", camera_pos);
            shader.uniform("debug_model_id",dbgpar.need_focus_model ? dbgpar.model_focused : -1);
            if (to_shadow)
                shader.uniform("opaqueness",0.3f);
            glMultiDrawElementsIndirectCountARB(GL_TRIANGLES, GL_UNSIGNED_INT, (void *)mdrd.cmd_buffer_offset,
                                                mdrd.current_types_offset, mdrd.max_models, mdrd.cmd_size);
        }
    };
    if (explicit_lod >= 0)
      rend_lod(types_descs.back());
    else
    {
      for (auto &type : types_descs)
      {
        rend_lod(type);
      }
    }
    //ts.resolve();
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    frames++;
}

void GroveRenderer::add_instance_model(LOD &lod, GrovePacked *source, InstancedBranch &branch, int up_to_level, bool need_leaves)
{
    if (branch.branches.empty())
        return;
    //clusterization process guarantees that type of all branches in instance
    //will be the same
    uint type = source->instancedCatalogue.get(branch.branches.front()).type_id;
    glm::vec3 pos = source->instancedCatalogue.get(branch.branches.front()).joints.front().pos;
    uint ind_offset = base_container->indices.size();
    uint verts = base_container->positions.size();
    for (int id : branch.branches)
    {
        int type_slice = types[type].wood_id;
        if (id < 0)
        {
            logerr("invalid id = %d", id);
            continue;//invalid id - TODO fix it
        }
        PackedBranch &b = source->instancedCatalogue.get(id);
        if (b.level <= up_to_level && !b.joints.empty())
            visualizer::packed_branch_to_model(b, base_container, false, 2, glm::vec2(type_slice, 0));
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
        int type_slice = types[type].leaf_id;
        for (int id : branch.branches)
        {
            if (id < 0)
            {
                logerr("invalid id = %d", id);
                continue;//invalid id - TODO fix it
            }
            PackedBranch &b = source->instancedCatalogue.get(id);
            if (!b.joints.empty())
                visualizer::packed_branch_to_model(b, base_container, true, 2, glm::vec2(type_slice, 0));
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
            mx = glm::max(mx,glm::vec3(base_container->positions[i],base_container->positions[i+1],base_container->positions[i+2]));
            mn = glm::min(mn,glm::vec3(base_container->positions[i],base_container->positions[i+1],base_container->positions[i+2]));
        }

        l_bbox.position = mn;
        l_bbox.sizes = mx - mn;
        l_bbox.a = glm::vec3(1,0,0);
        l_bbox.b = glm::vec3(0,1,0);
        l_bbox.c = glm::vec3(0,0,1);
    }
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

    if (need_leaves)
    {
        if (l_count > 0 && !branch.IDA.transforms.empty())
        {
            lod.leaves_instances.push_back(Instance2(nullptr,type,branch.IDA));   
            lod.leaves_instances.back().id = lod.instances.size() - 1;
            lod.leaves_instances.back().cmd.firstIndex = l_ind_offset;
            lod.leaves_instances.back().cmd.count = l_count;
            lod.leaves_instances.back().cmd.instanceCount = lod.leaves_instances.back().ida.transforms.size();
            lod.leaves_instances.back().bbox = l_bbox;
        }
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
    const int tex_size = 512;
    int num_texs = 0;
    std::map<GLuint, std::pair<Texture,int>> tex_map;
    for (auto &type : types)
    {
        auto it = tex_map.find(type.wood.texture);
        if (it == tex_map.end())
        {
            tex_map.emplace(type.wood.texture,std::pair<Texture,int>(type.wood,num_texs));
            type.wood_id = num_texs;
            num_texs++;
        }
        else
        {
            type.wood_id = it->second.second;
        }
        it = tex_map.find(type.leaf.texture);
        if (it == tex_map.end())
        {
            tex_map.emplace(type.leaf.texture,std::pair<Texture,int>(type.leaf,num_texs));
            type.leaf_id = num_texs;
            num_texs++;
        }
        else
        {
            type.leaf_id = it->second.second;
        }
    }
    if (num_texs <= 0)
        return;
    atlas = new TextureAtlas(tex_size,tex_size,num_texs);
    atlas->set_grid(tex_size,tex_size);
    Texture &t = atlas->tex(0);
    glBindTexture(t.type, t.texture);
    glTexParameteri(t.type, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(t.type, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(t.type, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

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
    for (auto &pr : tex_map)
    {
        k = atlas->add_tex();
        atlas->target(k,0);
        bm.construct(_c_mip);
        copy.use();
        copy.texture("tex", pr.second.first);
        bm.render(GL_TRIANGLES);

        for (auto &type : types)
        {
            if (type.leaf_id == pr.second.second)
                type.leaf_id = k+1000;
            
            if (type.wood_id == pr.second.second)
                type.wood_id = k+1000;
        }
    }
    for (auto &type : types)
        {
            type.leaf_id = -1000+type.leaf_id;
            type.wood_id = -1000+type.wood_id;
    }

    atlas->gen_mipmaps();
    engine::textureManager->save_png(atlas->tex(0), "wl_atlas");
}