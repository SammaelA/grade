#include "general_gltf_writer.h"
#include "common_utils/utility.h"
#include "generation/generation_settings.h"
namespace gltf
{
    void GeneralGltfWriter::clear()
    {
        models.clear();
        cameras.clear();

        for (int i=0;i<temp_models.size();i++)
            delete temp_models[i];
        temp_models.clear();
    }
    void GeneralGltfWriter::add_model(Model *m)
    {
        models.push_back(ModelData{m,nullptr,{m->model}});
    }
    void GeneralGltfWriter::add_packed_grove(GrovePacked &grove, GroveGenerationData &ggd)
    {
        Visualizer v;
        int start = temp_models.size();

        for (auto &inst : grove.instancedBranches)
        {
            //TODO: support ofr different type ids
            std::vector<int> type_ids = {inst.IDA.type_ids.front()};
            InstanceDataArrays &ida = inst.IDA;
            for (int type_id : type_ids)
            {
                //branches
                Model *m = new Model();
                
                for (int br_id : inst.branches)
                {
                    auto &br = grove.instancedCatalogue.get(br_id);
                    if (br.type_id != type_id)
                        continue;
                    v.packed_branch_to_model(br,m,false,3);
                }
                if (!m->positions.empty())
                {
                    temp_models.push_back(m);
                    ModelData md;
                    md.m = m;
                    md.transforms = ida.transforms;
                    md.t = &(ggd.types[type_id].wood);
                    models.push_back(md);
                }

                //leaves
                Model *m_l = new Model();
                
                for (int br_id : inst.branches)
                {
                    auto &br = grove.instancedCatalogue.get(br_id);
                    if (br.type_id != type_id)
                        continue;
                    v.packed_branch_to_model(br,m_l,true,3);
                }
                if (!m_l->positions.empty())
                {
                    temp_models.push_back(m_l);
                    ModelData md;
                    md.m = m_l;
                    md.transforms = ida.transforms;
                    md.t = &(ggd.types[type_id].leaf);
                    models.push_back(md);
                }
            }
        }
    }
    void GeneralGltfWriter::convert_to_gltf(std::string name)
    {
        FullData fullData;
        asset_name = name;

        fullData.gltf_file.main_scene = 0;
        fullData.gltf_file.scenes.emplace_back();

        auto &scene = fullData.gltf_file.scenes[fullData.gltf_file.main_scene];

        //create textures, samplers and materials
        std::vector<::Texture *> unique_textures;
        for (ModelData &md : models)
        {
            if (md.m && md.t && !md.transforms.empty() && !md.t->origin.empty())
            {
                bool exist = false;
                for (int i=0;i<unique_textures.size();i++)
                {
                    if (unique_textures[i]->origin == md.t->origin)
                    {
                        exist = true;
                        md.material_id = i;
                        break;
                    }
                }
                if (!exist)
                {
                    unique_textures.push_back(md.t);
                    md.material_id = unique_textures.size() - 1;
                }
            }
        }
        if (!unique_textures.empty())
        {
            fullData.gltf_file.samplers.push_back(Sampler());
            fullData.textures_files.resize(unique_textures.size());
            fullData.gltf_file.images.resize(unique_textures.size());
            fullData.gltf_file.textures.resize(unique_textures.size());
            fullData.gltf_file.materials.resize(unique_textures.size());

            for (int i = 0;i<unique_textures.size();i++)
            {
                auto *t = unique_textures[i];

                fullData.textures_files[i].existed_texture = true;
                fullData.textures_files[i].file_name = t->origin;

                fullData.gltf_file.images[i].picture = &fullData.textures_files[i];
                
                fullData.gltf_file.textures[i].image = i;
                fullData.gltf_file.textures[i].sampler = 0;

                fullData.gltf_file.materials[i].alpha_cutoff = 0.5;
                fullData.gltf_file.materials[i].alpha_mode = materialAlphaMode::MASK;
                fullData.gltf_file.materials[i].double_sided =true;
                fullData.gltf_file.materials[i].baseColorTex.texCoord = 0;
                fullData.gltf_file.materials[i].baseColorTex.texture_index = i;
            }
        }
        int bin_file_id = 0;
        int model_n = 0;
        int max_model = MIN(settings.max_models,models.size());
        while (model_n < max_model)
        {
            int start = model_n;
            int end = model_n;
            int ind_count = 0;
            int vert_count = 0;
            bool correct_patch = true;
            while (end < max_model)
            {
                Model *m = models[end].m;
                
                if (m->positions.empty())
                {
                    logerr("glTF_writer: cannot convert model, it is empty");
                    correct_patch = false;
                }
                if (!m->indexed || m->indices.empty())
                {
                    logerr("glTF_writer: cannot convert model, only indexed models are supported");
                    correct_patch = false;
                }
                if (m->normals.size() != 0 && m->normals.size() != m->positions.size())
                {
                    logerr("glTF_writer: cannot convert model, normals does not match positions");
                    correct_patch = false;
                }
                if (MAX(sizeof(GLfloat)*m->positions.size(), sizeof(GLuint)*m->indices.size()) > settings.max_binary_file_size)
                {
                    logerr("glTF_writer: cannot convert model, model is too large");
                    correct_patch = false;
                }

                ind_count += m->indices.size();
                vert_count += m->positions.size()/3;
                if (MAX(sizeof(GLuint)*ind_count,3*sizeof(GLfloat)*vert_count) > settings.max_binary_file_size)
                {
                    break;
                }
                end++;
            }

            model_n = end;
            
            if (correct_patch)
            {
                fullData.pos_binary_files.emplace_back();
                auto &pbf = fullData.pos_binary_files.back();
                pbf.max_size = 3*sizeof(GLfloat)*vert_count;
                pbf.cur_size = 0;
                pbf.data = new char[pbf.max_size];
                pbf.file_name = asset_name + "_" + std::to_string(bin_file_id) + "_pos.bin";

                fullData.norm_binary_files.emplace_back();
                auto &nbf = fullData.norm_binary_files.back();
                nbf.max_size = 3*sizeof(GLfloat)*vert_count;
                nbf.cur_size = 0;
                nbf.data = new char[nbf.max_size];
                nbf.file_name = asset_name + "_" + std::to_string(bin_file_id) + "_norm.bin";

                fullData.tc_binary_files.emplace_back();
                auto &tbf = fullData.tc_binary_files.back();
                tbf.max_size = 2*sizeof(GLfloat)*vert_count;
                tbf.cur_size = 0;
                tbf.data = new char[tbf.max_size];
                tbf.file_name = asset_name + "_" + std::to_string(bin_file_id) + "_tc.bin";

                fullData.ind_binary_files.emplace_back();
                auto &ibf = fullData.ind_binary_files.back();
                ibf.max_size = sizeof(GLuint)*ind_count;
                ibf.cur_size = 0;
                ibf.data = new char[ibf.max_size];
                ibf.file_name = asset_name + "_" + std::to_string(bin_file_id) + "_ind.bin";

                for (int i =start; i<end; i++)
                {
                    correct_patch = correct_patch && model_to_gltf(models[i].m,fullData, bin_file_id);
                }

                if (correct_patch)
                {
                    bool b = write_to_binary_file(pbf.data, pbf.cur_size, pbf.file_name);
                    correct_patch = correct_patch && b;

                    b = write_to_binary_file(nbf.data, nbf.cur_size, nbf.file_name);
                    correct_patch = correct_patch && b;

                    b = write_to_binary_file(tbf.data, tbf.cur_size, tbf.file_name);
                    correct_patch = correct_patch && b;

                    b = write_to_binary_file(ibf.data, ibf.cur_size, ibf.file_name);
                    correct_patch = correct_patch && b;
                }
                delete[] pbf.data;
                delete[] nbf.data;
                delete[] tbf.data;
                delete[] ibf.data;

                //create buffers and buffer views for patch
                if (correct_patch && settings.debug)
                {
                    debug("glTF_writer: successfully saved models binary data. Models %d - %d. %d verts, %d ind\n",
                          start, end - 1, vert_count, ind_count);
                    int pos_bv;
                    int norm_bv;
                    int ind_bv;
                    int tc_bv;

                    //positions
                    fullData.gltf_file.buffers.emplace_back();
                    fullData.gltf_file.buffers.back().data = &pbf;
                    fullData.gltf_file.buffers.back().byte_length = pbf.cur_size;

                    fullData.gltf_file.buffer_views.emplace_back();    
                    fullData.gltf_file.buffer_views.back().buffer = fullData.gltf_file.buffers.size() - 1;
                    fullData.gltf_file.buffer_views.back().byte_offset = 0;
                    fullData.gltf_file.buffer_views.back().byte_stride = 3*sizeof(float);
                    fullData.gltf_file.buffer_views.back().byte_length = fullData.gltf_file.buffers.back().byte_length;
                    fullData.gltf_file.buffer_views.back().target = BufferViewTargetType::ARRAY_BUFFER;  
                    pos_bv = fullData.gltf_file.buffer_views.size() - 1;
                    
                    //normals
                    fullData.gltf_file.buffers.emplace_back();
                    fullData.gltf_file.buffers.back().data = &nbf;
                    fullData.gltf_file.buffers.back().byte_length = nbf.cur_size;

                    fullData.gltf_file.buffer_views.emplace_back();    
                    fullData.gltf_file.buffer_views.back().buffer = fullData.gltf_file.buffers.size() - 1;
                    fullData.gltf_file.buffer_views.back().byte_offset = 0;
                    fullData.gltf_file.buffer_views.back().byte_stride = 3*sizeof(float);
                    fullData.gltf_file.buffer_views.back().byte_length = fullData.gltf_file.buffers.back().byte_length;
                    fullData.gltf_file.buffer_views.back().target = BufferViewTargetType::ARRAY_BUFFER;
                    norm_bv = fullData.gltf_file.buffer_views.size() - 1;
                                        
                    //texture coordinates
                    fullData.gltf_file.buffers.emplace_back();
                    fullData.gltf_file.buffers.back().data = &tbf;
                    fullData.gltf_file.buffers.back().byte_length = tbf.cur_size;

                    fullData.gltf_file.buffer_views.emplace_back();    
                    fullData.gltf_file.buffer_views.back().buffer = fullData.gltf_file.buffers.size() - 1;
                    fullData.gltf_file.buffer_views.back().byte_offset = 0;
                    fullData.gltf_file.buffer_views.back().byte_stride = 2*sizeof(float);
                    fullData.gltf_file.buffer_views.back().byte_length = fullData.gltf_file.buffers.back().byte_length;
                    fullData.gltf_file.buffer_views.back().target = BufferViewTargetType::ARRAY_BUFFER;
                    tc_bv = fullData.gltf_file.buffer_views.size() - 1;

                    //indices
                    fullData.gltf_file.buffers.emplace_back();
                    fullData.gltf_file.buffers.back().data = &ibf;
                    fullData.gltf_file.buffers.back().byte_length = ibf.cur_size;

                    fullData.gltf_file.buffer_views.emplace_back();    
                    fullData.gltf_file.buffer_views.back().buffer = fullData.gltf_file.buffers.size() - 1;
                    fullData.gltf_file.buffer_views.back().byte_offset = 0;
                    fullData.gltf_file.buffer_views.back().byte_length = fullData.gltf_file.buffers.back().byte_length;
                    fullData.gltf_file.buffer_views.back().target = BufferViewTargetType::ELEMENT_ARRAY_BUFFER;
                    ind_bv = fullData.gltf_file.buffer_views.size() - 1;
                    //create accessors and meshes for each model

                    int ind_byte_offset = 0;
                    int pos_byte_offset = 0;
                    int tc_byte_offset = 0;
                    for (int i =start; i<end; i++)
                    {
                        int verts = models[i].m->positions.size()/3;
                        int ind_acc_n, pos_acc_n, norm_acc_n, tc_acc_n;
                        
                        glm::vec3 max_bounds = glm::vec3(settings.max_bound);
                        glm::vec3 min_bounds = glm::vec3(-settings.max_bound);
                        if (settings.calc_exact_bbox && verts > 0)
                        {
                            max_bounds = glm::vec3(-1e6);
                            min_bounds = glm::vec3(1e6);

                            for (int j = 0;j<models[i].m->positions.size(); j+=3)
                            {
                                max_bounds.x = MAX(max_bounds.x, models[i].m->positions[j]);
                                max_bounds.y = MAX(max_bounds.y, models[i].m->positions[j+1]);
                                max_bounds.z = MAX(max_bounds.z, models[i].m->positions[j+2]);
                                
                                min_bounds.x = MIN(min_bounds.x, models[i].m->positions[j]);
                                min_bounds.y = MIN(min_bounds.y, models[i].m->positions[j+1]);
                                min_bounds.z = MIN(min_bounds.z, models[i].m->positions[j+2]);
                            }
                        }
                        //positions accessor
                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &pos_acc = fullData.gltf_file.accessors.back();
                        pos_acc.buffer_view = pos_bv;
                        pos_acc.byte_offset = pos_byte_offset;
                        pos_acc.count = verts;
                        pos_acc.componentType = AccessorComponentType::FLOAT;
                        pos_acc.type = AccessorType::VEC3;
                        pos_acc.max_values = {max_bounds.x, max_bounds.y, max_bounds.z};
                        pos_acc.min_values = {min_bounds.x, min_bounds.y, min_bounds.z};
                        pos_acc_n = fullData.gltf_file.accessors.size() - 1;
                        
                        //normals accessor
                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &norm_acc = fullData.gltf_file.accessors.back();
                        norm_acc.buffer_view = norm_bv;
                        norm_acc.byte_offset = pos_byte_offset;
                        norm_acc.count = verts;
                        norm_acc.componentType = AccessorComponentType::FLOAT;
                        norm_acc.type = AccessorType::VEC3;
                        norm_acc_n = fullData.gltf_file.accessors.size() - 1;
                        
                        //tc accessor
                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &tc_acc = fullData.gltf_file.accessors.back();
                        tc_acc.buffer_view = tc_bv;
                        tc_acc.byte_offset = tc_byte_offset;
                        tc_acc.count = verts;
                        tc_acc.componentType = AccessorComponentType::FLOAT;
                        tc_acc.type = AccessorType::VEC2;
                        tc_acc_n = fullData.gltf_file.accessors.size() - 1;

                        //indices accessor
                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &ind_acc = fullData.gltf_file.accessors.back();
                        ind_acc.buffer_view = ind_bv;
                        ind_acc.byte_offset = ind_byte_offset;
                        ind_acc.count = models[i].m->indices.size();
                        ind_acc.componentType = AccessorComponentType::UNSIGNED_INT;
                        ind_acc.type = AccessorType::SCALAR;
                        ind_acc.min_values = {0};
                        ind_acc.max_values = {(float)(verts - 1)};
                        ind_acc_n = fullData.gltf_file.accessors.size() - 1;
                        

                        fullData.gltf_file.meshes.emplace_back();
                        fullData.gltf_file.meshes.back().primitives.emplace_back();
                        auto &pr = fullData.gltf_file.meshes.back().primitives.back();
                        pr.indicies = ind_acc_n;
                        pr.attributes.emplace(primitiveAttributeType::POSITION,pos_acc_n);
                        pr.attributes.emplace(primitiveAttributeType::NORMAL,norm_acc_n);
                        pr.attributes.emplace(primitiveAttributeType::TEXCOORD_0,tc_acc_n);
                        
                        if (models[i].material_id >= 0)
                            pr.material = models[i].material_id;
                        
                        ind_byte_offset += sizeof(uint)*models[i].m->indices.size();
                        pos_byte_offset += 3*sizeof(float)*verts;
                        tc_byte_offset += 2*sizeof(float)*verts;
                    }
                }
                bin_file_id++;
            }
        }

        //create a node for each set of transforms
        for (int i=0;i<models.size();i++)
        {
            auto &t = models[i].transforms;
            if (t.size() == 0)
                continue;
            int mesh_id = i;
            fullData.gltf_file.nodes.emplace_back();
            Node &n = fullData.gltf_file.nodes.back();
            scene.nodes.push_back(fullData.gltf_file.nodes.size() - 1);

            if (t.size() == 1)
            {
                n.mesh = mesh_id;
                n.transform = t[0];
            }
            else
            {
                n.child_nodes = {};
                n.transform = glm::mat4(1.0f);
                n.scale = glm::vec3(1,1,1);
                int nn = fullData.gltf_file.nodes.size() - 1;
                for (auto &tr : t)
                {
                    fullData.gltf_file.nodes.emplace_back();
                    fullData.gltf_file.nodes[nn].child_nodes.push_back(fullData.gltf_file.nodes.size() - 1);
                    fullData.gltf_file.nodes.back().mesh = mesh_id;
                    fullData.gltf_file.nodes.back().transform = tr;
                }
            }
        }

        int camera_id = 0;
        for (Camera *c : cameras)
        {
            if (camera_to_gltf(c,fullData, camera_id))
                camera_id++;
        }

        GltfStructureWriter gsw;
        gsw.write_to_json(fullData,name);
    }
    bool GeneralGltfWriter::model_to_gltf(Model *m, FullData &full_data, int bin_file_id)
    {
        bool ok = true;

        const char *pos_data = reinterpret_cast<const char *>(m->positions.data());
        ok = ok && add_to_binary_file(pos_data, m->positions.size()*sizeof(GLfloat), full_data.pos_binary_files[bin_file_id]);

        const char *ind_data = reinterpret_cast<const char *>(m->indices.data());
        ok = ok && add_to_binary_file(ind_data, m->indices.size()*sizeof(GLuint), full_data.ind_binary_files[bin_file_id]);

        const char *norm_data = reinterpret_cast<const char *>(m->normals.data());
        ok = ok && add_to_binary_file(norm_data, m->normals.size()*sizeof(GLfloat), full_data.norm_binary_files[bin_file_id]);

        float *tc_vec2 = safe_new<float>(m->colors.size()/2, "tc_vec2");
        for (int i=0;i<m->colors.size();i+=4)
        {
            tc_vec2[i/2] = m->colors[i];
            tc_vec2[i/2 + 1] = m->colors[i+1];
        }
        const char *tc_data = reinterpret_cast<const char *>(tc_vec2);
        ok = ok && add_to_binary_file(tc_data, (m->colors.size()/2)*sizeof(GLfloat), full_data.tc_binary_files[bin_file_id]);
        safe_delete(tc_vec2, "tc_vec2");

        return ok;
    }
    bool GeneralGltfWriter::camera_to_gltf(Camera *c, FullData &full_data, int id)
    {

    }
    bool GeneralGltfWriter::add_to_binary_file(const char *data, int size, BinaryFile &b_file)
    {
        if (size + b_file.cur_size <= b_file.max_size)
        {
            memcpy(b_file.data + b_file.cur_size, data, size);
            b_file.cur_size += size;
            return true;
        }
        return false;
    }
    bool GeneralGltfWriter::write_to_binary_file(const char *data, int size, std::string file_name)
    {
        std::ofstream fs(file_name, std::ios::out | std::ios::binary); 
        fs.write(data, size);
        fs.close();
        return !fs.fail();
    }
}