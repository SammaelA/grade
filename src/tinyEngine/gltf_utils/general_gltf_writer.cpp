#include "general_gltf_writer.h"
#include "../utility.h"

namespace gltf
{
    void GeneralGltfWriter::clear()
    {
        models.clear();
        cameras.clear();
        transforms.clear();
    }
    void GeneralGltfWriter::add_model(Model *m)
    {
        models.push_back(m);
        transforms.push_back(std::pair<int, std::vector<glm::mat4>>(models.size() - 1, {m->model}));
    }
    void GeneralGltfWriter::convert_to_gltf(std::string name)
    {
        FullData fullData;
        asset_name = name;

        fullData.gltf_file.main_scene = 0;
        fullData.gltf_file.scenes.emplace_back();

        auto &scene = fullData.gltf_file.scenes[fullData.gltf_file.main_scene];

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
                Model *m = models[end];
                
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

                fullData.ind_binary_files.emplace_back();
                auto &ibf = fullData.ind_binary_files.back();
                ibf.max_size = sizeof(GLuint)*ind_count;
                ibf.cur_size = 0;
                ibf.data = new char[ibf.max_size];
                ibf.file_name = asset_name + "_" + std::to_string(bin_file_id) + "_ind.bin";

                for (int i =start; i<end; i++)
                {
                    correct_patch = correct_patch && model_to_gltf(models[i],fullData, bin_file_id);
                }

                if (correct_patch)
                {
                    bool b = write_to_binary_file(pbf.data, pbf.cur_size, pbf.file_name);
                    correct_patch = correct_patch && b;

                    b = write_to_binary_file(nbf.data, nbf.cur_size, nbf.file_name);
                    correct_patch = correct_patch && b;

                    b = write_to_binary_file(ibf.data, ibf.cur_size, ibf.file_name);
                    int *data = reinterpret_cast<int *>(ibf.data);
                    logerr("max = %d",ibf.max_size);
                    for (int i=0;i<72;i++)
                    {
                        logerr("[%d] %d",i,data[i]);
                    }
                    correct_patch = correct_patch && b;
                }
                delete[] pbf.data;
                delete[] nbf.data;
                delete[] ibf.data;

                //create buffers and buffer views for patch
                if (correct_patch && settings.debug)
                {
                    debug("glTF_writer: successfully saved models binary data. Models %d - %d. %d verts, %d ind\n",
                          start, end - 1, vert_count, ind_count);
                    int pos_bv;
                    int norm_bv;
                    int ind_bv;

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
                    /*
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
                    */

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
                    for (int i =start; i<end; i++)
                    {
                        int verts = models[i]->positions.size()/3;
                        int ind_acc_n, pos_acc_n, norm_acc_n;

                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &pos_acc = fullData.gltf_file.accessors.back();
                        pos_acc.buffer_view = pos_bv;
                        pos_acc.byte_offset = pos_byte_offset;
                        pos_acc.count = verts;
                        pos_acc.componentType = AccessorComponentType::FLOAT;
                        pos_acc.type = AccessorType::VEC3;
                        pos_acc.min_values = {-settings.max_bound, -settings.max_bound, -settings.max_bound};
                        pos_acc.max_values = {settings.max_bound, settings.max_bound, settings.max_bound};
                        pos_acc_n = fullData.gltf_file.accessors.size() - 1;
                        /*
                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &norm_acc = fullData.gltf_file.accessors.back();
                        norm_acc.buffer_view = norm_bv;
                        norm_acc.byte_offset = pos_byte_offset;
                        norm_acc.count = verts;
                        norm_acc.componentType = AccessorComponentType::FLOAT;
                        norm_acc.type = AccessorType::VEC3;
                        norm_acc_n = fullData.gltf_file.accessors.size() - 1;
                        */
                        fullData.gltf_file.accessors.emplace_back();
                        Accessor &ind_acc = fullData.gltf_file.accessors.back();
                        ind_acc.buffer_view = ind_bv;
                        ind_acc.byte_offset = ind_byte_offset;
                        ind_acc.count = models[i]->indices.size();
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
                        //pr.attributes.emplace(primitiveAttributeType::NORMAL,norm_acc_n);
                        
                        ind_byte_offset += sizeof(uint)*models[i]->indices.size();
                        pos_byte_offset += 3*sizeof(float)*verts;
                    }
                }
                bin_file_id++;
            }
        }

        //create a node for each set of transforms
        for (auto &t :transforms)
        {
            if (t.second.size() == 0)
                continue;
            int mesh_id = t.first;
            fullData.gltf_file.nodes.emplace_back();
            Node &n = fullData.gltf_file.nodes.back();
            scene.nodes.push_back(fullData.gltf_file.nodes.size() - 1);

            if (t.second.size() == 1)
            {
                n.mesh = mesh_id;
                n.transform = t.second[0];
            }
            else
            {
                n.child_nodes = {};
                int nn = fullData.gltf_file.nodes.size() - 1;
                for (auto &tr : t.second)
                {
                    fullData.gltf_file.nodes.emplace_back();
                    logerr("aa %d",n.child_nodes.size());
                    fullData.gltf_file.nodes[nn].child_nodes.push_back(fullData.gltf_file.nodes.size() - 1);
                    logerr("bb");
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

        /*scene.nodes.push_back(0);

        fullData.gltf_file.nodes.emplace_back();
        fullData.gltf_file.nodes.back().mesh = 0;

        fullData.gltf_file.meshes.emplace_back();
        fullData.gltf_file.meshes.back().primitives.emplace_back();
        auto &p = fullData.gltf_file.meshes.back().primitives.back();
        p.attributes.emplace(primitiveAttributeType::POSITION,1);
        p.indicies = 0;

        fullData.gltf_file.buffers.emplace_back();
        fullData.gltf_file.buffers.back().data = &(fullData.pos_binary_files[0]);
        fullData.gltf_file.buffers.back().byte_length = fullData.pos_binary_files[0].cur_size;
        fullData.gltf_file.buffers.emplace_back();
        fullData.gltf_file.buffers.back().data = &(fullData.ind_binary_files[0]);
        fullData.gltf_file.buffers.back().byte_length = fullData.ind_binary_files[0].cur_size;

        fullData.gltf_file.buffer_views.emplace_back();
        fullData.gltf_file.buffer_views.back().buffer = 0;
        fullData.gltf_file.buffer_views.back().byte_offset = 0;
        fullData.gltf_file.buffer_views.back().byte_length = 96;
        fullData.gltf_file.buffer_views.back().target = BufferViewTargetType::ELEMENT_ARRAY_BUFFER;
        fullData.gltf_file.buffer_views.emplace_back();
        fullData.gltf_file.buffer_views.back().buffer = 1;
        fullData.gltf_file.buffer_views.back().byte_offset = 0;
        fullData.gltf_file.buffer_views.back().byte_length = 144;
        fullData.gltf_file.buffer_views.back().target = BufferViewTargetType::ARRAY_BUFFER;

        fullData.gltf_file.accessors.emplace_back();
        fullData.gltf_file.accessors.back().buffer_view = 1;
        fullData.gltf_file.accessors.back().byte_offset = 0;
        fullData.gltf_file.accessors.back().componentType = AccessorComponentType::UNSIGNED_INT;
        fullData.gltf_file.accessors.back().count = 36;
        fullData.gltf_file.accessors.back().type = AccessorType::SCALAR;
        fullData.gltf_file.accessors.back().max_values = {7};
        fullData.gltf_file.accessors.back().min_values = {0};
        fullData.gltf_file.accessors.emplace_back();
        fullData.gltf_file.accessors.back().buffer_view = 0;
        fullData.gltf_file.accessors.back().byte_offset = 0;
        fullData.gltf_file.accessors.back().componentType = AccessorComponentType::FLOAT;
        fullData.gltf_file.accessors.back().count = 8;
        fullData.gltf_file.accessors.back().type = AccessorType::VEC3;*/

        GltfStructureWriter gsw;
        gsw.write_to_json(fullData,"test");
    }
    bool GeneralGltfWriter::model_to_gltf(Model *m, FullData &full_data, int bin_file_id)
    {
        logerr("model to gltf");
        bool ok = true;

        const char *pos_data = reinterpret_cast<const char *>(m->positions.data());
        ok = ok && add_to_binary_file(pos_data, m->positions.size()*sizeof(GLfloat), full_data.pos_binary_files[bin_file_id]);
        //std::string pos_bin_name = asset_name + "_" + std::to_string(id) + "_pos.bin";
        //ok = ok && write_to_binary_file(pos_data, m->positions.size()*sizeof(GLfloat), pos_bin_name);

        const char *ind_data = reinterpret_cast<const char *>(m->indices.data());
        ok = ok && add_to_binary_file(ind_data, m->indices.size()*sizeof(GLuint), full_data.ind_binary_files[bin_file_id]);

        const char *norm_data = reinterpret_cast<const char *>(m->normals.data());
        ok = ok && add_to_binary_file(norm_data, m->normals.size()*sizeof(GLfloat), full_data.norm_binary_files[bin_file_id]);

        /*float *tc_vec2 = new float[m->colors.size()/2];
        for (int i=0;i<m->colors.size();i+=4)
        {
            tc_vec2[i/2] = m->colors[i];
            tc_vec2[i/2 + 1] = m->colors[i+1];
        }*/

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