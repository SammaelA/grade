#include "general_gltf_writer.h"
#include "../utility.h"

namespace gltf
{
    void GeneralGltfWriter::clear()
    {
        models.clear();
        cameras.clear();
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
                    correct_patch = correct_patch && b;
                }

                if (correct_patch && settings.debug)
                {
                    debug("glTF_writer: successfully saved models binary data. Models %d - %d. %d verts, %d ind\n",
                          start, end - 1, vert_count, ind_count);
                }

                bin_file_id++;
            }
        }

        int camera_id = 0;
        for (Camera *c : cameras)
        {
            if (camera_to_gltf(c,fullData, camera_id))
                camera_id++;
        }
    }
    bool GeneralGltfWriter::model_to_gltf(Model *m, FullData &full_data, int bin_file_id)
    {
        bool ok = true;

        const char *pos_data = reinterpret_cast<const char *>(m->positions.data());
        ok = ok && add_to_binary_file(pos_data, m->positions.size()*sizeof(GLfloat), full_data.pos_binary_files[bin_file_id]);
        //std::string pos_bin_name = asset_name + "_" + std::to_string(id) + "_pos.bin";
        //ok = ok && write_to_binary_file(pos_data, m->positions.size()*sizeof(GLfloat), pos_bin_name);

        const char *ind_data = reinterpret_cast<const char *>(m->indices.data());
        ok = ok && add_to_binary_file(ind_data, m->indices.size()*sizeof(GLfloat), full_data.ind_binary_files[bin_file_id]);

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
            memcpy(b_file.data, data, size);
            b_file.cur_size += size;
            return true;
        }
        return false;
    }
    bool GeneralGltfWriter::write_to_binary_file(const char *data, int size, std::string file_name)
    {
        std::ofstream fs(file_name, std::ios::out | std::ios::binary | std::ios::app); 
        fs.write(data, size);
        fs.close();
        return !fs.fail();
    }
}