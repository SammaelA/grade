#pragma once

#include "gltf_structure.h"
#include "gltf_structure_writer.h"
#include "tinyEngine/model.h"
#include "tinyEngine/texture.h"
#include "tinyEngine/camera.h"
#include "visualizer.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
namespace gltf
{
class GeneralGltfWriter
{
public:
    struct Settings
    {
        unsigned max_binary_file_size = 64 * (2 << 20);//64 Mb
        unsigned max_models = 8192;
        unsigned max_cameras = 16;
        float max_bound = 25000;
        bool debug = true;
        bool calc_exact_bbox = true;
    } settings;
    struct ModelData
    {
        Model *m;
        ::Texture *t;
        std::vector<glm::mat4> transforms;
        int material_id = -1;
    };
    void add_model(Model *m);
    void add_packed_grove(GrovePacked &grove, GroveGenerationData &ggd);
    void add_camera(Camera *c) {cameras.push_back(c);}
    void clear();
    void convert_to_gltf(std::string name);
private:
    bool model_to_gltf(Model *m, FullData &full_data, int id);
    bool camera_to_gltf(Camera *c, FullData &full_data, int id);
    bool write_to_binary_file(const char *data, int size, std::string file_name);
    bool add_to_binary_file(const char *data, int size, BinaryFile &b_file);

    std::vector<ModelData> models;
    std::vector<Camera *> cameras;
    std::vector<Model *> temp_models;
    std::string asset_name;
};
}