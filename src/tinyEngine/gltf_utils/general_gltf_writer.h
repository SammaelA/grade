#pragma once

#include "gltf_structure.h"
#include "../utility/model.h"
#include "../camera.h"
#include <cstring>
#include <fstream>
#include <iostream>
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
        bool debug = true;
    } settings;
    void add_model(Model *m) {models.push_back(m);}
    void add_camera(Camera *c) {cameras.push_back(c);}
    void clear();
    void convert_to_gltf(std::string name);
private:
    bool model_to_gltf(Model *m, FullData &full_data, int id);
    bool camera_to_gltf(Camera *c, FullData &full_data, int id);
    bool write_to_binary_file(const char *data, int size, std::string file_name);
    bool add_to_binary_file(const char *data, int size, BinaryFile &b_file);
    std::vector<Model *> models;
    std::vector<Camera *> cameras;
    std::string asset_name;
};
}