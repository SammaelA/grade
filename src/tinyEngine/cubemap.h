#pragma once
#include "third_party/stb_image.h"
#include "texture.h"
#include "shader.h"
#include <vector>
#include <string>
#include "camera.h"
#include "model.h"
class Cubemap
{
public:
    Cubemap(int w = 100, int h = 100);
    ~Cubemap();
    void loadCubemap(std::vector<std::string> faces);
    Texture &get_tex(){ return tex;}
    void render(glm::mat4 &projection, glm::mat4 &view, Camera &camera);
private:
    int width, height;
    GLuint cubeFBO;
    Texture tex, cube;
    Shader cube_shader;
    Cube cube_model;
};