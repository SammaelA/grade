#pragma once
#include "texture.h"
#include "shader.h"
#include "model.h"
#include <map>
class PostFx
{
    using slist = std::initializer_list<std::string>;
public:
    PostFx(std::string pixel_shader);
    Shader &get_shader() {return shader;}
    void use();
    void render();
    ~PostFx();
    static glm::vec4 tc_tr_mult(glm::vec4 tc_tr_1, glm::vec4 tc_tr_2)
    {
      return glm::vec4(tc_tr_2.x + tc_tr_1.x * tc_tr_2.z,
                       tc_tr_2.y + tc_tr_1.y * tc_tr_2.w,
                       tc_tr_1.z * tc_tr_2.z,
                       tc_tr_1.w * tc_tr_2.w);
    }

private:
    Model m;
    Shader shader;
};