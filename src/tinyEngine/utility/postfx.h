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
private:
    Model m;
    Shader shader;
};