#include "postfx.h"
PostFx::PostFx(std::string pixel_shader):
 shader({"postfx.vs", pixel_shader}, {"in_Position", "in_Tex"})
{
    std::vector<float> vertexes = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0};
    std::vector<float> tc = {0,0,0,0, 1,0,0,0, 0,1,0,0, 1,1,0,0};
    std::vector<GLuint> indices = {0, 1, 3, 2, 0, 3};

    std::function<void(Model *)> _c_mip = [&](Model *h) 
    {
        m.positions = vertexes;
        m.colors = tc;
        m.indices = indices;
    };
    m.construct(_c_mip);
}
void PostFx::use()
{
    shader.use();
}
void PostFx::render()
{
    m.render(GL_TRIANGLES);
}
PostFx::~PostFx()
{

}