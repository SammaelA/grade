#include "grass_renderer.h"
#include "texture_manager.h"
GrassRenderer::GrassRenderer():
 shader({"grass.vs", "grass.fs"}, {"in_Position","in_Normal", "in_Tex"}),
 grass_base(textureManager.get("grass")),
 grass_tall(textureManager.get("grass")),
 perlin(textureManager.get("noise")),
 noise(textureManager.get("colored_noise"))
{
    std::vector<float> vertexes = {-0.5,0,0, -0.5,1,0, 0.5,0,0, 0.5,1,0,   0,0,-0.5, 0,1,-0.5, 0,0,0.5, 0,1,0.5};
    std::vector<float> tc = {0,1,0,0, 0,0,0,0, 1,1,0,0, 1,0,0,0, 0,1,0,0, 0,0,0,0, 1,1,0,0, 1,0,0,0};
    std::vector<float> normals = {0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,0,0, 1,0,0, 1,0,0, 1,0,0};
    std::vector<GLuint> indices = {0, 1, 3, 2, 0, 3, 4,5,7, 6,4,7};

    std::function<void(Model *)> _c_mip = [&](Model *h) 
    {
        m.positions = vertexes;
        m.colors = tc;
        m.normals = normals;
        m.indices = indices;
    };
    m.construct(_c_mip);
}

void GrassRenderer::render(glm::mat4 prc, glm::mat4 shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                           HeightmapTex &heightmap_tex, DirectedLight &light)
{
    shader.use();
    shader.uniform("projectionCamera", prc);
    shader.texture("hmap",heightmap_tex.get());
    shader.texture("perlin",perlin);
    shader.texture("noise",noise);
    shader.texture("tex1",grass_base);
    shader.texture("tex2",grass_tall);
    glBindVertexArray(m.vao);
    if(m.indexed)
    {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ibo);
      glDrawElementsInstanced(GL_TRIANGLES, m.SIZE, GL_UNSIGNED_INT, 0, 40000);
    }
}