#include "grass_renderer.h"
#include "texture_manager.h"
#include "camera.h"
GrassRenderer::GrassRenderer():
 grass({"grass.vs", "grass.fs"}, {"in_Position","in_Normal", "in_Tex"}),
 grassShadow({"grass.vs", "depth_billboard.fs"}, {"in_Position","in_Normal", "in_Tex"}),
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
                           HeightmapTex &heightmap_tex, DirectedLight &light, bool to_shadow)
{
    Shader &shader = to_shadow ? grassShadow : grass;
    shader.use();
    shader.uniform("projectionCamera", prc);
    shader.texture("hmap",heightmap_tex.get());
    shader.texture("perlin",perlin);
    shader.texture("noise",noise);
    shader.texture("shadowMap",shadow_tex);
    shader.uniform("sts_inv",1.0f/light.shadow_map_size);
    shader.uniform("need_shadow",shadow_tex != 0);
    shader.uniform("lightSpaceMatrix",shadow_tr);
    shader.uniform("dir_to_sun", light.dir);
    shader.uniform("light_color", light.color*light.intensity);
    shader.uniform("camera_pos", camera_pos);
    shader.uniform("ambient_diffuse_specular", glm::vec3(light.ambient_q,light.diffuse_q,light.specular_q));
    glBindVertexArray(m.vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ibo);

    shader.texture("tex",grass_base);
    glDrawElementsInstanced(GL_TRIANGLES, m.SIZE, GL_UNSIGNED_INT, 0, 40000);
}