#include "grass_renderer.h"
#include "graphics_utils/texture_manager.h"
#include "tinyEngine/camera.h"
#include "graphics_utils/modeling.h"

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

void GrassRenderer::render(glm::mat4 &projection, glm::mat4 &view, glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                           HeightmapTex &heightmap_tex, DirectedLight &light, bool to_shadow)
{
    Shader &shader = to_shadow ? grassShadow : grass;
    shader.use();
    shader.uniform("projection", projection);
    shader.uniform("view", view);
    shader.texture("hmap",heightmap_tex.get());
    shader.texture("perlin",perlin);
    shader.texture("noise",noise);
    if (to_shadow)
        shader.uniform("opaqueness",0.4f);
    glBindVertexArray(m.vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ibo);

    shader.texture("tex",grass_base);
    glDrawElementsInstanced(GL_TRIANGLES, m.SIZE, GL_UNSIGNED_INT, 0, 40000);
}

GrassRenderer2::GrassRenderer2(GrassPacked &data):
grass_atlas(data.grass_textures),
grass({"grass2.vs", "grass2.fs"}, {"in_Position","in_Normal", "in_Tex"}),
grassShadow({"grass2.vs", "grass2_shadow.fs"}, {"in_Position","in_Normal", "in_Tex"})
{
    glm::vec4 tex_transform = glm::vec4(1,1,0,0);

    int total_instances = 0;
    Texture null = textureManager.empty();
    int type_n = 0;
    for (auto &p : data.grass_instances)
    {
        tex_transform = data.grass_textures.tc_transform(p.first);
        models.push_back(model_loader::create_model_by_name(data.used_grass_types[type_n].model_name,null));
        for (int i=0;i<models.back()->colors.size();i+=4)
        {
            models.back()->colors[i] = tex_transform.x*(models.back()->colors[i] + tex_transform.z);
            models.back()->colors[i + 1] = tex_transform.y*(models.back()->colors[i + 1] + tex_transform.w);
        }
        models.back()->update();
        inst_offsets.push_back(total_instances);
        inst_counts.push_back(p.second.size());
        total_instances += p.second.size();
        type_n++;
    }

    glm::mat4 *matrices = new glm::mat4[total_instances];
    int i=0;
    for (auto &p : data.grass_instances)
    {
        for (auto &in : p.second)
        {
            matrices[i] = glm::scale(
                          glm::rotate(glm::translate(glm::mat4(1.0f),glm::vec3(in.pos)),
                                      in.rot_y,glm::vec3(0,1,0)), glm::vec3(in.size));
            i++;
        }
    }
    instances_buffer = create_buffer();
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 11, instances_buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::mat4)*total_instances, matrices, GL_STATIC_DRAW);

    delete[] matrices;
}

GrassRenderer2::~GrassRenderer2()
{
    delete_buffer(instances_buffer);
    for (auto *m : models)
        delete m;
}

void GrassRenderer2::render(glm::mat4 &projection, glm::mat4 &view, glm::mat4 &shadow_tr, GLuint shadow_tex, glm::vec3 camera_pos,
                            HeightmapTex &heightmap_tex, DirectedLight &light, bool to_shadow)
{
    Shader &shader = to_shadow ? grassShadow : grass;
    shader.use();
    shader.uniform("projection", projection);
    shader.uniform("view", view);
    if (to_shadow)
        shader.uniform("opaqueness",0.67f);
    shader.texture("tex",grass_atlas.tex(0));
    for (int i=0;i<models.size();i++)
    {
        shader.uniform("inst_buf_offset",inst_offsets[i]);

        glBindVertexArray(models[i]->vao);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, models[i]->ibo);
        //logerr("set IBO %d size %d",models[i]->ibo,models[i]->SIZE);
        //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 11, instances_buffer);
        glDrawElementsInstanced(GL_TRIANGLES, models[i]->SIZE, GL_UNSIGNED_INT, 0, inst_counts[i]);
    }
}