#define STB_IMAGE_IMPLEMENTATION
#include "cubemap.h"
#include "graphics_utils/texture_manager.h"
#include "image.h"

Cubemap::Cubemap(int w, int h):
cube_model(),
cube_shader({"cubemap.vs", "cubemap.fs"}, {"in_Position"})
{
    width = w;
    height = h;
    /*
    std::vector<std::string> faces
    {
        "Skybox/hills_bk.tga",
        "Skybox/hills_dn.tga",
        "Skybox/hills_ft.tga",
        "Skybox/hills_lf.tga",
        "Skybox/hills_rt.tga",
        "Skybox/Untitled.png"
    };*/
    std::vector<std::string> faces
    {
        "clouds1/clouds1_south.bmp",
        "clouds1/clouds1_down.bmp",
        "clouds1/clouds1_north.bmp",
        "clouds1/clouds1_east.bmp",
        "clouds1/clouds1_west.bmp",
        "clouds1/clouds1_up.bmp"
    };
    loadCubemap(faces);
    
    cubeFBO = create_framebuffer();

    tex = textureManager.create_texture(w, h, GL_RGB8, 1, nullptr, GL_RGB, GL_UNSIGNED_BYTE);

    glBindFramebuffer(GL_FRAMEBUFFER, cubeFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex.texture, 0);

}
void Cubemap::render(glm::mat4 &projection, glm::mat4 &view, Camera &camera)
{
    glBindFramebuffer(GL_FRAMEBUFFER, cubeFBO);
    glViewport(0, 0, width, height);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    view = glm::mat4(glm::mat3(view));
    glDepthMask(GL_FALSE);
    cube_shader.use();
    cube_shader.textureCube("skybox",cube.texture);
    cube_shader.uniform("projection",projection);
    cube_shader.uniform("view",view);
    glBindVertexArray(cube_model.vao);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glDepthMask(GL_TRUE);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
Cubemap::~Cubemap()
{
  textureManager.delete_tex(cube);
  textureManager.delete_tex(tex);
  delete_framebuffer(cubeFBO);
}
void Cubemap::loadCubemap(std::vector<std::string> faces)
{
    int width, height, nrChannels;
    int nums[6] = {3,4,5,1,2,0};
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        std::string path = image::base_img_path + faces[nums[i]];
        unsigned char *data = stbi_load(path.c_str(), &width, &height, &nrChannels, 0);
        if (i ==0)
        {
          cube = textureManager.create_texture_cube(width, height, GL_RGB8, 1);
          glBindTexture(GL_TEXTURE_CUBE_MAP, cube.texture);
        }
        if (data)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 
                         0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data
            );
            stbi_image_free(data);
        }
        else
        {
            std::cout << "Cubemap texture failed to load at path: " << path << std::endl;
            stbi_image_free(data);
        }
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
}