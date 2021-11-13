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
    std::vector<std::string> faces
    {
        "Skybox/hills_bk.tga",
        "Skybox/hills_dn.tga",
        "Skybox/hills_ft.tga",
        "Skybox/hills_lf.tga",
        "Skybox/hills_rt.tga",
        "Skybox/Untitled.png"
    };
    cube = loadCubemap(faces);
    
    glGenFramebuffers(1, &cubeFBO);

    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindFramebuffer(GL_FRAMEBUFFER, cubeFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);

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
    cube_shader.texture("skybox",cube);
    cube_shader.uniform("projection",projection);
    cube_shader.uniform("view",view);
    glBindVertexArray(cube_model.vao);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glDepthMask(GL_TRUE);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
Cubemap::~Cubemap()
{
    glDeleteTextures(1,&cube);
    glDeleteTextures(1,&tex);
    glDeleteFramebuffers(1, &cubeFBO);
}
unsigned int Cubemap::loadCubemap(std::vector<std::string> faces)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

    int width, height, nrChannels;
    int nums[6] = {3,4,5,1,2,0};
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        std::string path = image::base_img_path + faces[nums[i]];
        unsigned char *data = stbi_load(path.c_str(), &width, &height, &nrChannels, 0);
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

    return textureID;
}