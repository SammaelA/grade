#include "cities_generator/ImageLoader.h"
#include "third_party/stb_image.h"

GLuint ImageLoader::StbLoadImage(std::string path)
{
    GLuint tex;

    int imageW, imageH;
    unsigned char* image = stbi_load(path.c_str(), &imageW, &imageH, 0, 4);

    if (!image)
    {
        debug("INCORRECT PATH: ", path);
        throw std::exception{};
    }

    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    glGenerateMipmap(GL_TEXTURE_2D);
    stbi_image_free(image);
    glBindTexture(GL_TEXTURE_2D, 0);

    return tex;
}