#include "BasicTexturedObjectDataHolder.h"
#include "RenderableObjectDataHolder.h"

BasicTexturedObjectDataHolder::~BasicTexturedObjectDataHolder()
{
    glDeleteTextures(1, &diffuseTexture);
}

BasicTexturedObjectDataHolder::BasicTexturedObjectDataHolder()
{
    diffuseTexture = BAD_TEXTURE;
    diffuseColor = DEFAULT_COLOR;
}

void BasicTexturedObjectDataHolder::LoadInShader_BasicTexturedObject(int& texCount)
{
    VerifyUnoformPos(isUsingDiffuseTexture_Unifrom, "is_isotropic_color");
    VerifyUnoformPos(diffuseColor_Unifrom, "diffuse_color");
    VerifyUnoformPos(diffuseTexture_Unifrom, "diffuse_texture");
    VerifyUnoformPos(roughnessFactor_Unifrom, "roughness_factor");

    glUniform1i(isUsingDiffuseTexture_Unifrom, !isUsingDiffuseTexture);
    glUniform4fv(diffuseColor_Unifrom, 1, glm::value_ptr(diffuseColor));
    glUniform1f(roughnessFactor_Unifrom, roughnessFactor);
    glActiveTexture(GL_TEXTURE0 + texCount);
    glBindTexture(GL_TEXTURE_2D, diffuseTexture);
    glUniform1i(diffuseTexture_Unifrom, texCount);
    texCount++;
}

void BasicTexturedObjectDataHolder::SetColor(glm::vec4 col)
{
    diffuseColor = col;
}