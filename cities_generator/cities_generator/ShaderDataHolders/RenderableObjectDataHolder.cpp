#include "RenderableObjectDataHolder.h"
#include "cities_generator/shaderClass.h"
#include "cities_generator/ECSClass.h"

RenderableObjectDataHolder::RenderableObjectDataHolder(std::string Name)
{
    m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer(Name);
}

void RenderableObjectDataHolder::UseShader()
{
    m_shader->Use();
}

void RenderableObjectDataHolder::Render(RENDER_MODE rmode)
{
    UseShader();
    glUniform1i(glGetUniformLocation(m_shader->Program, "render_mode"), rmode);
    switch (rmode)
    {
        case DEFUALT:
            VertexDataHolder::Render(rmode);
            break;

        case CLIPPING:
            glEnable(GL_CLIP_DISTANCE0); 
            VertexDataHolder::Render(rmode);
            glDisable(GL_CLIP_DISTANCE0); 
            break;

        default:
            VertexDataHolder::Render(rmode);
            break;
    };
}

void RenderableObjectDataHolder::SetClippingSettings(glm::vec4 plane)
{
    UseShader();
    glUniform4fv(glGetUniformLocation(m_shader->Program, "clipping_plane"), 1, glm::value_ptr(plane));
}