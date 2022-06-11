#pragma once
#include "cities_generator/global.h"
#include "cities_generator/shaderClass.h"
#include "VertexDataHolder.h"

class RenderableObjectDataHolder : public VertexDataHolder
{
    friend class Entity;
    public:
        Shader* m_shader;
        RenderableObjectDataHolder() = default;
        RenderableObjectDataHolder(std::string shaderName);
        virtual ~RenderableObjectDataHolder() = default;
        virtual void Render(RENDER_MODE rmode);
        virtual void SetClippingSettings(glm::vec4 plane);
        void UseShader();
};
