#pragma once
#include "cities_generator/global.h"
#include "ShaderDataHolders/RenderableObjectDataHolder.h"

class Skybox : public RenderableObjectDataHolder
{
    public:
        Skybox(std::string texPath, glm::vec3 sunPos, float rotation = 0.f);
        virtual void Render(RENDER_MODE mode = RENDER_MODE::DEFUALT);
        glm::vec3 GetSunPosition();

    private:
        GLuint cubeMap;
        float rotationRadians;
        glm::mat4 rotaionMatrix;
        glm::vec3 sunPos;
};