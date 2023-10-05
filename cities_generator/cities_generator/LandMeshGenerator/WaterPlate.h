#pragma once
#include "cities_generator/global.h"
#include "../ShaderDataHolders/RenderableObjectDataHolder.h"
#include "cities_generator/ECSClass.h"

class WaterPlate : public RenderableObjectDataHolder, public Scene3DObjectDataHolder
{
    public:
        WaterPlate(glm::vec3 pos, glm::vec2 size);
        virtual void Render(RENDER_MODE mode);
        ~WaterPlate();    

        GLuint FBO;
        GLuint aboveWaterTexture, underWaterTexture, depthTexture, bloomTexture;

    private:
        Entity* entityPtr;
        GLuint DUDVTexture, normalMapTexture;
        float timer;
};