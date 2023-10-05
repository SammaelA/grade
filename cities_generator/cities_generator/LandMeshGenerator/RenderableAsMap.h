#pragma once
#include "cities_generator/global.h"

class RenderableAsMap
{
    friend class Landscape;
    protected:
        GLuint m_textureId;
        virtual bool IsTextureReady() = 0;
        virtual void SaveAsTexture() = 0;

    public:
        virtual GLuint GetTexture() = 0;
};