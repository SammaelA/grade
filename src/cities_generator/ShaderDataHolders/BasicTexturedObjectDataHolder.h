#pragma once
#include "cities_generator/global.h"
#include "BaseDataHolder.h"
#include "DataHolderGlobal.h"
#include "cities_generator/Renderer.h"

class BasicTexturedObjectDataHolder : public BaseDataHolder
{
    friend class Renderer;

    public:
        void SetColor(glm::vec4);
        void LoadInShader_BasicTexturedObject(int& texCount);
        BasicTexturedObjectDataHolder();
        virtual ~BasicTexturedObjectDataHolder();

    protected:
        GLuint diffuseTexture;
        glm::vec4 diffuseColor;
        bool isUsingDiffuseTexture;
        float roughnessFactor;

        int diffuseTexture_Unifrom = MISSING_UNIFORM_POS;
        int diffuseColor_Unifrom = MISSING_UNIFORM_POS;
        int isUsingDiffuseTexture_Unifrom = MISSING_UNIFORM_POS;
        int roughnessFactor_Unifrom = MISSING_UNIFORM_POS;
};
