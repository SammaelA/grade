#pragma once
#include "cities_generator/global.h"
#include "../Scene3DObjectDataHolder.h"
#include "../BasicTexturedObjectDataHolder.h"
#include "../RenderableObjectDataHolder.h"
#include "../DataHolderGlobal.h"
#include <vector>

class BasicGLTFModel : 
    public Scene3DObjectDataHolder,
    public BasicTexturedObjectDataHolder,
    public RenderableObjectDataHolder
{
    public:
        static BasicGLTFModel* LoadGltfModel
            (std::string path, std::string obj_name, std::string shaderName);
        static std::vector<glm::vec3> GetGltfVertexes(std::string path, std::string obj_name);
        virtual void Render(RENDER_MODE mode);

    protected:
        int viewProjMatrix_Uniform = MISSING_UNIFORM_POS;  
        int sunPos_Uniform = MISSING_UNIFORM_POS;
        int sunColor_Uniform = MISSING_UNIFORM_POS;
        int sunViewProjMatrix_Uniform = MISSING_UNIFORM_POS;
        int isDepthOnly_Uniform = MISSING_UNIFORM_POS;
        int cameraPosition_Uniform = MISSING_UNIFORM_POS;
        int shadowTexture_Unifrom = MISSING_UNIFORM_POS;
};
