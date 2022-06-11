#pragma once
#include "cities_generator/global.h"
#include "DataHolderGlobal.h"
#include "BaseDataHolder.h"

class Scene3DObjectDataHolder : public BaseDataHolder
{
    public:
        glm::vec3 position;
        glm::vec3 scale;
        glm::mat4 rotation;
        
        glm::mat4 GetMatrix();
        glm::mat4 GetViewMatrixFromThisPoint();
        Scene3DObjectDataHolder();
        void LoadInShader_Scene3DObject();

        static glm::vec3 GetDefaultViewDirection() {return glm::vec3(1,0,0);};

        virtual ~Scene3DObjectDataHolder() = default;

    protected:
        int m_TransformUniformPos = MISSING_UNIFORM_POS;
};
