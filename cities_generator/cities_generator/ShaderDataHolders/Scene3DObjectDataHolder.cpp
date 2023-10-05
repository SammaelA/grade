#include "Scene3DObjectDataHolder.h"
#include "RenderableObjectDataHolder.h"

Scene3DObjectDataHolder::Scene3DObjectDataHolder()
{
    position = glm::vec3(0,0,0);
    rotation = glm::mat4(1.0);
    scale = glm::vec3(1,1,1);
}

glm::mat4 Scene3DObjectDataHolder::GetMatrix()
{
    glm::mat4 res = glm::mat4(1.0f);
    res = glm::translate(res, position);
    res *= rotation;
    res = glm::scale(res, scale);
    return res;
}

glm::mat4 Scene3DObjectDataHolder::GetViewMatrixFromThisPoint()
{
    glm::vec4 direction = glm::vec4(Scene3DObjectDataHolder::GetDefaultViewDirection(), 0) * rotation;
    return glm::lookAt(position, position + glm::vec3(direction), glm::vec3(0,1,0));
}

void Scene3DObjectDataHolder::LoadInShader_Scene3DObject()
{
    VerifyUnoformPos(m_TransformUniformPos, "transform_matrix");

    glUniformMatrix4fv(m_TransformUniformPos, 1, GL_FALSE, glm::value_ptr(GetMatrix()));
}