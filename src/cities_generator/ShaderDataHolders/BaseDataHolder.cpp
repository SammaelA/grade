#include "BaseDataHolder.h"
#include "DataHolderGlobal.h"
#include "RenderableObjectDataHolder.h"

void BaseDataHolder::VerifyUnoformPos(int& pos, std::string shaderVariableName)
{
    if (pos == MISSING_UNIFORM_POS)
    {
        if (IS_DERIVED(this, RenderableObjectDataHolder))
        {
            pos = glGetUniformLocation(
                dynamic_cast<RenderableObjectDataHolder*>(this)->m_shader->Program,
                shaderVariableName.c_str()
            );
        }
        else
        {
            debug("NO CORRECT ANCESTOR!");
            throw std::exception{};
        }
    }
}
