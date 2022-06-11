#pragma once
#include "cities_generator/global.h"

class BaseDataHolder
{
    protected:
        void VerifyUnoformPos(int& pos, std::string shaderVariableName);
        virtual ~BaseDataHolder() = default;  
};
