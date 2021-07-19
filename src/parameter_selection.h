#pragma once
#include <functional>
#include "parameter.h"
class ParameterSelector
{
public:
    ParameterSelector(std::function<void(TreeStructureParameters &)> &generate);
private:
    std::function<void(TreeStructureParameters &)> &generate;
};