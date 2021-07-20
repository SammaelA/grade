#pragma once
#include <functional>
#include "parameter.h"
#include "grove.h"
enum SelectionType
{
    BruteForce
};
enum Metric
{
    CompressionRatio
};
class ParameterSelector
{
public:
    ParameterSelector(std::function<void(TreeStructureParameters &, GrovePacked &)> &_generate):
    generate(_generate)
    {

    }
    void select(TreeStructureParameters &param, SelectionType sel_type, Metric metric);
private:
    std::function<void(TreeStructureParameters &, GrovePacked &)> &generate;
};