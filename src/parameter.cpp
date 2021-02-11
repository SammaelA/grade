#include "parameter.h"
Parameter<int> TreeStructureParameters::from_float(Parameter<float> source)
{
    std::vector<int> a;
    for (int i=0;i<source.stateParams.size();i++)
        a.push_back(source.stateParams[i]);
    
    Parameter<int> par(source.baseValue,a,source.randomnessLevel,source.randomizer,source.minValue,source.maxValue);
    return par;
}