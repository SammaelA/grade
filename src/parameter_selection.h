#pragma once
#include <functional>
#include "parameter.h"
#include "grove.h"
enum SelectionType
{
    BruteForce,
    SimulatedAnnealing,
    PolycentricSimulatedAnnealing
};
enum MetricType
{
    Dummy,
    CompressionRatio,
    ImpostorSimilarity
};
enum SelectionSchedule
{
    AllInOne,//Metric for all possible values in set is calculated, values are sorted by metric and 
             //simulated annealing is repeated with new center in every value for best to worst
    SetbySet,//the same process independently for every set
    SetbySetRandomized,//the same process independently for every set. Values in every set are taken 
                      //in a random order
    UnitbyUnit//the same process independently for every unit in every set
};
struct SelectionUnit
{
    std::string parameter_name;
    std::vector<float> base_values_set;
};
typedef std::vector<SelectionUnit> SelectionSet;
struct SetSelectionProgram
{
    std::vector<SelectionSet> selections;
    SelectionSchedule schedule;
};
struct ExitConditions
{
    float metric_reached = 10;
    int max_iters = INT_MAX;
    float max_time_seconds = 24*60*60;
    float part_of_set_covered = 1;
};
class ParameterSelector
{
public:
    ParameterSelector(std::function<void(TreeStructureParameters &, GrovePacked &)> &_generate):
    generate(_generate)
    {

    }
    void select(TreeStructureParameters &param, SelectionType sel_type, MetricType metric);
private:
    std::function<void(TreeStructureParameters &, GrovePacked &)> &generate;
};