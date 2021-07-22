#include "parameter_selection.h"
#include "metric.h"

float brute_force_selection(TreeStructureParameters &param, Metric *metric,
                            std::function<void(TreeStructureParameters &, GrovePacked &)> &generate)
{
    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param.get_mask_and_data(mask, data);
    GrovePacked tree;
    float max_metr = 0;

    for (int i = 0; i < 10; i++)
    {
        GrovePacked tree = GrovePacked();
        data[0]--;
        param.load_from_mask_and_data(mask, data);
        generate(param, tree);
        float metr = metric->get(tree);
        if (metr > max_metr)
        {
            max_metr = metr;
            data_max = data;
        }
    }
    param.load_from_mask_and_data(mask, data_max);
    return max_metr;
}
void ParameterSelector::select(TreeStructureParameters &param, SelectionType sel_type, MetricType metric_type)
{
    Metric *metric = nullptr;
    CompressionMetric default_m;
    metric = &default_m;
    if (metric_type == ImpostorSimilarity)
    {
        Texture ref = textureManager.get("reference_tree_test");
        ImpostorMetric im = ImpostorMetric(ref);
        metric = &im;
    }

    if (sel_type == BruteForce)
    {
        float m = brute_force_selection(param,metric,generate);
        logerr("bruteforce parameter selection finished with max_metric %f", m);
    } 
}