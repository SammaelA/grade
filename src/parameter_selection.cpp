#include "parameter_selection.h"
#include "metric.h"
#include "distribution.h"
#include <iostream>
#include <functional>

float brute_force_selection(TreeStructureParameters &param, Metric *metric,
                            std::function<void(TreeStructureParameters &, GrovePacked &)> &generate)
{
    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param.get_mask_and_data(mask, data);
    data_max = data;
    GrovePacked tree;
    float max_metr = 0;

    for (int i = 0; i < 2; i++)
    {
        GrovePacked tree = GrovePacked();
        data[0] = (i + 1.0)/12;
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

float simulated_annealing_selection(TreeStructureParameters &param, Metric *metric,
                                    std::function<void(TreeStructureParameters &, GrovePacked &)> &generate)
{
    double first_run, second_run, third_run;        //(first, second and third run) are defined for the purpose of comparing the resulting
    double  alpha = 0.975;                         //alpha is used for the cooling schedule of the temperature            
    const double e = 2.718281828;


    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param.get_mask_and_data(mask, data);
    GrovePacked tree;

    std::function<float(std::vector<double> &)> f = [&](std::vector<double> &d)
    {   
        textureManager.set_textures_tag(1);
        GrovePacked tree = GrovePacked();
        param.load_from_mask_and_data(mask, data);
        generate(param, tree);
        float metr = metric->get(tree);
        textureManager.clear_unnamed_with_tag(1);
        return metr;
    };

    data_max = data;
    double L = f(data);
    double max_T = 80;
    int k = 0;
    bool *pos_changed = safe_new<bool>(data.size(),"pos_changed_arr");
    float *prev_val = safe_new<float>(data.size(),"prev_val_arr");
    bool single_param = false;
    for (double T = max_T; T > 1; T *= alpha)
    {
        k++;
        logerr("T = %f", T);

        if (single_param)
        {
            int och = urand(0, data.size());
            pos_changed[och] = true;
            float val_changed = sqrt(T/max_T)*urand(-1,1);
            prev_val[och] = data[och];
            data[och] = CLAMP(data[och] + val_changed, 0, 1);
            logerr("pos changed %d %f --> %f", och, prev_val[och], data[och]);
        }
        else
        {
            float pos_change_chance = CLAMP(sqrt(T/max_T), 0.025, 0.025);
            for (int i=0;i<data.size();i++)
            {
                if (urand()<pos_change_chance)
                {
                    pos_changed[i] = true;
                    float val_changed = sqrt(T/max_T)*urand(-1,1);
                    prev_val[i] = data[i];
                    data[i] = CLAMP(data[i] + val_changed, 0, 1);
                    logerr("pos changed %d %f --> %f", i, prev_val[i], data[i]);
                }
                else
                {
                    pos_changed[i] = false;
                }
            }
        }

        double LNew = f(data);

        if (LNew > L || (rand() / (double)RAND_MAX) <= pow(e, (LNew - L) / T))
        {
            if (LNew > L)
                data_max = data;
            L = LNew;
        }
        else
        {
            for (int i=0;i<data.size();i++)
            {
                if (pos_changed[i])
                    data[i] = prev_val[i];
            }
        }
    }
    safe_delete<bool>(pos_changed,"pos_changed_arr");
    safe_delete<float>(prev_val,"prev_val_arr");
    param.load_from_mask_and_data(mask, data_max);
    logerr("base r %s",param.base_r.to_string().c_str());
    logerr("r_split save power %s",param.r_split_save_pow.to_string().c_str());
    return L;
}
void ParameterSelector::select(TreeStructureParameters &param, SelectionType sel_type, MetricType metric_type)
{
    Metric *metric = nullptr;
    DummyMetric default_m;
    metric = &default_m;
    if (metric_type == CompressionRatio)
    {
        CompressionMetric cm;
        metric = &cm;
    }
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
    else if (sel_type == SimulatedAnnealing)
    {
        float m = simulated_annealing_selection(param,metric,generate);
        logerr("simulated annealing parameter selection finished with max_metric %f", m);
    }
}