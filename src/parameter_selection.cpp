#include "parameter_selection.h"
#include "metric.h"
#include <iostream>
#include <functional>

float brute_force_selection(TreeStructureParameters &param, Metric *metric,
                            std::function<void(TreeStructureParameters &, GrovePacked &)> &generate)
{
    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param.get_mask_and_data(mask, data);
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
    time_t systime;                                // solutions of the three runs will be chosen as the final solution
    time(&systime);
    srand((unsigned int)systime);
    double  alpha = 0.75;                         //alpha is used for the cooling schedule of the temperature            
    const double e = 2.718281828;


    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param.get_mask_and_data(mask, data);
    GrovePacked tree;

    std::function<float(std::vector<double> &)> f = [&](std::vector<double> &d)
    {
        GrovePacked tree = GrovePacked();
        param.load_from_mask_and_data(mask, data);
        generate(param, tree);
        float metr = metric->get(tree);
        //textureManager.clear_unnamed();
        return metr;
    };
    //std::cout << "Initial State = " << x << "\t, and F(x)= " << f(x) << std::endl;

    double L = f(data);
    int k = 0;
    for (double T = 80; T > 0.08; T *= alpha) //T = T * alpha which used as a cooling schedule 
    {
        k++;
        logerr("T = %f", T);
        for (int i = 0; i<2; i++) //This loop is for the process of iteration (or searching for new states)
        {
            int pos_changed = rand() % data.size();
            float val_changed = (rand() / (double)RAND_MAX);
            float prev_val = data[pos_changed];
            data[pos_changed] = CLAMP(data[pos_changed] + val_changed, 0, 1);
            //double xNew = x + ((rand() / (double)RAND_MAX) * 2 - 1);
            double LNew = f(data);

            if (LNew > L || (rand() / (double)RAND_MAX) <= pow(e, (LNew - L) / T))
            {
                L = LNew;
                if (LNew > L)
                    data_max = data;
            }
            else
            {
                data[pos_changed] = prev_val;
            }
        }
        if (k > 25)
            break;
    }

    param.load_from_mask_and_data(mask, data_max);
    return L;
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
    else if (sel_type == SimulatedAnnealing)
    {
        float m = simulated_annealing_selection(param,metric,generate);
        logerr("simulated annealing parameter selection finished with max_metric %f", m);
    }
}