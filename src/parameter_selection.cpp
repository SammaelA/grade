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
float brute_force_selection(TreeStructureParameters &param, Metric *metric,
                            std::function<void(TreeStructureParameters &, GrovePacked &)> &generate,
                            SetSelectionProgram &set_selection_program)
{
    if (set_selection_program.selections.empty() || !metric)
        return 0;
    if (set_selection_program.schedule == AllInOne && set_selection_program.selections.size() > 1)
    {
        logerr("AllInOne schedule for multiple sets is not yet implemented");
        set_selection_program.schedule = SetbySet;
    }

    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param.get_mask_and_data(mask, data);
    data_max = data;
    GrovePacked tree;
    float max_metr = 0;

    for (auto &set : set_selection_program.selections)
    {
        int par_count = set.size();
        if (par_count == 0)
            continue;
        std::map<std::string, int> par_data_pos_by_name;
        int cur_pos = 0;
        for (auto &desk : mask)
        {
            if (desk.mask != ParameterMaskValues::CONSTANT)
            {
                par_data_pos_by_name.emplace(desk.name,cur_pos);
            }
            cur_pos += desk.var_count;
        }
        struct ParamTmp
        {
            int offset = 0;
            int count = 0;
            int data_pos = 0; //pos of modified parameter in data vector
            std::vector<float> values;
            std::string name = "unknown";
        };
        std::vector<ParamTmp> valid_params;
        int next_offset = 1;
        for (int i=0;i<par_count;i++)
        {
            if (set[i].base_values_set.empty())
                continue;
            auto it = par_data_pos_by_name.find(set[i].parameter_name);
            if (it != par_data_pos_by_name.end())
            {
                ParamTmp ptmp;
                ptmp.offset = next_offset;
                ptmp.data_pos = it->second;
                ptmp.values = set[i].base_values_set;
                ptmp.count = set[i].base_values_set.size();
                ptmp.name = it->first;
                valid_params.push_back(ptmp);
                next_offset *= ptmp.count;
            }
            else
            {
                debugl(4, "warning: parameter %s in selection set is constant or does not exist\n",
                       set[i].parameter_name);
            }
        }
        int val_count = next_offset;
        debugl(4, "Brute force parameter selection set up with %d parameters totally\n",val_count);

        const int print_iter = 100;

        for (int i = 0; i<val_count;i++)
        {
            for (auto &ptmp : valid_params)
            {
                int num = i / ptmp.offset % ptmp.count;
                data[ptmp.data_pos] = ptmp.values[num];
            }

            textureManager.set_textures_tag(1);
            GrovePacked tree = GrovePacked();
            param.load_from_mask_and_data(mask, data);
            generate(param, tree);
            float metr = metric->get(tree);
            if (metr > max_metr)
            {
                max_metr = metr;
                data_max = data;
            }
            textureManager.clear_unnamed_with_tag(1);

            if (i % print_iter == 0 || i == val_count - 1)
            {
                debugl(4,"iter %d/%d max_metric %f\n",i+1,val_count,max_metr);
            }
        }
    }

    param.load_from_mask_and_data(mask, data_max);
    return max_metr;
}
float simulated_annealing_selection(TreeStructureParameters &param, Metric *metric,
                                    std::function<void(TreeStructureParameters &, GrovePacked &)> &generate)
{
    double first_run, second_run, third_run;        //(first, second and third run) are defined for the purpose of comparing the resulting
    const double  alpha = 0.995;                         //alpha is used for the cooling schedule of the temperature            
    const double e = 2.718281828;
    const double min_metric = 0.1;

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
    double max_L = L;
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
            float pos_change_chance = CLAMP(sqrt(T/max_T), 0.025, 0.25);
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

        if (LNew > L)
        {
            data_max = data;
            L = LNew;
            max_L = LNew;
        }
        else if (LNew > min_metric && ((rand() / (double)RAND_MAX) <= pow(e, (LNew - L) / T)))
        {
            L = LNew;
        }
        else
        {
            //restore data
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
    return max_L;
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
        //Texture ref = textureManager.get("leaf1");
       //ImpostorMetric im = ImpostorMetric(ref);
       TreeSilhouette sil;
       sil.trunk_down_r = 0.04;
       sil.trunk_height = 0.6;
       sil.trunk_up_r = 0.02;
       sil.crown_center_height = 0.6;
       sil.crown_height_r = 0.3;
       sil.crown_width_r = 0.2;
       sil.crown_ellipsoid_power = 2;
       ImpostorMetric im = ImpostorMetric(sil);
        metric = &im;
    }
    if (sel_type == BruteForce)
    {
        float m = brute_force_selection(param,metric,generate);
        logerr("bruteforce parameter selection finished with max_metric %f", m);
    } 
    else if (sel_type == SimulatedAnnealing)
    {
        //float m = simulated_annealing_selection(param,metric,generate);
        SetSelectionProgram set_p;
        set_p.schedule = SelectionSchedule::SetbySet;
        SelectionSet set;
        SelectionUnit u1;
        u1.parameter_name = "base_branch_feed";
        u1.base_values_set = {50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150};
        SelectionUnit u3;
        u3.parameter_name = "base_seg_feed";
        u3.base_values_set = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150};
        set = {u1,u3};
        set_p.selections = {set};
        float m = brute_force_selection(param, metric, generate, set_p);
        logerr("simulated annealing parameter selection finished with max_metric %f", m);
    }
}