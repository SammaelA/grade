#include "parameter_selection.h"
#include "metric.h"
#include "distribution.h"
#include <iostream>
#include <functional>
#include <chrono>
#include <set>

float brute_force_selection(ParametersSet *param, Metric *metric,
                            std::function<void(ParametersSet *, GrovePacked &)> &generate)
{
    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param->get_mask_and_data(mask, data);
    data_max = data;
    GrovePacked tree;
    float max_metr = 0;

    for (int i = 0; i < 2; i++)
    {
        GrovePacked tree = GrovePacked();
        data[0] = (i + 1.0)/12;
        param->load_from_mask_and_data(mask, data);
        generate(param, tree);
        float metr = metric->get(tree);
        if (metr > max_metr)
        {
            max_metr = metr;
            data_max = data;
        }
    }
    param->load_from_mask_and_data(mask, data_max);
    return max_metr;
}
float brute_force_selection(ParametersSet *param, Metric *metric,
                            std::function<void(ParametersSet *, GrovePacked &)> &generate,
                            SetSelectionProgram &set_selection_program)
{
    std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();
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
    param->get_mask_and_data(mask, data);
    data_max = data;
    GrovePacked tree;
    float max_metr = 0;
    int global_iters = 0;
    std::vector<SelectionSet> &selections = set_selection_program.selections;
    if (set_selection_program.schedule == UnitbyUnit)
    {
        //make a new set from every unit;
        std::vector<SelectionSet> new_selections;
        for (auto &set : set_selection_program.selections)
        {
            for (auto &unit : set)
            {
                new_selections.push_back({unit});
            }
        }
        selections = new_selections;
    }
    int set_count = 0;
    bool finished = false;
    for (auto &set : selections)
    {
        if (finished)
            break;
        bool local_finished = false;
        std::vector<double> local_max;
        float local_max_metr = 0;
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
                       set[i].parameter_name.c_str());
            }
        }
        int val_count = next_offset;
        set_count++;
        debugl(4, "Brute force parameter selection set %d with %d parameters totally\n",set_count,val_count);

        const int print_iter = 10;
        std::vector<bool> visited = {};
        if (set_selection_program.schedule == SetbySetRandomized)
        {
            visited = std::vector<bool>(val_count, false);
        }
        for (int i = 0; i<val_count;i++)
        {
            int val_id = i;
            if (set_selection_program.schedule == SetbySetRandomized)
            {
                int rnd_id = urand(0, val_count);
                if (visited[rnd_id])
                {
                    bool all_visited = true;
                    int add = 1;
                    while (all_visited)
                    {
                        if (rnd_id - add >= 0 && !visited[rnd_id - add])
                        {
                            all_visited = false;
                            val_id = rnd_id - add;
                        }
                        else
                        if (rnd_id + add < val_count && !visited[rnd_id + add])
                        {
                            all_visited = false;
                            val_id = rnd_id + add;
                        }
                        add++;
                    }
                }
                else
                {
                    val_id = rnd_id;
                }
                visited[val_id] = true;
                //logerr("%d visited", val_id);
            }
            for (auto &ptmp : valid_params)
            {
                int num = val_id / ptmp.offset % ptmp.count;
                data[ptmp.data_pos] = ptmp.values[num];
            }

            textureManager.set_textures_tag(1);
            GrovePacked tree = GrovePacked();
            param->load_from_mask_and_data(mask, data);
            generate(param, tree);
            global_iters++;
            float metr = metric->get(tree);
            if (metr > local_max_metr)
            {
                local_max_metr = metr;
                local_max = data;
            }
            textureManager.clear_unnamed_with_tag(1);

            std::chrono::steady_clock::time_point t_cur = std::chrono::steady_clock::now();
            auto delta_t = t_cur - t_start;
            float dt_sec = 1e-9*delta_t.count();
            if (dt_sec > set_selection_program.exit_conditions.max_time_seconds)
            {
                finished = true;
                debugl(4, "Maximum selection time exceeded. Finishing selection\n");
            }
            else if (global_iters > set_selection_program.exit_conditions.max_iters)
            {
                finished = true;
                debugl(4, "Maximum selection iteration count exceeded. Finishing selection\n");
            }
            else if (local_max_metr > set_selection_program.exit_conditions.metric_reached)
            {
                finished = true;
                debugl(4, "Desired quality level reached. Finishing selection\n");
            }
            else if  ((float)i/val_count > set_selection_program.exit_conditions.part_of_set_covered)
            {
                local_finished = true;
                debugl(4, "Specified part of values tested. Moving to next set\n");
            }
            if ((i + 1) % print_iter == 0 || i == val_count - 1)
            {
                debugl(4,"iter %d/%d max_metric %f time spent %f\n",i+1,val_count,local_max_metr,dt_sec);
            }

            if (finished || local_finished)
                break;
        }
        if (local_max_metr > max_metr)
        {
            data_max = local_max;
            max_metr = local_max_metr;
        }
        data = data_max;
        debugl(4, "parameter set %d processed local_max = %f, global_max = %f\n",
               set_count, local_max_metr, max_metr);
    }

    param->load_from_mask_and_data(mask, data_max);
    return max_metr;
}
float simulated_annealing_selection(ParametersSet *param, Metric *metric,
                                    std::function<void(ParametersSet *, GrovePacked &)> &generate,
                                    int *iterations = nullptr,
                                    double _alpha = 0.9,
                                    double _min_metric = 0.1,
                                    double val_interval = 0.5,
                                    double pos_change_chance_min = 0.05,
                                    double pos_change_chance_max = 0.25)
{
    double first_run, second_run, third_run;        //(first, second and third run) are defined for the purpose of comparing the resulting
    const double  alpha = _alpha;                         //alpha is used for the cooling schedule of the temperature            
    const double e = 2.718281828;
    const double min_metric = _min_metric;

    std::vector<ParameterDesc> mask;
    std::vector<double> data;
    std::vector<double> data_max;
    param->get_mask_and_data(mask, data);
    GrovePacked tree;

    std::function<float(std::vector<double> &)> f = [&](std::vector<double> &d)
    {   
        textureManager.set_textures_tag(1);
        GrovePacked tree = GrovePacked();
        param->load_from_mask_and_data(mask, data);
        generate(param, tree);
        float metr = metric->get(tree);
        textureManager.clear_unnamed_with_tag(1);
        return metr;
    };

    data_max = data;
    double L = f(data);
    double max_L = L;
    double max_T = 100;
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
            //logerr("pos changed %d %f --> %f", och, prev_val[och], data[och]);
        }
        else
        {
            float pos_change_chance = CLAMP(sqrt(T/max_T), pos_change_chance_min, pos_change_chance_max);
            for (int i=0;i<data.size();i++)
            {
                if (urand()<pos_change_chance)
                {
                    pos_changed[i] = true;
                    float val_changed = val_interval*sqrt(T/max_T)*urand(-1,1);
                    prev_val[i] = data[i];
                    data[i] = CLAMP(data[i] + val_changed, 0, 1);
                    //logerr("pos changed %d %f --> %f", i, prev_val[i], data[i]);
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
    param->load_from_mask_and_data(mask, data_max);
    if (iterations)
        *iterations = k;
    return max_L;
}

float simulated_annealing_selection(ParametersSet *param, Metric *metric,
                            std::function<void(ParametersSet *, GrovePacked &)> &generate,
                            SetSelectionProgram &set_selection_program)
{
    struct SACenter
    {
        int val_id;
        float metr;
        SACenter(int id, float met)
        {
            val_id = id;
            metr = met;
        }
    };
    struct compare
    {
        bool operator()(const SACenter &j1, const SACenter &j2) const
        {
            return j1.metr > j2.metr;
        }
    };

    struct ParamTmp
    {
        int offset = 0;
        int count = 0;
        int data_pos = 0; //pos of modified parameter in data vector
        std::vector<float> values;
        std::string name = "unknown";
    };

    std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();
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
    param->get_mask_and_data(mask, data);
    data_max = data;
    GrovePacked tree;
    float max_metr = 0;
    int global_iters = 0;
    int global_centers = 0;
    std::vector<SelectionSet> &selections = set_selection_program.selections;
    if (set_selection_program.schedule == UnitbyUnit)
    {
        //make a new set from every unit;
        std::vector<SelectionSet> new_selections;
        for (auto &set : set_selection_program.selections)
        {
            for (auto &unit : set)
            {
                new_selections.push_back({unit});
            }
        }
        selections = new_selections;
    }
    int set_count = 0;
    bool finished = false;
    for (auto &set : selections)
    {
        if (finished)
            break;
        bool local_finished = false;
        std::vector<double> local_max;
        float local_max_metr = 0;
        int par_count = set.size();
        if (par_count == 0)
            continue;
        struct PdTmpData
        {
            int pos;
            float mn;
            float mx;

        };
        std::map<std::string, PdTmpData> par_data_pos_by_name;
        int cur_pos = 0;
        for (auto &desk : mask)
        {
            if (desk.mask != ParameterMaskValues::CONSTANT)
            {
                par_data_pos_by_name.emplace(desk.name,PdTmpData{cur_pos,desk.original.get_min(),desk.original.get_max()});
            }
            cur_pos += desk.var_count;
        }

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
                ptmp.data_pos = it->second.pos;
                ptmp.values = set[i].base_values_set;
                ptmp.count = set[i].base_values_set.size();
                ptmp.name = it->first;
                valid_params.push_back(ptmp);
                next_offset *= ptmp.count;
                for (float &val : ptmp.values)//normalize values
                {
                    val = (val - it->second.mn)/(it->second.mx - it->second.mn);
                }
            }
            else
            {
                debugl(4, "warning: parameter %s in selection set is constant or does not exist\n",
                       set[i].parameter_name.c_str());
            }
        }
        int val_count = next_offset;
        set_count++;
        debugl(4, "Polycentric simulated annealing parameter selection set %d with %d parameters totally\n",set_count,val_count);

        const int print_iter = 10;
        std::multiset<SACenter, compare> centers;
        std::vector<bool> visited = {};
        if (set_selection_program.schedule == SetbySetRandomized)
        {
            visited = std::vector<bool>(val_count, false);
        }

        for (int i = 0; i<val_count;i++)
        {
            int val_id = i;
            if (set_selection_program.schedule == SetbySetRandomized)
            {
                int rnd_id = urand(0, val_count);
                if (visited[rnd_id])
                {
                    bool all_visited = true;
                    int add = 1;
                    while (all_visited)
                    {
                        if (rnd_id - add >= 0 && !visited[rnd_id - add])
                        {
                            all_visited = false;
                            val_id = rnd_id - add;
                        }
                        else
                        if (rnd_id + add < val_count && !visited[rnd_id + add])
                        {
                            all_visited = false;
                            val_id = rnd_id + add;
                        }
                        add++;
                    }
                }
                else
                {
                    val_id = rnd_id;
                }
                visited[val_id] = true;
                //logerr("%d visited", val_id);
            }
            for (auto &ptmp : valid_params)
            {
                int num = val_id / ptmp.offset % ptmp.count;
                data[ptmp.data_pos] = ptmp.values[num];
            }

            textureManager.set_textures_tag(1);
            GrovePacked tree = GrovePacked();
            //param->load_from_mask_and_data(mask, data);
            generate(param, tree);
            global_iters++;
            float metr = metric->get(tree);
            centers.emplace(SACenter(val_id,metr));
            /*if (metr > local_max_metr)
            {
                local_max_metr = metr;
                local_max = data;
            }*/
            textureManager.clear_unnamed_with_tag(1);

            std::chrono::steady_clock::time_point t_cur = std::chrono::steady_clock::now();
            auto delta_t = t_cur - t_start;
            float dt_sec = 1e-9*delta_t.count();
                        
            if ((i + 1) % print_iter == 0 || i == val_count - 1)
            {
                debugl(4,"preparing center %d/%d time spent %f\n",i+1,val_count,dt_sec);
            }
            if (dt_sec > set_selection_program.exit_conditions.max_time_seconds)
            {
                finished = true;
                debugl(4, "Maximum selection time exceeded. Finishing selection\n");
            }
            else if (global_iters > set_selection_program.exit_conditions.max_iters)
            {
                finished = true;
                debugl(4, "Maximum selection iteration count exceeded. Finishing selection\n");
            }
            else if (local_max_metr > set_selection_program.exit_conditions.metric_reached)
            {
                finished = true;
                debugl(4, "Desired quality level reached. Finishing selection\n");
            }
            if (finished || local_finished)
                break;
        }
        
        int centers_count = centers.size();
        int centers_processed = 0;
        for (const SACenter &center : centers)
        {
            if (local_finished || finished)
                break;
            //param->load_from_mask_and_data(mask, data);
            int sa_iter = 0;
            float metr = simulated_annealing_selection(param, metric, generate, &sa_iter, 0.6, 0.1, 0.15);
            global_iters += sa_iter;
            if (metr > local_max_metr)
            {
                local_max.clear();
                std::vector<ParameterDesc> tmp_mask;
                param->get_mask_and_data(tmp_mask, local_max);
                local_max_metr = metr;

            }       

            std::chrono::steady_clock::time_point t_cur = std::chrono::steady_clock::now();
            auto delta_t = t_cur - t_start;
            float dt_sec = 1e-9*delta_t.count();
            debugl(4, "processed center %d, metr %f --> %f time spent %f\n",center.val_id, center.metr, metr, dt_sec);
            if (dt_sec > set_selection_program.exit_conditions.max_time_seconds)
            {
                finished = true;
                debugl(4, "Maximum selection time exceeded. Finishing selection\n");
            }
            else if (global_iters > set_selection_program.exit_conditions.max_iters)
            {
                finished = true;
                debugl(4, "Maximum selection iteration count exceeded. Finishing selection\n");
            }
            else if (local_max_metr > set_selection_program.exit_conditions.metric_reached)
            {
                finished = true;
                debugl(4, "Desired quality level reached. Finishing selection\n");
            }
            if ((float)centers_processed/centers_count > set_selection_program.exit_conditions.part_of_set_covered)
            {
                local_finished = true;
                debugl(4, "Specified part of values tested. Moving to next set\n");
            }
            centers_processed++;
        }

        if (local_max_metr > max_metr)
        {
            data_max = local_max;
            max_metr = local_max_metr;
        }
        data = data_max;
        global_centers += centers_processed;
        debugl(4, "parameter set %d processed local_max = %f, global_max = %f\n",
               set_count, local_max_metr, max_metr);
    }

    std::chrono::steady_clock::time_point t_cur = std::chrono::steady_clock::now();
    auto delta_t = t_cur - t_start;
    float dt_sec = 1e-9*delta_t.count();

    debugl(4, "\nPolycentric simulated allealing finished\n");
    debugl(4, "Max metric: %f\n",max_metr);
    debugl(4, "Parameter sets processed: %d\n",(int)selections.size());
    debugl(4, "Centers tested: %d\n",(int)global_centers);
    debugl(4, "Trees generated: %d\n",(int)global_iters);
    debugl(4, "Time spent: %d h %d m %d s\n",(int)dt_sec/3600, ((int)dt_sec)/60 % 60, ((int)dt_sec)%60);


    //param->load_from_mask_and_data(mask, data_max);
    return max_metr;
}
void ParameterSelector::select(ParametersSet *param, std::string selection_program_name)
{
    BlkManager man;
    Block blk;
    man.load_block_from_file("parameter_selection.blk",blk);
    SetSelectionProgram set_p;
    load_from_blk(set_p,selection_program_name,blk);
    Metric *metric = nullptr;
    DummyMetric default_m;
    metric = &default_m;
    if (set_p.metr_type == CompressionRatio)
    {
        CompressionMetric cm;
        metric = &cm;
    }
    if (set_p.metr_type == ImpostorSimilarity)
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
    if (set_p.sel_type == BruteForce)
    {
        float m = brute_force_selection(param,metric,generate);
        logerr("bruteforce parameter selection finished with max_metric %f", m);
    } 
    else if (set_p.sel_type == PolycentricSimulatedAnnealing)
    {
        float m = simulated_annealing_selection(param, metric, generate, set_p);
        logerr("simulated annealing parameter selection finished with max_metric %f", m);
    }
}

void ParameterSelector::save_to_blk(SetSelectionProgram &prog, std::string name, Block &blk)
{
    Block *bl = new Block();
    bl->add_string("schedule",ToString(prog.schedule));
    bl->add_string("selection_type",ToString(prog.sel_type));
    bl->add_string("metric_type",ToString(prog.metr_type));
    for (SelectionSet &set : prog.selections)
    {
        Block *set_bl = new Block();
        for (SelectionUnit &unit : set)
        {
            set_bl->add_arr(unit.parameter_name,unit.base_values_set);
        }
        bl->add_block("selection_set",set_bl);
        
    }
    Block *ec_bl = new Block();
    ec_bl->add_double("max_iters",prog.exit_conditions.max_iters);
    ec_bl->add_double("max_time_seconds",prog.exit_conditions.max_time_seconds);
    ec_bl->add_double("metric_reached",prog.exit_conditions.metric_reached);
    ec_bl->add_double("part_of_set_covered",prog.exit_conditions.part_of_set_covered);
    bl->add_block("exit_conditions",ec_bl);
    blk.add_block(name,bl);
}
void ParameterSelector::load_from_blk(SetSelectionProgram &prog, std::string name, Block &blk)
{
    prog.selections.clear();
    prog.schedule = SelectionSchedule::UnitbyUnit;

    Block *bl = blk.get_block(name);
    std::string schedule_name = bl ? bl->get_string("schedule","UnitbyUnit") : "";
    if (!bl || schedule_name.empty())
    {
        logerr("cannot read SetSelectionProgram %s from block",name);
        return;
    }
    std::string sel_type = bl->get_string("selection_type","BruteForce");
    std::string metr_type = bl->get_string("metric_type","Dummy");
    for (int i=0;i<=(int)SelectionSchedule::UnitbyUnit;i++)
    {
        if (ToString((SelectionSchedule)i) == schedule_name)
        {
            prog.schedule = (SelectionSchedule)i;
            break;
        }
    }
    for (int i=0;i<=(int)SelectionType::PolycentricSimulatedAnnealing;i++)
    {
        if (ToString((SelectionType)i) == sel_type)
        {
            prog.sel_type = (SelectionType)i;
            break;
        }
    }
    for (int i=0;i<=(int)MetricType::ImpostorSimilarity;i++)
    {
        if (ToString((MetricType)i) == metr_type)
        {
            prog.metr_type = (MetricType)i;
            break;
        }
    }
    Block *ec_bl = bl->get_block("exit_conditions");
    if (ec_bl)
    {
        prog.exit_conditions.max_iters = ec_bl->get_double("max_iters",prog.exit_conditions.max_iters);
        prog.exit_conditions.max_time_seconds = ec_bl->get_double("max_time_seconds",prog.exit_conditions.max_time_seconds);
        prog.exit_conditions.metric_reached = ec_bl->get_double("metric_reached",prog.exit_conditions.metric_reached);
        prog.exit_conditions.part_of_set_covered = ec_bl->get_double("part_of_set_covered",prog.exit_conditions.part_of_set_covered);
    }
    int id = bl->get_id("selection_set");
    while (id >= 0)
    {
        Block *sel_set = bl->get_block(id);
        if (!sel_set)
        {
            id = -1;
        }
        else
        {
            int params = sel_set->size();

            if (params <= 0)
            {
                logerr("cannot read SetSelectionProgram %s from block. selection set is empty",name);
                id = -1;
            }
            else
            {
                prog.selections.emplace_back();
                SelectionSet &set = prog.selections.back();
                for (int i=0;i<params;i++)
                {
                    set.emplace_back();
                    SelectionUnit &unit = set.back();
                    unit.parameter_name = sel_set->get_name(i);
                    sel_set->get_arr(i,unit.base_values_set);
                    if (unit.base_values_set.empty())
                    {
                        logerr("cannot read SetSelectionProgram %s from block. selection unit %s is empty",name,
                                unit.parameter_name );
                        id = -1;
                    }
                }
            }
        }
        if (id >= 0)
            id = sel_set->get_next_id("selection_set",id+1);
    }
    if (prog.selections.empty())
    {
        logerr("cannot read SetSelectionProgram %s from block. Empty selection program",name);
    }
}