#include "genetic_algorithm.h"
#include "common_utils/distribution.h"
#include <algorithm>
#include <atomic>
#include "GA_utils.h"
#include <set>

std::atomic<int> next_id(0);
int di_cnt[256];
long double di_res[256];
long double di_dev[256];
void GeneticAlgorithm::perform(std::vector<float> &param_list, my_opt::MetaParameters *params, my_opt::ExitConditions exit_conditions,
                               const my_opt::OptFunction &opt_f,
                               std::vector<std::pair<float, std::vector<float>>> &best_results,
                               std::vector<std::vector<float>> &initial_types)
{
    t_start = std::chrono::steady_clock::now();

    opt_function = opt_f;
    MetaParameters *ga_p = dynamic_cast<MetaParameters*>(params);
    if (ga_p)
        metaParams = *ga_p;
    exitConditions = exit_conditions;
    free_parameters_cnt = param_list.size();
    all_parameters_cnt = param_list.size();

    metaParams.heaven_size = MAX(metaParams.heaven_size, metaParams.best_genoms_count + 1);
    metaParams.min_population_size = MAX(MAX(metaParams.min_population_size,metaParams.elite), 
                                     MAX(metaParams.heaven_size, metaParams.best_genoms_count));
    metaParams.max_population_size = MAX(metaParams.max_population_size, metaParams.min_population_size);
    next_id.store(0);
    save_load_function_stat(true);

    
    if (metaParams.type == GA_Type::ISLANDS_GA)
        islands_GA(initial_types);
    else
        tree_GA(initial_types);
    prepare_best_params(best_results);
    save_load_function_stat(false);

    std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    debug("GA result: best metric %f took %.2f seconds and %d tries to find\n", best_metric_ever, time / 1000, func_called);
}
int base_size_by_depth(int base_size, int cur_depth)
{
    return base_size;
}
void GeneticAlgorithm::tree_GA(std::vector<std::vector<float>> &initial_types)
{
    PopulationBackup result;
    int base_size = metaParams.tree_GA_iters;
    int base_width = 2;
    int tries = exitConditions.function_calculated;
    int cur_tries = 0;
    int cur_depth = 0;
    while (cur_tries < tries)
    {
        cur_depth++;
        cur_tries = base_width*cur_tries + (base_size_by_depth(base_size, cur_depth) - 1)*metaParams.max_population_size*metaParams.weaks_to_kill + 
                    metaParams.initial_population_size;
    }
    int sz = base_size * tries/ cur_tries;
    logerr("starting tree GA depth %d, size %d, tries %d", cur_depth, sz, sz*cur_tries/base_size);
    TreeGA_stat = std::vector<glm::vec2>(cur_depth+1, glm::vec2(0,0));
    tree_GA_internal(cur_depth, sz, base_width, initial_types, result);
    logerr("tree GA stat:");
    for (auto &s : TreeGA_stat)
    {
        logerr("%.3f %d", s[1] > 0 ? s[0]/s[1] : -1, (int)s[1]);
    }
    population = result.pop;
}

float MM_func(glm::vec4 model, float x)
{
    float x_p = model.z < 0  ? 1/(-model.z) : model.z;
    x = pow(x, x_p);
    return model.x*x/(model.y + x);
}

void MM_regression_grid(glm::vec4 model_start, glm::vec4 model_end, const std::vector<float> &values, 
                        float &least_dist, glm::vec4 &best_model, int level)
{
    int steps = 20;
    glm::vec4 step = (model_end - model_start)/(float)steps;
    for (int i=0;i<steps;i++)
    {
        for (int j=0;j<steps;j++)
        {
            for (int k=0;k<steps;k++)
            {
                glm::vec4 model = model_start + glm::vec4(step.x*(float)(i + 0.5), step.y*(float)(j+0.5),step.z*(float)(k+0.5),0);
                float sq_dist = 0;
                for (int x=0;x<values.size();x++)
                {
                    sq_dist += SQR(values[x] - MM_func(model, x));
                }
                //logerr("model %f %f sq_dist %f",model.x, model.y, sq_dist);
                if (sq_dist < least_dist)
                {
                    least_dist = sq_dist;
                    best_model = model;
                }
            }
        }
    }
    if (level > 1)
    {
        model_start = best_model - 0.5f*step;
        model_end = best_model + 0.5f*step;
        //logerr("%f %f -- %f %f", model_start.x, model_start.y, model_end.x, model_end.y);
        MM_regression_grid(model_start, model_end, values, least_dist, best_model, level - 1);
    }
}

glm::vec4 MM_regression(const std::vector<float> &values)
{
    int steps = 100;
    float least_dist = 1e9;
    glm::vec4 best_model = glm::vec4(0,0,0,0);
    MM_regression_grid(glm::vec4(0,0,-5,0), glm::vec4(1,50,5,0), values, least_dist, best_model, 3);

    return best_model;
}

void GeneticAlgorithm::tree_GA_internal(int depth, int iters, int width, std::vector<std::vector<float>> &initial_types, PopulationBackup &result)
{
    bool regression_test = false;
    metaParams.n_islands = 1;
    if (depth == 1)
    {
        initialize_population(initial_types);
    }
    else
    {
        std::vector<PopulationBackup> backups = std::vector<PopulationBackup>(width, PopulationBackup());
        int best_pos = 0;
        std::vector<Creature> all_pop;
        for (int i=0;i<width;i++)
        {
            tree_GA_internal(depth - 1, iters, width, initial_types, backups[i]);
            all_pop.insert(all_pop.end(), backups[i].pop.begin(), backups[i].pop.end());
            if (backups[i].predicted_best_value > backups[best_pos].predicted_best_value)
            {
                best_pos = i;
            }
        }
        std::sort(all_pop.begin(), all_pop.end(), 
              [&](const Creature& a, const Creature& b) -> bool{return a.alive*a.metric > b.alive*b.metric;});
        //logerr("%d %f %f", all_pop.size(), all_pop.front().metric, all_pop.back().metric);
        population = backups[best_pos].pop;
        for (int i=0;i<population.size();i++)
        {
            population[i] = all_pop[i];
        }

    }
    calculate_metric(1, false);
    recalculate_fitness();
    pick_best_to_heaven();
    //debug("iteration 0 Pop: %d Best: %.4f\n", current_population_size, best_metric_current);
    debug("[%d] start %.3f\n", depth, best_metric_current);

    iteration_n = 0;
    int cur_iters = base_size_by_depth(iters, depth);
    std::vector<float> values;
    std::vector<float> predictions;
    values.push_back(0);
    while (!should_exit() && iteration_n < cur_iters)
    {
        int c_id = next_id / 1000 * 1000;
        next_id.store(1000 + c_id);
        kill_old();
        kill_weak(current_population_size - metaParams.max_population_size * (1 - metaParams.weaks_to_kill));
        if (current_population_size < 5)
        {
            logerr("population extincted");
            return;
        }
        else if (current_population_size < metaParams.min_population_size)
        {
            logerr("population %d is about to extinct", current_population_size);
        }
        int space_left = metaParams.max_population_size - current_population_size;
        std::vector<std::pair<int, int>> pairs;
        find_pairs(space_left, pairs);
        make_new_generation(pairs);

        calculate_metric(1, true);
        recalculate_fitness();
        pick_best_to_heaven();
        iteration_n++;
        if (metaParams.evolution_stat)
        {
            float max_val = 0;
            int max_pos = 0;
            for (int i = 0; i < heaven.size(); i++)
            {
                if (heaven[i].metric > max_val)
                {
                    max_pos = i;
                    max_val = heaven[i].metric;
                }
            }
            crio_camera.push_back(heaven[max_pos]);
        }
        //debug("iteration %d Pop: %d Best: %.4f\n", iteration_n, current_population_size, best_metric_current);
        values.push_back(best_metric_current);
        if (regression_test)
        {
            glm::vec4 model = MM_regression(values);
            logerr("model %f %f %f %f", model.x, model.y, model.z);
            predictions.push_back(MM_func(model, cur_iters - 1));
        }
    }
    glm::vec4 model = MM_regression(values);

    result.pop = population;
    result.best_value = 0;
    result.model = model;
    int predicted_iter = regression_test ? cur_iters : 3*cur_iters;
    result.predicted_best_value = MM_func(model, predicted_iter);
    for (auto &p : result.pop)
    {
        result.best_value = MAX(result.best_value, p.metric);
    }
    TreeGA_stat[depth].x += result.best_value;
    TreeGA_stat[depth].y++;
    debug("[%d] end %.3f\n", depth, result.best_value);
    if (regression_test)
    {
        logerr("predicted %f",result.predicted_best_value);
        for (int i=0;i<predictions.size();i++)
        {
            logerr("%d pred %f real %f quality %f",i,predictions[i], result.best_value, 1 - abs(predictions[i]- result.best_value)/result.best_value);
        }
    }
}

void GeneticAlgorithm::islands_GA(std::vector<std::vector<float>> &initial_types)
{
    for (int i=0;i<metaParams.n_islands;i++)
    {
        sub_population_infos.emplace_back();
        sub_population_infos[i].best_result_iter = 0;
        sub_population_infos[i].no_progress_time = -1;
    }
    initialize_population(initial_types);
    calculate_metric(1, false);
    recalculate_fitness();
    pick_best_to_heaven();
    debug("iteration 0 Pop: %d Best: %.4f\n", current_population_size, best_metric_ever);
    while (!should_exit())
    {
        next_id.store(1000*(iteration_n + 1));
        kill_old();
        kill_weak(current_population_size - metaParams.max_population_size*(1-metaParams.weaks_to_kill));
        if (current_population_size < 5)
        {
            logerr("population extincted");
            return;
        }
        else if (current_population_size < metaParams.min_population_size)
        {
            logerr("population %d is about to extinct", current_population_size);
        }
        int space_left = metaParams.max_population_size - current_population_size;
        std::vector<std::pair<int, int>> pairs;
        find_pairs(space_left, pairs);
        make_new_generation(pairs);
        
        calculate_metric(1, true);
        recalculate_fitness();
        pick_best_to_heaven();
        if (iteration_n % metaParams.migration_interval == 0 && iteration_n)
            migration();
        for (int i=0;i<metaParams.n_islands;i++)
        {
            sub_population_infos[i].best_results_history.push_back(0);
        }
        for (int i=0;i<population.size();i++)
        {
            if (population[i].alive)
            {
                population[i].age++;
                auto &sp_info = sub_population_infos[population[i].sub_population_n];
                sp_info.best_results_history.back() = MAX(sp_info.best_results_history.back(), population[i].metric);
            }
        }
        for (int i=0;i<metaParams.n_islands;i++)
        {
            if (sub_population_infos[i].best_results_history.back() >= 
                sub_population_infos[i].best_results_history[sub_population_infos[i].best_result_iter] + SubPopulationInfo::progress_thr)
            {
                sub_population_infos[i].best_result_iter = sub_population_infos[i].best_results_history.size() - 1;
                sub_population_infos[i].no_progress_time = 0;
            }
            else
            {
                sub_population_infos[i].no_progress_time++;
            }
            //logerr("sp %d best %f no_progress %d",i,sub_population_infos[i].best_results_history[sub_population_infos[i].best_result_iter],sub_population_infos[i].no_progress_time);
        }
        if (iteration_n % 10 == 0)
        {
            for (int i=0;i<metaParams.n_islands;i++)
            {
                auto &info = sub_population_infos[i];
                if (info.no_progress_time < 10)
                {
                    //create backup 
                    info.backups.emplace_back();
                    auto &backup = info.backups.back();
                    backup.backup_uses = 0;
                    backup.best_value = info.best_results_history[info.best_result_iter];
                    for (int j=0;j<population.size();j++)
                    {
                        if (population[j].alive && population[j].sub_population_n == i)
                        {
                          backup.pop.push_back(population[j]);  
                        }
                    }
                    logerr("subpop %d backup %d saved %d", i, info.backups.size(), backup.pop.size());
                }   
                else
                {
                    //restore population
                    int stall_limit = MAX(10, 50*(info.best_results_history[info.best_result_iter]));
                    while (!info.backups.empty() && info.backups.back().backup_uses >= 2)
                    {
                        info.backups.pop_back();
                    }
                    if (false && !info.backups.empty() && info.no_progress_time < stall_limit)
                    {
                        int cnt = 0;
                        for (int j=0;j<population.size();j++)
                        {
                            if (population[j].alive && population[j].sub_population_n == i)
                            {
                                if (cnt >= info.backups.back().pop.size())
                                    break;
                                population[j] = info.backups.back().pop[cnt];
                                cnt++;
                            }
                        }
                        info.backups.back().backup_uses++;
                        logerr("subpop %d backup %d restored %d", i, info.backups.size(), cnt);
                    }
                    else if (info.no_progress_time >= stall_limit)
                    {
                        //create new random population here
                        info = SubPopulationInfo();
                        int cnt = 0;
                        for (int j=0;j<population.size();j++)
                        {
                            if (population[j].alive && population[j].sub_population_n == i)
                            {
                                population[j] = Creature();
                                population[j].sub_population_n = i;
                                population[j].main_genome = random_genes();
                                for (int jj=1;jj<metaParams.n_ploid_genes;jj++)
                                    population[j].other_genomes.push_back(random_genes());
                                population[j].alive = true;
                                population[j].max_age = metaParams.max_age;
                                population[j].id = next_id.fetch_add(1);
                            }
                        }
                        logerr("subpop %d recreated %d", i, cnt);
                        calculate_metric(1, true);
                        recalculate_fitness();
                        pick_best_to_heaven();
                    }
                }
            }
        }
        iteration_n++;
        if (metaParams.evolution_stat)
        {
            float max_val = 0;
            int max_pos = 0;
            for (int i=0;i<heaven.size();i++)
            {
                if (heaven[i].metric > max_val)
                {
                    max_pos = i;
                    max_val = heaven[i].metric;
                }
            }
            crio_camera.push_back(heaven[max_pos]);
        }
        debug("iteration %d Pop: %d Best: %.4f\n", iteration_n, current_population_size, best_metric_ever);
    }
}

bool GeneticAlgorithm::should_exit()
{
    iter_start = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(iter_start - t_start).count();
    time_spent_sec = 0.001*time;
    return (time_spent_sec > exitConditions.time_elapsed_seconds || iteration_n >= exitConditions.generations
    || func_called >= exitConditions.function_calculated || best_metric_ever >= exitConditions.function_reached);
}

void add_all_ids(std::set<int> &ids, int id, std::vector<glm::ivec3> &all_births)
{
    ids.emplace(id);
    for (auto &v : all_births)
    {
        if (v.z == id)
        {
            add_all_ids(ids, v.x, all_births);
            add_all_ids(ids, v.y, all_births);
        }
    }
}

void GeneticAlgorithm::prepare_best_params(std::vector<std::pair<float, std::vector<float>>> &best_results)
{
    for (int i=0;i<population.size();i++)
    {
        if (population[i].alive)
            kill_creature(i);
    }
    calculate_metric(metaParams.heaven_recalc_n, false);
    std::sort(heaven.begin(), heaven.end(), 
              [&](const Creature& a, const Creature& b) -> bool{return a.metric > b.metric;});
    debug("heaven popultion %d best metric %.4f\n", heaven.size(), heaven[0].metric);
    debug("heaven [ ");
    for (auto &creature : heaven)
    {
        debug("%.4f ", creature.metric);
    }
    debug("]\n");
    for (int i=0;i<MIN(metaParams.best_genoms_count, heaven.size());i++)
    {
        best_results.emplace_back();
        best_results.back().first = heaven[i].metric;
        best_results.back().second = heaven[i].main_genome;
        //best_results.back().second.print();
    }
    for (auto &creature : crio_camera)
    {
        best_results.emplace_back();
        best_results.back().first = creature.metric;
        best_results.back().second = creature.main_genome;
        //best_results.back().second.print();
    }

    if (metaParams.debug_graph)
    {
        DebugGraph g_test;
        g_test.add_node(DebugGraph::Node(glm::vec2(0,0), glm::vec3(1,0,0),0.1));
        g_test.add_node(DebugGraph::Node(glm::vec2(-1,1), glm::vec3(0,1,0),0.15));
        g_test.add_node(DebugGraph::Node(glm::vec2(1,1), glm::vec3(0,0,1),0.2));
        g_test.add_edge(DebugGraph::Edge(0,1,0.02));
        g_test.add_edge(DebugGraph::Edge(1,2,0.03));
        g_test.add_edge(DebugGraph::Edge(0,2,0.04));
        g_test.save_as_image("graph_test", 300, 300);

        DebugGraph stat_g;
        std::map<int, int> pos_by_id;
        int cnt = 0;
        for (auto &p : stat.all_results)
        {
            int iter = p.first/1000;
            int n = p.first % 1000;
            stat_g.add_node(DebugGraph::Node(glm::vec2(iter, 0.1*n), glm::vec3(1 - p.second, p.second, 0), 0.045));
            pos_by_id.emplace(p.first, cnt);
            cnt++;
        }

        int n = 0;
        for (auto &c : heaven)
        {
            DebugGraph stat_gn = stat_g;
            std::set<int> ids;
            add_all_ids(ids, c.id, stat.all_births);
            glm::vec3 color = glm::vec3(n/4%2, n/2%2, n%2);
            if (n % 8 == 0)
            {
                color = glm::vec3(0.5, 0.5, 0.5);
            }
            for (auto &id : ids)
            {
                for (auto &v : stat.all_births)
                {
                    if (v.z == id)
                    {
                        auto itA = pos_by_id.find(v.x);
                        auto itB = pos_by_id.find(v.y);
                        auto itC = pos_by_id.find(v.z);
                        if (itA != pos_by_id.end() && itB != pos_by_id.end() && itC != pos_by_id.end())
                        {
                            stat_gn.add_edge(DebugGraph::Edge(itA->second, itC->second, 0.005, color));
                            stat_gn.add_edge(DebugGraph::Edge(itB->second, itC->second, 0.005, color));
                        }                        
                    }
                }
            }
            n++;
            stat_gn.save_as_image("graph_gen_" + std::to_string(n), 4000, 4000);
        }

        
        for (auto &v : stat.all_births)
        {
            logerr("aaa %d %d %d",v.x, v.y, v.z);
            auto itA = pos_by_id.find(v.x);
            auto itB = pos_by_id.find(v.y);
            auto itC = pos_by_id.find(v.z);
            if (itA != pos_by_id.end() && itB != pos_by_id.end() && itC != pos_by_id.end())
            {
                stat_g.add_edge(DebugGraph::Edge(itA->second, itC->second, 0.005));
                stat_g.add_edge(DebugGraph::Edge(itB->second, itC->second, 0.005));
            }

        }
        
        stat_g.save_as_image("graph_gen", 5000, 5000);
    }
}

GeneticAlgorithm::Genome GeneticAlgorithm::random_genes()
{
    Genome g = Genome(all_parameters_cnt, 0);
    Genome best_g = g;
    float best_mark = -1000;
    int tries = 10;
    for (int t =0;t<tries;t++)
    {
        for(int i=0;i<all_parameters_cnt;i++)
            g[i] = urand();
        float mark = get_mark(g);
        if (mark > best_mark)
        {
            best_g = g;
            best_mark = mark;
        }
    }
    return best_g;
}

void GeneticAlgorithm::initialize_population(std::vector<std::vector<float>> &initial_types)
{
    population = std::vector<Creature>(MAX(metaParams.max_population_size, metaParams.initial_population_size), Creature());
    new_population = std::vector<Creature>(MAX(metaParams.max_population_size, metaParams.initial_population_size), Creature());
    
    int max_iters = 10;
    int max_cnt = 0.5*metaParams.initial_population_size;
    int i0=0;
    int iter = 0;
    
    while (i0 < max_cnt && iter < max_iters)
    {
        for (auto &t : initial_types)
        {
            population[i0].main_genome = t;
            for (int j=1;j<metaParams.n_ploid_genes;j++)
                population[i0].other_genomes.push_back(population[i0].main_genome);
            population[i0].alive = true;
            population[i0].max_age = metaParams.max_age;
            population[i0].id = next_id.fetch_add(1);
            population[i0].sub_population_n = metaParams.n_islands*i0/metaParams.initial_population_size;
            i0++;
            if (i0 >= max_cnt)
                break;
        }
        iter++;
    }

    for (int i=i0;i<metaParams.initial_population_size;i++)
    {
        population[i].main_genome = random_genes();
        for (int j=1;j<metaParams.n_ploid_genes;j++)
            population[i].other_genomes.push_back(random_genes());
        population[i].alive = true;
        population[i].max_age = metaParams.max_age;
        population[i].id = next_id.fetch_add(1);
        population[i].sub_population_n = metaParams.n_islands*i/metaParams.initial_population_size;
    }
    current_population_size = metaParams.initial_population_size;
}

float GeneticAlgorithm::get_mark(Genome &G)
{
    float mark = 0;
    for (int i=0;i<G.size();i++)
    {
        int bucket = CLAMP(G[i]*function_stat.Q_NUM,0,function_stat.Q_NUM-1);
        mark += function_stat.marks[i][bucket];
    }

    return mark;
}

void GeneticAlgorithm::mutation(Genome &G, float mutation_power, int mutation_genes_count, int *single_mutation_pos)
{
    std::vector<float> mutation_weights = std::vector<float>(G.size(), 1);
    float w_sum = G.size();
    if (!single_mutation_pos && metaParams.use_genes_importance)
    {
        w_sum = 0;
        for (int i=0;i<G.size();i++)
        {
            mutation_weights[i] = MAX(0.05, di_res[i]/MAX(1, di_cnt[i]));
            w_sum += mutation_weights[i];
        }
    }
    Genome GBase = G;
    Genome GBest = G;
    float best_mark = -100;
    int g_tries = function_stat.tries > 10000 ? 100 : 1;
    g_tries = 1;
    for (int g_try = 0;g_try < g_tries;g_try++)
    {
        for (int gene = 0; gene < mutation_genes_count; gene++)
        {
            bool found = false;
            int tries = 0;
            while (!found && tries < 100)
            {
                float rnd = urand(0, w_sum);  
                int pos = 0;
                for (int i=0;i<G.size();i++)
                {
                    if (rnd < mutation_weights[i])
                    {
                        pos = i;
                        found = true;
                        break;
                    }   
                    else
                    {
                        rnd -= mutation_weights[i];
                    }
                }
                int g_pos = pos;
                if (single_mutation_pos)
                {
                    G[g_pos] = G[g_pos] + mutation_power*urand(-1,1);
                    G[g_pos] = CLAMP(G[g_pos],0,1);
                }
                else 
                    G[g_pos] = urand(0,1);
                tries++;
                if (single_mutation_pos && found == true && gene == mutation_genes_count-1)
                {
                    //logerr("mutation %d", g_pos);
                    *single_mutation_pos = g_pos;
                }
            }
        }
        float mark = get_mark(G);
        if (mark > best_mark)
        {
            best_mark = mark;
            GBest = G;
        }
    }
    //logerr("%d mutation. Best mark = %f",function_stat.tries, best_mark);
    G = GBest;
}

void GeneticAlgorithm::find_pairs(int cnt, std::vector<std::pair<int, int>> &pairs)
{
    float fitsum = 0;
    for (auto &c : population)
    {
        if (c.alive)
            fitsum += c.fitness;
    }
    std::vector<float> islands_fitsum(metaParams.n_islands, 0);
    for (auto &c : population)
    {
        if (c.alive)
            islands_fitsum[c.sub_population_n] += c.fitness;
    }
    float f_cnt = (float)cnt/metaParams.n_islands;
    float cnt_err = 0;
    for (int island=0;island<metaParams.n_islands;island++)
    {
        int i_cnt = f_cnt + cnt_err;
        cnt_err = f_cnt - i_cnt;
        for (int n = 0; n < i_cnt; n++)
        {
            int p0 = -1;
            int p1 = -1;

            bool found = false;
            int tries = 0;
            while (!found && tries < 100)
            {
                float v = urand(0, islands_fitsum[island]);
                for (int i = 0; i < population.size(); i++)
                {
                    if (population[i].alive && population[i].sub_population_n == island)
                    {
                        if (v < population[i].fitness)
                        {
                            p0 = i;
                            found = true;
                            break;
                        }
                        else
                        {
                            v -= population[i].fitness;
                        }
                    }
                }    
                tries++;           
            }
            found = false;
            while (!found && tries < 100)
            {
                float v = urand(0, islands_fitsum[island]);
                for (int i = 0; i < population.size(); i++)
                {
                    if (population[i].alive && population[i].sub_population_n == island)
                    {
                        if (v < population[i].fitness)
                        {
                            p1 = i;
                            found = true;
                            break;
                        }
                        else
                        {
                            v -= population[i].fitness;
                        }
                    }
                }
                tries++;
            }
            if (p0 >=0 && p1>=0)
                pairs.push_back(std::pair<int, int>(p0, p1));
        }
    }
    return;
    for (int n = 0; n < cnt; n++)
    {
        int p0 = -1;
        int p1 = -1;

        bool found = false;
        while (!found)
        {
            float v = urand(0, fitsum);
            for (int i = 0; i < population.size(); i++)
            {
                if (population[i].alive)
                {
                    if (v < population[i].fitness)
                    {
                        p0 = i;
                        found = true;
                        break;
                    }
                    else
                    {
                        v -= population[i].fitness;
                    }
                }
            }
        }
        int p0_island = population[p0].sub_population_n;
        found = false;
        while (!found)
        {
            float v = urand(0, islands_fitsum[p0_island]);
            for (int i = 0; i < population.size(); i++)
            {
                if (population[i].alive && population[i].sub_population_n == p0_island)
                {
                    if (v < population[i].fitness)
                    {
                        p1 = i;
                        found = true;
                        break;
                    }
                    else
                    {
                        v -= population[i].fitness;
                    }
                }
            }
        }
        pairs.push_back(std::pair<int, int>(p0, p1));
        //logerr("parents %d %d chosen", p0, p1);
    }
}

void GeneticAlgorithm::kill_creature(int n)
{
    population[n].alive = false;
    population[n].fitness = -1;
    current_population_size--;
}

void GeneticAlgorithm::kill_old()
{
    for (int i=0;i<population.size();i++)
    {
        if (population[i].age >= population[i].max_age && current_population_size > metaParams.min_population_size)
            kill_creature(i);
    }
}

void GeneticAlgorithm::kill_weak(int cnt_to_kill)
{
    std::sort(population.begin(), population.end(), 
              [&](const Creature& a, const Creature& b) -> bool
                 {
                     if (a.alive == b.alive)
                     {
                         if (a.sub_population_n == b.sub_population_n)
                            return a.fitness < b.fitness;
                         else
                            return a.sub_population_n < b.sub_population_n;
                     }
                     else
                        return a.alive > b.alive;
                 });
    
    int killed = 0;
    std::vector<int> islands_start = std::vector<int>(metaParams.n_islands,0);
    std::vector<int> islands_end = std::vector<int>(metaParams.n_islands,0);
    for (int i=0;i<population.size();i++)
    {
        //logerr("%d)%d island %d fitness %f",i, population[i].alive, population[i].sub_population_n, population[i].fitness);
        if (population[i].alive)
        {
            if (i == 0)
                islands_start[population[i].sub_population_n] = i;
            else if (population[i].sub_population_n != population[i-1].sub_population_n)
            {
                islands_start[population[i].sub_population_n] = i;
                islands_end[population[i-1].sub_population_n] = i-1;
            }
        }
        else if (i > 0 && population[i-1].alive)
            islands_end[population[i-1].sub_population_n] = i-1;
    }
    if (population.back().alive && islands_end[population.back().sub_population_n] == 0)
        islands_end[population.back().sub_population_n] = population.size() - 1;
    float percent = (float)cnt_to_kill/current_population_size;
    float np_err = 0;
    for (int i=0;i<metaParams.n_islands;i++)
    {
        //logerr("st en %d %d", islands_start[i], islands_end[i]);
        //float np_f = percent*(islands_end[i] - islands_start[i]) + np_err;
        float dest_cnt = metaParams.max_population_size/(2*metaParams.n_islands);
        int new_pop_end = MAX(islands_start[i], islands_end[i] - dest_cnt);
        //np_err += np_f - (new_pop_end - islands_start[i]);
        for (int j=islands_start[i];j<new_pop_end;j++)
        {
            kill_creature(j); 
            killed++;
            if (killed >= cnt_to_kill)
                break;
        }
        if (killed >= cnt_to_kill)
            break;
    }
}

void GeneticAlgorithm::crossingover(Genome &A, Genome &B, Genome &res)
{
    float rnd = urand();
    Genome &p1 = rnd < 0.5 ? A : B;
    Genome &p2 = rnd < 0.5 ? B : A; 
    res = p1;

        if (urand() < metaParams.mix_chance)
        {
            for (int i=0;i<A.size();i++)
            {
                if (urand() > 0.5)
                    res[i] = p2[i];
            }
        }
        else
        {
            float pos = urandi(1, A.size() - 1);
            for (int i=pos;i<A.size();i++)
            {
                res[i] = p2[i];
            }
        }
}

float GeneticAlgorithm::genes_dist(Creature &A, Creature &B)
{
    if (!A.alive || !B.alive || A.main_genome.size() != B.main_genome.size())
        return 1;
    float dist = 0;
    for (int i=0;i<A.main_genome.size();i++)
    {
        dist += abs(A.main_genome[i] - B.main_genome.size());
    }
    return dist/A.main_genome.size();
}

float GeneticAlgorithm::closest_neighbour(Creature &C, std::vector<Creature> &population)
{
    int i=0;
    float min_dist = 1;
    while (i < population.size() && population[i].alive)
    {
        if (population[i].id != C.id && population[i].sub_population_n == C.sub_population_n)
        {
            float dist =  genes_dist(C, population[i]);
            min_dist = MIN(min_dist,dist);
        }
        i++;
    }
    return min_dist;
}

void GeneticAlgorithm::make_child(Creature &A, Creature &B, Creature &C)
{
    if (A.sub_population_n != B.sub_population_n)
    {
        logerr("something went wrong %d %d", A.sub_population_n, B.sub_population_n);
    }
    C = Creature();
    C.alive = true;
    C.age = 0;
    C.max_age = metaParams.max_age;
    C.id = next_id.fetch_add(1);
    C.sub_population_n = A.sub_population_n;
    if (metaParams.evolution_stat)
    {
        stat.all_births.push_back(glm::ivec3(A.id, B.id, C.id));
    }

    if (metaParams.n_ploid_genes == 1)
    {
        crossingover(A.main_genome, B.main_genome, C.main_genome);
    }
    else if (metaParams.n_ploid_genes == 2)
    {
        C.other_genomes.push_back(A.other_genomes[0]);
        int c = rand() % 4;
        switch (c)
        {
        case 0:
            crossingover(A.main_genome, B.main_genome, C.main_genome);
            crossingover(A.other_genomes[0], B.other_genomes[0], C.other_genomes[0]);
            break;
        case 1:
            crossingover(A.main_genome, B.other_genomes[0], C.main_genome);
            crossingover(A.other_genomes[0], B.main_genome, C.other_genomes[0]);
            break;
        case 2:
            crossingover(A.other_genomes[0], B.main_genome, C.main_genome);
            crossingover(A.main_genome, B.other_genomes[0], C.other_genomes[0]);
            break;
        case 3:
            crossingover(A.other_genomes[0], B.other_genomes[0], C.main_genome);
            crossingover(A.main_genome, B.main_genome, C.other_genomes[0]);
            break;
        default:
            break;
        }
    }
    else
    {
        //it's shit
    }
    float it = 1 + 0.5*sqrt(iteration_n);
    float mutation_power = CLAMP(0.2*(1 + 0.1*sub_population_infos[C.sub_population_n].no_progress_time),0,1);
    mutation_power = 1 - MIN(A.metric, B.metric);
    //mutation_power = 0.15;
    if (urand() < metaParams.m_ch_min || iteration_n < 30)
    {
        mutation(C.main_genome, 1, urandi(1, 1 + metaParams.m_ch_genes*free_parameters_cnt));
    }
    for (auto &g : C.other_genomes)
    {
        if (urand() < 0.33/it)
            mutation(C.main_genome, 1.0, urandi(1, free_parameters_cnt));
        else
            mutation(C.main_genome, 1/it, urandi(1, 0.33*free_parameters_cnt));
    }


    while (closest_neighbour(C, new_population) < metaParams.clone_thr)
    {
        //logerr("remutation %d, it has clone in population",C.id);
        mutation(C.main_genome, 1/it, urandi(1, free_parameters_cnt));
        //mutation(C.main_genome, 0.33/it, urandi(1, 0.5*free_parameters_cnt));
    }
}

void GeneticAlgorithm::make_new_generation(std::vector<std::pair<int, int>> &pairs)
{
    int pos = 0;
    for (int i=0;i<population.size();i++)
    {
        if (population[i].alive)
        {
            new_population[pos] = population[i];
            pos++;
        }
    }

    for (auto &p : pairs)
    {
        //logerr("children %d %d", p.first, p.second);
        make_child(population[p.first], population[p.second], new_population[pos]);
        pos++;
        //if (pos >= metaParams.max_population_size)
        //    break;
    }
    current_population_size += pairs.size();
    population = new_population;
}

void GeneticAlgorithm::print_function_stat()
{
    logerr("d_xi: ");
    for (int k = 0; k < population[0].main_genome.size(); k++)
    {
        logerr("(%.3f, %.3f, %d)", (float)di_res[k] / MAX(1, di_cnt[k]), (float)di_dev[k] / MAX(1, di_cnt[k]), di_cnt[k]);
    }
}

void GeneticAlgorithm::calculate_metric(int heaven_n, bool elite_fine_tuning)
{
    std::vector<std::vector<float>> params = {};
    enum ParOrigin
    {
        NEW_ENTITY,
        ELITE_ORIGIN,
        ELITE_MODIFIED,
        CREATURE_MODIFIED
    };
    std::vector<std::pair<ParOrigin,int>> positions;
    int fine_tune_pos[100] = {-1};
    int i=0;
    for (auto &p : population)
    {
        if (p.alive)
        {
            if (p.metric_calc_n == 0)
            {
                params.push_back(p.main_genome);
                positions.push_back(std::pair<ParOrigin,int>(NEW_ENTITY, i));
            }
        }
        i++;
    }

    int recnt = 0;
    float delta = 0.2;
    for (int i=population.size()-1;i>=0;i--)
    {
        auto &p = population[i];
        if (p.alive)
        {
            if (p.metric >= 0 && iteration_n > 20 && recnt < 0)
            {
                Genome modified = p.main_genome; 
                mutation(modified, delta, 1, &(fine_tune_pos[recnt]));
                params.push_back(modified);
                positions.push_back(std::pair<ParOrigin, int>(CREATURE_MODIFIED, i));   
                recnt++;           
            }
        }
    }
    std::vector<float> metrics = opt_function.f(params);
    func_called += metrics.size();
    int ftp = 0;
    for (int i=0;i<metrics.size();i++)
    {
        if (positions[i].first == NEW_ENTITY)
        {
            auto &c = population[positions[i].second];
            c.metric = metrics[i];
            c.metric_calc_n++;
        }
        else if (positions[i].first == ELITE_ORIGIN)
        {
            auto &c = heaven[positions[i].second];
            c.metric = (metrics[i] + c.metric_calc_n*c.metric)/(c.metric_calc_n+1);
            c.metric_calc_n++;
        }
        else if (positions[i].first == CREATURE_MODIFIED)
        {
            auto &c = population[positions[i].second];
            int var_num = fine_tune_pos[ftp];
            if (var_num >= 0 && var_num < 256)
            {
                di_cnt[var_num]++;
                di_res[var_num]+= abs(c.metric - metrics[i])/delta;
                di_dev[var_num]+= SQR(abs(c.metric - metrics[i])/delta);
            }

            if (c.metric < metrics[i])
            {
                c.metric = metrics[i];
                c.main_genome = params[i];
                c.metric_calc_n++;
            }
            ftp++;
        }
    }
    for (int i=0;i<metrics.size();i++)
    {
        if (positions[i].first == ELITE_MODIFIED)
        {
            auto &c = heaven[positions[i].second];
            if (metrics[i] > c.metric)
            {
                //logerr("elite individual %d %d improved %f --> %f", positions[i].second, c.id, c.metric, metrics[i]);
                c.main_genome = params[i];
                c.metric = metrics[i];
                c.metric_calc_n = 3;
            }
        }
    }
}

bool GeneticAlgorithm::better(Creature &A, Creature &B)
{
    return A.metric > B.metric;
}

void GeneticAlgorithm::pick_best_to_heaven()
{
    best_metric_current = 0;
    for (int n=0;n<population.size();n++)
    {
        if (!population[n].alive)
            continue;
        best_metric_current = MAX(best_metric_current, population[n].metric);
        if (heaven.size() < MAX(metaParams.best_genoms_count, metaParams.heaven_size))
        {
            heaven.push_back(population[n]);
            best_metric_ever = MAX(best_metric_ever, population[n].metric);
        }
        else
        {
            int worst_pos = -1;
            float worst_val = 1e9;
            bool same = false;

            //creature should be better that someone in heaven
            for (int i=0;i<heaven.size();i++)
            {
                if (heaven[i].id == population[n].id)
                {
                    same = true;
                    worst_pos = i;
                    break;
                }
                if (heaven[i].metric < population[n].metric && heaven[i].metric < worst_val)
                {
                   worst_val = heaven[i].metric;
                   worst_pos = i;
                }
            }
            if (!same && worst_pos >= 0)
            {
                heaven[worst_pos] = population[n];
            }
            else if (same && worst_pos >= 0 && heaven[worst_pos].metric < population[n].metric)
            {
                heaven[worst_pos] = population[n];
            }
        }
    }
    
    best_metric_ever = 0;
    for (auto &c : heaven)
    {
        best_metric_ever = MAX(best_metric_ever, c.metric);
    }
    if (metaParams.evolution_stat)
    {
        for (int i=0;i<population.size();i++)
        {
            auto &c = population[i];
            if (c.alive)
            logerr("%d metric[%d] = %.3f %.3f(%d calc)",c.id, i, c.metric, c.fitness, c.metric_calc_n);
        }
        for (int i=0;i<heaven.size();i++)
        {
            auto &c = heaven[i];
            logerr("%d metric[%d] = %.3f (%d calc)",c.id, -i-1, c.metric, c.metric_calc_n);
        }
    }
   if (metaParams.evolution_stat)
   {
       for (auto &c : population)
       {
           if (c.alive)
           {
               stat.all_results.emplace(c.id, c.metric);
           }
       }
   }
}

void GeneticAlgorithm::migration()
{
    if (metaParams.n_islands == 1)
        return;
    for (auto &p : population)
    {
        if (urand() < metaParams.migration_chance)
        {
            p.sub_population_n =  (p.sub_population_n+1)% metaParams.n_islands;
        }
    }
}

void GeneticAlgorithm::recalculate_fitness()
{
    std::sort(population.begin(), population.end(), 
              [&](const Creature& a, const Creature& b) -> bool{return a.metric < b.metric;});
    for (auto &p : population)
    {
        if (p.alive && p.metric < metaParams.dead_at_birth_thr && current_population_size > metaParams.min_population_size)
        {
            p.alive = false;
            current_population_size--;
        }
        if (p.alive)
            p.fitness = pow(p.metric, 1 + 0.1*iteration_n) + 1e-4;
        else
            p.fitness = -1;
    }
}

void GeneticAlgorithm::save_load_function_stat(bool load)
{
    
    std::string stat_block_name = "GA_function_stats.blk";
    std::string function_stat_block_name = "function_values_stat.blk";
    if (load == true)
    {
        Block f_stat_block;
        load_block_from_file(function_stat_block_name, f_stat_block);
        Block *b2 = f_stat_block.get_block(opt_function.name);
        function_stat = my_opt::FunctionStat(all_parameters_cnt);
        if (b2 && opt_function.version == b2->get_int("version",-1))
            function_stat.save_load_blk(*b2, false);
    }
    Block stat_block;
    load_block_from_file(stat_block_name, stat_block);
    Block *bl = stat_block.get_block(opt_function.name);
    if (!bl)
    {
        debug("creating new GA function stat block for function [%s]\n", opt_function.name.c_str());
        bl = new Block();
        stat_block.add_block(opt_function.name, bl);
    }
    int version = bl->get_int("version",-1);
    int var_cnt = all_parameters_cnt;
    std::vector<int> di_cnt_h;
    std::vector<float> di_res_h;
    if (version != opt_function.version)
    {
        debug("new version of GA function %d --> %d\n", version, opt_function.version);
        bl->set_int("version", opt_function.version);
        di_cnt_h = std::vector<int>(var_cnt, 0);
        di_res_h = std::vector<float>(var_cnt, 0);
    }
    else
    {
        bl->get_arr("di_cnt", di_cnt_h);
        bl->get_arr("di_res", di_res_h);
        if (di_cnt_h.size() != var_cnt || di_res_h.size() != var_cnt)
        {
            debug("reset function stat %d %d %d\n", di_cnt_h.size(), di_res_h.size(), var_cnt);
            di_cnt_h = std::vector<int>(var_cnt, 0);
            di_res_h = std::vector<float>(var_cnt, 0);
        }
    }
    if (load)
    {
        for (int i=0;i<var_cnt;i++)
        {
            int cnt = di_cnt_h[i] + di_cnt[i];
            di_res[i] = (double)(di_res_h[i]*di_cnt_h[i] + di_res[i]);
            di_cnt[i] = cnt;
        }
    }
    else
    {
        for (int i=0;i<var_cnt;i++)
        {
            di_cnt_h[i] = di_cnt[i];
            di_res_h[i] = di_res[i]/MAX(di_cnt[i], 1);
        }
        bl->set_arr("di_cnt", di_cnt_h);
        bl->set_arr("di_res", di_res_h);

        save_block_to_file(stat_block_name, stat_block);
    }
}