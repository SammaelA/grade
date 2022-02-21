#include "genetic_algorithm.h"
#include <algorithm>
#include <atomic>
#include "GA_utils.h"
#include <set>

std::atomic<int> next_id(0);

void GeneticAlgorithm::perform(ParameterList &param_list, MetaParameters params, ExitConditions exit_conditions,
                               const std::function<std::vector<float>(std::vector<ParameterList> &)> &f,
                               std::vector<std::pair<float, ParameterList>> &best_results,
                               std::vector<ParameterList> &initial_types)
{
    t_start = std::chrono::steady_clock::now();

    function = f;
    original_param_list = param_list;
    metaParams = params;
    exitConditions = exit_conditions;
    free_parameters_cnt = 0;

    for (auto &p : original_param_list.categorialParameters)
    {
        parametersMask.categorialParameters.push_back(p);
        if (!p.second.fixed())
            free_parameters_cnt++;
    }
    for (auto &p : original_param_list.ordinalParameters)
    {
        parametersMask.ordinalParameters.push_back(p);
        if (!p.second.fixed())
            free_parameters_cnt++;
    }
    for (auto &p : original_param_list.continuousParameters)
    {
        parametersMask.continuousParameters.push_back(p);
        if (!p.second.fixed())
            free_parameters_cnt++;
    }

    metaParams.heaven_size = MAX(metaParams.heaven_size, metaParams.best_genoms_count + 1);
    metaParams.min_population_size = MAX(MAX(metaParams.min_population_size,metaParams.elite), 
                                     MAX(metaParams.heaven_size, metaParams.best_genoms_count));
    metaParams.max_population_size = MAX(metaParams.max_population_size, metaParams.min_population_size);
    next_id.store(0);
    initialize_population(initial_types);
    calculate_metric(1, false);
    recalculate_fitness();
    pick_best_to_heaven();
    debug("iteration 0 Pop: %d Best: %.4f\n", current_population_size, best_metric_ever);
    while (!should_exit())
    {
        next_id.store(1000*(iteration_n + 1));
        kill_old();
        kill_weak(metaParams.weaks_to_kill);
        if (current_population_size < 5)
        {
            logerr("population extincted");
            return;
        }
        else if (current_population_size < metaParams.min_population_size)
        {
            logerr("population %d is about to extinct", current_population_size);
        }
        int space_left = metaParams.max_population_size - heaven.size();
        std::vector<std::pair<int, int>> pairs;
        find_pairs(space_left, pairs);
        make_new_generation(pairs);
        
        calculate_metric(1, true);
        recalculate_fitness();
        pick_best_to_heaven();
        if (iteration_n % metaParams.migration_interval == 0)
            migration();
        for (int i=0;i<population.size();i++)
        {
            if (population[i].alive)
                population[i].age++;
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

    prepare_best_params(best_results);
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

void GeneticAlgorithm::prepare_best_params(std::vector<std::pair<float, ParameterList>> &best_results)
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
        best_results.back().second = original_param_list;
        best_results.back().second.from_simple_list(heaven[i].main_genome);
        //best_results.back().second.print();
    }
    for (auto &creature : crio_camera)
    {
        best_results.emplace_back();
        best_results.back().first = creature.metric;
        best_results.back().second = original_param_list;
        best_results.back().second.from_simple_list(creature.main_genome);
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
    Genome g;
    for (auto &p : parametersMask.categorialParameters)
        g.push_back(p.second.possible_values[(int)urandi(0, p.second.possible_values.size())]);
    for (auto &p : parametersMask.ordinalParameters)
        g.push_back((int)(p.second.min_val + urand()*(p.second.max_val - p.second.min_val)));
    for (auto &p : parametersMask.continuousParameters)
        g.push_back((p.second.min_val + urand()*(p.second.max_val - p.second.min_val)));
    return g;
}

void GeneticAlgorithm::initialize_population(std::vector<ParameterList> &initial_types)
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
            t.to_simple_list(population[i0].main_genome);
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

void GeneticAlgorithm::mutation(Genome &G, float mutation_power, int mutation_genes_count)
{
    for (int gene = 0; gene < mutation_genes_count; gene++)
    {
        bool found = false;
        int tries = 0;
        while (!found && tries < 100)
        {
            int pos = urandi(0, G.size());
            int g_pos = pos;
            if (pos < parametersMask.categorialParameters.size())
            {
                if (!parametersMask.categorialParameters[pos].second.fixed())
                {
                    if (urand() < mutation_power)
                    {
                        int val_pos = urandi(0, parametersMask.categorialParameters[pos].second.possible_values.size());
                        G[g_pos] = parametersMask.categorialParameters[pos].second.possible_values[val_pos];
                    }
                    found = true;
                }
            }
            else if (pos < parametersMask.categorialParameters.size() + parametersMask.ordinalParameters.size())
            {
                pos -= parametersMask.categorialParameters.size();
                if (!parametersMask.ordinalParameters[pos].second.fixed())
                {
                    float len = mutation_power*(parametersMask.ordinalParameters[pos].second.max_val - parametersMask.ordinalParameters[pos].second.min_val);
                    float from = CLAMP(G[g_pos] - len, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    float to = CLAMP(G[g_pos] + len, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    float val = urand(from,to);
                    val = CLAMP(val, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    G[g_pos] = val;
                    found = true;
                }
            }
            else
            {
                pos -= parametersMask.categorialParameters.size() + parametersMask.ordinalParameters.size();
                //logerr("continuous pos %d mutation", pos);
                if (!parametersMask.continuousParameters[pos].second.fixed())
                {
                    float len = mutation_power*(parametersMask.continuousParameters[pos].second.max_val - parametersMask.continuousParameters[pos].second.min_val);
                    float from = CLAMP(G[g_pos] - len, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    float to = CLAMP(G[g_pos] + len, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    float val = urand(from,to);
                    val = CLAMP(val, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    G[g_pos] = val;
                    //logerr("applied pos val %d %f in [%f %f]",pos, val, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    found = true;
                }
            }
            tries++;
        }
    }
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

void GeneticAlgorithm::kill_weak(float percent)
{
    std::sort(population.begin(), population.end(), 
              [&](const Creature& a, const Creature& b) -> bool{return a.fitness < b.fitness;});
    
    int cnt_to_kill = percent*current_population_size;
    int killed = 0;
    for (int i=0;i<population.size();i++)
    {
        //logerr("%d) %d fitness %f",i,population[i].alive, population[i].fitness);
    }
    for (int i=0;i<population.size();i++)
    {
        if (population[i].alive)
        {
           //logerr("killed weak %d metric %f", i, population[i].metric);
           kill_creature(i); 
           killed++;
        }
        if (killed >= cnt_to_kill)
            break;
    }
}

void GeneticAlgorithm::crossingover(Genome &A, Genome &B, Genome &res)
{
    res = A;
    for (int i=0;i<A.size();i++)
    {
        if (urand() > metaParams.mix_chance)
        {
            if (urand() > 0.5)
                res[i] = B[i];
        }
        else
        {
            float rnd = urand();
            res[i] = rnd*A[i] + (1-rnd)*B[i];
        }
    }
}

float GeneticAlgorithm::genes_dist(Creature &A, Creature &B)
{
    if (!A.alive || !B.alive)
        return 1;
    ParameterList par1 = original_param_list;
    par1.from_simple_list(A.main_genome);
    ParameterList par2 = original_param_list;
    par2.from_simple_list(B.main_genome);
    return par1.diff(par2); 
}

float GeneticAlgorithm::closest_neighbour(Creature &C, std::vector<Creature> &population)
{
    int i=0;
    float min_dist = 1;
    while (population[i].alive && i < population.size())
    {
        if (population[i].id != C.id)
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
    float strong_chance = 0.4 - 0.35*(MAX(0.1*iteration_n*iteration_n,100)/100);
    float medium_chance = 0.4 - 0.3*(MAX(iteration_n - 25,100)/100);
    float rnd = urand();
    if (rnd < strong_chance)
        mutation(C.main_genome, 1.0, 0.5*free_parameters_cnt);
    else if (rnd < strong_chance + medium_chance)
        mutation(C.main_genome, 0.75, 5);
    else
        mutation(C.main_genome, 0.5, 3);
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
    //assume that the population is sorted by fitness;
    int elite_count = metaParams.elite > 0 ? heaven.size() : 0;
    if (elite_count)
    {
        for (int i=0;i<heaven.size();i++)
        {
            new_population[i] = heaven[i];
        }
    }
    int pos = elite_count;
    for (auto &p : pairs)
    {
        //logerr("children %d %d", p.first, p.second);
        make_child(population[p.first], population[p.second], new_population[pos]);
        pos++;
    }

    current_population_size = elite_count + pairs.size();
    population = new_population;
}

void GeneticAlgorithm::calculate_metric(int heaven_n, bool elite_fine_tuning)
{
    std::vector<ParameterList> params = {};
    enum ParOrigin
    {
        NEW_ENTITY,
        ELITE_ORIGIN,
        ELITE_MODIFIED
    };
    std::vector<std::pair<ParOrigin,int>> positions;

    for (int i=0;i<heaven_n;i++)
    {
        int h = 0;
        for (auto &p : heaven)
        {
            if (p.metric_calc_n < 1 || !elite_fine_tuning)
            {
                params.push_back(original_param_list);
                params.back().from_simple_list(p.main_genome);
                positions.push_back(std::pair<ParOrigin, int>(ELITE_ORIGIN, h));
            }
            if (elite_fine_tuning)
            {
                int ft_cnt = 0;
                if (iteration_n < 10)
                    ft_cnt = 0;
                else if (iteration_n < 40)
                    ft_cnt = 1;
                else if (iteration_n < 100)
                    ft_cnt = 2;
                else
                    ft_cnt = 3;
                ft_cnt *= metaParams.heaven_fine_tuning_count;
                for (int j = 0; j < ft_cnt; j++)
                {
                    Genome modified = p.main_genome;
                    mutation(modified, 0.25, 3);
                    params.push_back(original_param_list);
                    params.back().from_simple_list(modified);
                    positions.push_back(std::pair<ParOrigin, int>(ELITE_MODIFIED, h));
                }
            }
            h++;
        }
    }

    int i=0;
    for (auto &p : population)
    {
        if (p.alive && p.metric < 0)
        {
            params.push_back(original_param_list);
            params.back().from_simple_list(p.main_genome);
            positions.push_back(std::pair<ParOrigin,int>(NEW_ENTITY, i));
        }
        i++;
    }
    std::vector<float> metrics = function(params);
    func_called += metrics.size();
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
    }
    for (int i=0;i<metrics.size();i++)
    {
        if (positions[i].first == ELITE_MODIFIED)
        {
            auto &c = heaven[positions[i].second];
            if (metrics[i] > c.metric)
            {
                //logerr("elite individual %d %d improved %f --> %f", positions[i].second, c.id, c.metric, metrics[i]);
                params[i].to_simple_list(c.main_genome);
                c.metric = metrics[i];
                c.metric_calc_n = 3;
            }
        }
    }

    for (auto &p : population)
    {
        for (auto &h : heaven)
        {
            if (p.id == h.id)
            {
                p.metric = h.metric;
                p.metric_calc_n = h.metric_calc_n;
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

    for (int n=0;n<population.size();n++)
    {
        if (!population[n].alive)
            continue;
        if (heaven.size() < MAX(metaParams.best_genoms_count, metaParams.heaven_size))
        {
            heaven.push_back(population[n]);
            best_metric_ever = MAX(best_metric_ever, population[n].metric);
        }
        else
        {
            int worst_pos = -1;
            float worst_val = 1;
            bool same = false;

            //creature should be better that someone in heaven
            for (int i=0;i<heaven.size();i++)
            {
                if (heaven[i].id == population[n].id)
                {
                    same = true;
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
                //but also it should be better than all it's clones
                std::vector<int> clones;
                bool better_all_clones = true;
                for (int i=0;i<heaven.size();i++)
                {
                    if (genes_dist(population[n], heaven[i]) < 2*metaParams.clone_thr)
                    {
                        clones.push_back(i);
                        if (!better(population[n], heaven[i]))
                        {
                            better_all_clones = false;
                            break;
                        }
                    }
                }
                if (better_all_clones)
                {
                    //no need to keep clones that are worse than our creature
                    for (int &c_n : clones)
                    {
                        heaven[c_n] = Creature();
                    }
                }
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
            p.sub_population_n = urandi(0, metaParams.n_islands);
        }
    }
}

void GeneticAlgorithm::recalculate_fitness()
{
    std::sort(population.begin(), population.end(), 
              [&](const Creature& a, const Creature& b) -> bool{return a.metric < b.metric;});
    for (auto &p : population)
    {
        //logerr("%d(%d) metr %f", p.id, p.sub_population_n, p.metric);
        if (p.alive && p.metric < metaParams.dead_at_birth_thr && current_population_size > metaParams.min_population_size)
        {
            p.alive = false;
            current_population_size--;
        }
        if (p.alive)
            p.fitness = pow(p.metric, 1 + 0.2*sqrt(iteration_n)) + 1e-4;
        else
            p.fitness = -1;
    }

    float min_v = 0.01;
    float step = 4.0 / SQR(metaParams.max_population_size);
    int n = 0;
    for (int i=population.size()-1;i>=0;i--)
    {
        if (population[i].alive)
        {
            population[i].fitness = MAX(min_v, 1 - step*n*n);
            if (population[i].metric_calc_n == 1)
                n++;
        }
    }

    /*
    std::sort(population.begin(), population.end(), 
              [&](const Creature& a, const Creature& b) -> bool{return a.fitness < b.fitness;});
    float min_v = 0.1;
    float step = 0.025;
    int n = 0;

    for (auto &p : population)
    {
        if (p.alive && p.fitness > 0)
        {
            p.fitness = MAX(min_v, 1 - step*n);
            n++;
        }
    }
    */
   /*
   int p_cnt = 0;
   double p_dist = 0;
   debug("dist matrix:\n");
   for (auto &p1 : heaven)
   {
       if (p1.alive)
       {
            for (auto &p2 : heaven)
            {
                if (p2.alive)
                {
                    p_cnt++;
                    ParameterList par1 = original_param_list;
                    par1.from_simple_list(p1.main_genome);
                    ParameterList par2 = original_param_list;
                    par2.from_simple_list(p2.main_genome);
                    float dist =  par1.diff(par2);
                    p_dist += dist;
                    debug("%.2f ", dist);
                }
            }
            debugnl();
        }
   }
   logerr("average population diff %f", (float)(p_dist/p_cnt));
   */
}