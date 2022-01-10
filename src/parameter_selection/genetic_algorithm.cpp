#include "genetic_algorithm.h"
#include <algorithm>

void GeneticAlgorithm::perform(ParameterList &param_list, MetaParameters params, ExitConditions exit_conditions,
                               const std::function<std::vector<float>(std::vector<ParameterList> &)> &f,
                               std::vector<std::pair<float, ParameterList>> &best_results)
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

    initialize_population();
    calculate_metric();
    recalculate_fitness();

    while (!should_exit())
    {
        kill_old();
        kill_weak(metaParams.weaks_to_kill);
        if (current_population_size == 0)
        {
            logerr("population extincted");
            return;
        }
        else if (current_population_size < 0.1*metaParams.initial_population_size)
        {
            logerr("population %d is about to extinct", current_population_size);
        }
        int space_left = metaParams.max_population_size - current_population_size;
        std::vector<std::pair<int, int>> pairs;
        find_pairs(0.67*space_left, pairs);
        for (auto &p : pairs)
        {
            logerr("children %d %d", p.first, p.second);
            make_children(population[p.first], population[p.second], 1);
        }
        
        calculate_metric();
        recalculate_fitness();
        for (int i=0;i<population.size();i++)
        {
            if (population[i].alive)
                population[i].age++;
        }
        iteration_n++;
    }

    prepare_best_params(best_results);
}

bool GeneticAlgorithm::should_exit()
{
    iter_start = std::chrono::steady_clock::now();
    float time = std::chrono::duration_cast<std::chrono::milliseconds>(iter_start - t_start).count();
    time_spent_sec += 0.001*time;
    return (time_spent_sec > exitConditions.time_elapsed_seconds || iteration_n >= exitConditions.generations
    || func_called >= exitConditions.function_calculated || best_metric_ever >= exitConditions.function_reached);
}

void GeneticAlgorithm::prepare_best_params(std::vector<std::pair<float, ParameterList>> &best_results)
{
    for (int i=0;i<population.size();i++)
    {
        if (population[i].alive)
            kill_creature(i);
    }
    logerr("heaven popultion %d", heaven.size());
    for (auto &creature : heaven)
    {
        best_results.emplace_back();
        best_results.back().first = creature.metric;
        best_results.back().second = original_param_list;
        best_results.back().second.from_simple_list(creature.main_genome);
        best_results.back().second.print();
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

void GeneticAlgorithm::initialize_population()
{
    population = std::vector<Creature>(metaParams.max_population_size, Creature());
    for (int i=0;i<metaParams.initial_population_size;i++)
    {
        population[i].main_genome = random_genes();
        for (int j=1;j<metaParams.n_ploid_genes;j++)
            population[i].other_genomes.push_back(random_genes());
        population[i].alive = true;
        population[i].max_age = metaParams.max_age;
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
            //found = true;

            int pos = urandi(0, G.size());
            logerr("mutation pos %d G size %d %f", pos, G.size(), (float)urandi(0, G.size()));
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
                    float val = G[g_pos] + len*urand(-1,1);
                    val = CLAMP(val, parametersMask.ordinalParameters[pos].second.min_val, parametersMask.ordinalParameters[pos].second.max_val);
                    G[g_pos] = val;
                    found = true;
                }
            }
            else
            {
                pos -= parametersMask.categorialParameters.size() + parametersMask.ordinalParameters.size();
                logerr("continuous pos %d mutation", pos);
                if (!parametersMask.continuousParameters[pos].second.fixed())
                {
                    float len = mutation_power*(parametersMask.continuousParameters[pos].second.max_val - parametersMask.continuousParameters[pos].second.min_val);
                    float val = G[g_pos] + len*urand(-1,1);
                    val = CLAMP(val, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
                    G[g_pos] = val;
                    logerr("applied pos val %d %f in [%f %f]",pos, val, parametersMask.continuousParameters[pos].second.min_val, parametersMask.continuousParameters[pos].second.max_val);
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
    int vals[2];
    for (int n=0;n<cnt;n++)
    {
        vals[0] = -1;
        for (int k=0;k<2;k++)
        {
            bool found = false;
            while (!found)
            {
                float v = urand(0, fitsum);
                for (int i=0;i<population.size();i++)
                {
                    if (population[i].alive)
                    {
                        if (v < population[i].fitness && i != vals[0])
                        {
                            vals[k] = i;
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
        }
        pairs.push_back(std::pair<int, int>(vals[0], vals[1]));
        logerr("parents %d %d chosen", vals[0], vals[1]);
    }
}

void GeneticAlgorithm::kill_creature(int n)
{
    if (heaven.size() < metaParams.best_genoms_count)
    {
        heaven.push_back(population[n]);
    }
    else
    {
        for (int i=0;i<heaven.size();i++)
        {
            if (heaven[i].metric < population[n].metric)
            {
                heaven[i] = population[n];
                break;
            }
        }
    }
    population[n].alive = false;
    population[n].fitness = -1;
    current_population_size--;
}

void GeneticAlgorithm::kill_old()
{
    for (int i=0;i<population.size();i++)
    {
        if (population[i].age >= population[i].max_age)
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
        logerr("%d) %d fitness %f",i,population[i].alive, population[i].fitness);
    }
    for (int i=0;i<population.size();i++)
    {
        if (population[i].alive)
        {
           logerr("killed weak %d metric %f", i, population[i].metric);
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
        if (urand() > 0.5)
            res[i] = B[i];
    }
}

void GeneticAlgorithm::make_children(Creature &A, Creature &B, int count)
{
    for (int ch=0;ch<count;ch++)
    {
        int pos = -1;
        for (int i=0;i<population.size();i++)
        {
            if (!population[i].alive)
            {
                pos = i;
                population[i] = Creature();
                break;
            }
        }
        population[pos].alive = true;
        population[pos].age = 0;
        population[pos].max_age = metaParams.max_age;
        logerr("child %d created", pos);
        current_population_size++;

        if (metaParams.n_ploid_genes == 1)
        {
            crossingover(A.main_genome, B.main_genome, population[pos].main_genome);
        }
        else
        {
            //TODO: implement diploid genomes
        }
        mutation(population[pos].main_genome, 0.4, urandi(1, 0.5*free_parameters_cnt));
        for (auto &g : population[pos].other_genomes)
        {
            mutation(g, 0.5, urandi(1, 0.5*free_parameters_cnt));
        }
    }
}

void GeneticAlgorithm::calculate_metric()
{
    std::vector<ParameterList> params;
    std::vector<int> positions;
    int i=0;
    for (auto &p : population)
    {
        if (p.alive && p.metric < 0)
        {
            params.push_back(original_param_list);
            params.back().from_simple_list(p.main_genome);
            positions.push_back(i);
        }
        i++;
    }

    std::vector<float> metrics = function(params);
    for (int i=0;i<metrics.size();i++)
    {
        population[positions[i]].metric = metrics[i];
        logerr("metric[%d] = %f", positions[i], metrics[i]);
    }
}
void GeneticAlgorithm::recalculate_fitness()
{
    for (auto &p : population)
    {
        if (p.alive && p.metric < metaParams.dead_at_birth_thr)
        {
            p.alive = false;
            current_population_size--;
        }
        if (p.alive)
            p.fitness = p.metric;
        else
            p.fitness = -1;
    }
}