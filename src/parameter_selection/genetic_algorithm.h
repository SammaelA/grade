#pragma once
#include "impostor_similarity.h"
#include "tree_generators/abstract_generator.h"
#include <chrono>

class GeneticAlgorithm
{
public:
    struct MetaParameters
    {
        int initial_population_size = 40;
        int max_population_size = 50;
        int best_genoms_count = 1;
        float weaks_to_kill = 0.33;
        float dead_at_birth_thr = 0.0;
        int n_ploid_genes = 1;
        int max_age = 100;
    };
    struct ExitConditions
    {
        float time_elapsed_seconds = 60*60;
        float function_reached = 1;
        int function_calculated = 1000;
        int generations = 30;
    };
    void perform(ParameterList &param_list, MetaParameters params, ExitConditions exit_conditions,
                 const std::function<std::vector<float>(std::vector<ParameterList> &)> &f,
                 std::vector<std::pair<float,ParameterList>> &best_results);
private:
    std::function<std::vector<float>(std::vector<ParameterList> &)> function;
    ParameterList original_param_list;
    struct ParametersMask
    {
        std::vector<std::pair<std::string, CategorialParameter>> categorialParameters;
        std::vector<std::pair<std::string, OrdinalParameter>> ordinalParameters;
        std::vector<std::pair<std::string, ContinuousParameter>> continuousParameters;
    } parametersMask;
    MetaParameters metaParams;
    ExitConditions exitConditions;
    int free_parameters_cnt = 0;

    typedef std::vector<float> Genome;
    struct Creature
    {
        Genome main_genome;
        std::vector<Genome> other_genomes;//recessive
        float metric = -1;
        float fitness = -1;
        int age = 0;
        int max_age = 1000;
        int children_cnt = 0;
        bool alive = false;
    };
    std::vector<Creature> population;
    std::vector<Creature> heaven;

    std::chrono::steady_clock::time_point t_start, iter_start;
    int current_population_size = 0;
    float best_metric_ever = 0;
    float time_spent_sec = 0; 
    int func_called = 0;
    int iteration_n = 0;

    bool should_exit();
    void prepare_best_params(std::vector<std::pair<float,ParameterList>> &best_results);
    void initialize_population();
    void mutation(Genome &G, float mutation_power, int mutation_genes_count);
    void find_pairs(int cnt, std::vector<std::pair<int, int>> &pairs);
    void kill_creature(int n);
    void kill_old();
    void kill_weak(float percent);
    void crossingover(Genome &A, Genome &B, Genome &res);
    void make_children(Creature &A, Creature &B, int count);
    void calculate_metric();
    void recalculate_fitness();
    Genome random_genes();
};