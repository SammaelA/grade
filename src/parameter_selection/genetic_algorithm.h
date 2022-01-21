#pragma once
#include "impostor_similarity.h"
#include "tree_generators/abstract_generator.h"
#include <chrono>

class GeneticAlgorithm
{
public:
    struct MetaParameters
    {
        int initial_population_size = 100;
        int max_population_size = 50;
        int min_population_size = 10;
        int best_genoms_count = 3;
        int heaven_size = 7;
        int heaven_recalc_n = 4;
        int elite = 5;
        float weaks_to_kill = 0.5;
        float dead_at_birth_thr = 0.001;
        float mix_chance = 0.25;
        int n_ploid_genes = 1;
        int max_age = 100;
        bool evolution_stat = false;
        bool debug_graph = false;
    };
    struct ExitConditions
    {
        float time_elapsed_seconds = 45*60;
        float function_reached = 1;
        int function_calculated = 10000;
        int generations = 10000;
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
        int metric_calc_n = 0;
        int age = 0;
        int max_age = 1000;
        int children_cnt = 0;
        bool alive = false;
        int id = 0;
    };
    std::vector<Creature> population;
    std::vector<Creature> new_population;
    std::vector<Creature> heaven;
    std::vector<Creature> crio_camera;


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
    void make_new_generation(std::vector<std::pair<int, int>> &pairs);
    void make_child(Creature &A, Creature &B, Creature &C);
    void calculate_metric(int heaven_n = 1);
    void recalculate_fitness();
    Genome random_genes();

    struct Stat
    {
        std::map<int, float> all_results;
        std::vector<glm::ivec3> all_births;//mother, father, child
    } stat;
};