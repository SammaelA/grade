#pragma once
#include "impostor_similarity.h"
#include "tree_generators/abstract_generator.h"
#include "generic_optimization_algorithm.h"
#include "function_stat.h"
#include <chrono>

class GeneticAlgorithm : public my_opt::Optimizer
{
public:
    enum GA_Type
    {
        ISLANDS_GA,
        TREE_GA
    };
    struct MetaParameters : public my_opt::MetaParameters
    {
        virtual void RW_vector(std::vector<float> &res) {};
        int initial_population_size = 100;
        int max_population_size = 50;
        int min_population_size = 5;
        int best_genoms_count = 3;
        int heaven_size = 1;
        int heaven_recalc_n = 4;
        int heaven_fine_tuning_count = 0;
        int elite = 5;
        float weaks_to_kill = 0.5;
        float dead_at_birth_thr = -1e9;
        float mix_chance = 1;
        int n_ploid_genes = 1;
        int max_age = 100;
        float clone_thr = 0.01;
        int n_islands = 1;
        int migration_interval = 5;
        float migration_chance = 0.2;
        int tree_GA_iters = 75;
        float m_ch_min = 0.5;
        float m_ch_genes = 0.15;
        bool use_genes_importance = true;
        bool use_function_stat = true;
        GA_Type type = GA_Type::TREE_GA;
        bool evolution_stat = false;
        bool debug_graph = false;
    };
    void perform(std::vector<float> &param_list, my_opt::MetaParameters *params, my_opt::ExitConditions exit_conditions,
                 const my_opt::OptFunction &my_opt_f,
                 std::vector<std::pair<float,std::vector<float>>> &best_results,
                 std::vector<std::vector<float>> &initial_types);
    virtual ~GeneticAlgorithm() {};
private:
    my_opt::OptFunction opt_function;
    MetaParameters metaParams;
    my_opt::ExitConditions exitConditions;
    int free_parameters_cnt = 0;
    int all_parameters_cnt = 0;
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
        int sub_population_n = -1;
    };
    struct Stat
    {
        std::map<int, float> all_results;
        std::vector<glm::ivec3> all_births;//mother, father, child
    } stat;

    struct PopulationBackup
    {
        std::vector<Creature> pop;
        float best_value = 0;
        glm::vec4 model = glm::vec4(0,0,0,0);
        float predicted_best_value = 0;
        int backup_uses = 0;
    };

    struct SubPopulationInfo
    {
        static constexpr float progress_thr = 0.01;
        std::vector<float> best_results_history;
        std::vector<PopulationBackup> backups;
        int best_result_iter = 0;
        int no_progress_time = 0;
    };

    std::vector<Creature> population;
    std::vector<Creature> new_population;
    std::vector<Creature> heaven;
    std::vector<Creature> crio_camera;
    std::vector<glm::vec2> TreeGA_stat;

    std::chrono::steady_clock::time_point t_start, iter_start;
    int current_population_size = 0;
    float best_metric_ever = 0;
    float best_metric_current = 0;
    float time_spent_sec = 0; 
    int func_called = 0;
    int iteration_n = 0;

    bool should_exit();
    void prepare_best_params(std::vector<std::pair<float,std::vector<float>>> &best_results);
    void initialize_population(std::vector<std::vector<float>> &initial_types);
    float get_mark(Genome &G);
    void mutation(Genome &G, float mutation_power, int mutation_genes_count, int *single_mutation_pos = nullptr);
    void find_pairs(int cnt, std::vector<std::pair<int, int>> &pairs);
    void kill_creature(int n);
    void kill_old();
    void kill_weak(int count);
    void crossingover(Genome &A, Genome &B, Genome &res);
    void make_new_generation(std::vector<std::pair<int, int>> &pairs);
    void make_child(Creature &A, Creature &B, Creature &C);
    void calculate_metric(int heaven_n = 1, bool elite_fine_tuning = true);
    void recalculate_fitness();
    void pick_best_to_heaven();
    void migration();
    Genome random_genes();
    float closest_neighbour(Creature &C, std::vector<Creature> &population);
    float genes_dist(Creature &A, Creature &B);
    bool better(Creature &A, Creature &B);
    void print_function_stat();
    void save_load_function_stat(bool load);
    void islands_GA(std::vector<std::vector<float>> &initial_types);
    void tree_GA(std::vector<std::vector<float>> &initial_types);
    void tree_GA_internal(int depth, int iters, int width, std::vector<std::vector<float>> &initial_types, PopulationBackup &result);

    std::vector<SubPopulationInfo> sub_population_infos;
    my_opt::FunctionStat function_stat;
};