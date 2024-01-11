#pragma once
#include "tensor_processor.h"
#include <cassert>
#include <map>
#include <string>

namespace nn
{
  struct TensorToken;
  class TensorCompiler
  {
  public:
    friend struct TensorToken;
    constexpr static int MAX_DIM = 4;

    void start_program();
    TensorProgram finish_program(bool print_program = false);
    void input(const TensorToken &t, std::string name);
    void output(const TensorToken &t, std::string name);

  private:
    struct Variable
    {
      unsigned Dim;
      unsigned offset;
      unsigned total_size;
      unsigned sizes[4];

      unsigned cmd_start;
      unsigned cmd_end;
      bool is_input = false;
      bool is_output = false;

      bool is_alias = false;
      unsigned alias_master_id = 0;
      unsigned alias_range_from = 0;
      unsigned alias_range_to = 0;
      std::vector<unsigned> aliases;
    };

    static bool have_same_scheme(const Variable &a, const Variable &b);
    static bool is_self_applicable_command(TensorProgram::CommandType type);

    unsigned add_var(const TensorToken &t);
    void ftt(unsigned id, float val);
    void add_command(TensorProgram::CommandType type, unsigned A    = 0, unsigned B    = 0, unsigned C    = 0, 
                                                      unsigned arg0 = 0, unsigned arg1 = 0, unsigned arg2 = 0);
    void compactify();
    void remove_noop();
    void calculate_variable_usage_intervals();
    void replace_output_var(unsigned old_id, unsigned new_id);
    void reset_alias_rec(unsigned alias_id, unsigned master_id, unsigned base_offset);
    void set_alias(unsigned alias_id, unsigned master_id, unsigned from, unsigned to);
    bool optimize_unused_cycle();
    void optimize_renaming_moves();
    void optimize_self_applicable_commands();
    void optimize_copy_to_aliases();
    void optimize_program();
    unsigned calculate_memory_layout();

    std::vector<Variable> vars;
    std::vector<TensorProgram::Command> commands;
    
    std::vector<float> constants;
    std::map<std::string, unsigned> input_vars;  //name -> var_id
    std::map<std::string, unsigned> output_vars; //name -> var_id

    std::vector<unsigned> fm, lm, fu, lu; //variable usage intervals
  };

  struct TensorToken
  {
    static TensorCompiler *tp;

    //constructors and assignments
    TensorToken();
    TensorToken(float val);
    TensorToken(const std::vector<unsigned> &shape);
    TensorToken(const unsigned _sizes[TensorCompiler::MAX_DIM]);
    explicit TensorToken(int sz_0, int sz_1 = 0, int sz_2 = 0, int sz_3 = 0);
    explicit TensorToken(unsigned sz_0, unsigned sz_1 = 0, unsigned sz_2 = 0, unsigned sz_3 = 0);
    TensorToken(const TensorToken &other);
    TensorToken &operator=(const TensorToken &other);

    //helper functions
    unsigned total_size() const;
    void check_dimensions_for_arithmetics(const TensorToken &other) const;
    TensorToken reshape(std::vector<unsigned> new_shape) const;
    TensorToken flatten() const;

    //data manipulation
    TensorToken get(unsigned n) const;
    TensorToken get(std::pair<unsigned, unsigned> range) const;
    void set(unsigned n, const TensorToken &t);
    void set_flatten(std::pair<unsigned, unsigned> range, const TensorToken &t);
    void copy_to(std::pair<unsigned, unsigned> to_range, const TensorToken &t, std::pair<unsigned, unsigned> from_range);
    void set(std::pair<unsigned, unsigned> range, const TensorToken &t);
    void fill(float val);

    //ariphmetic
    TensorToken &operator+=(const TensorToken &other);
    TensorToken operator+(const TensorToken &other) const;
    TensorToken &operator*=(const TensorToken &other);
    TensorToken operator*(const TensorToken &other) const;
    TensorToken &operator-=(const TensorToken &other);
    TensorToken operator-(const TensorToken &other) const;
    TensorToken &operator/=(const TensorToken &other);
    TensorToken operator/(const TensorToken &other) const;

    //per-element operations
    static TensorToken exp(const TensorToken &A);
    static TensorToken pow(const TensorToken &A, const TensorToken &B);
    static TensorToken sin(const TensorToken &A);
    static TensorToken cos(const TensorToken &A);
    static TensorToken log(const TensorToken &A);

    //aggregation
    TensorToken sum(int Dims = -1) const;
    TensorToken outer_sum() const;

    //linear algebra operations
    TensorToken transpose() const;
    static TensorToken vector_outer_product(const TensorToken &A, const TensorToken &B);
    static TensorToken vector_outer_product_sum(const TensorToken &A, const TensorToken &B);
    static TensorToken mat_mul_t(const TensorToken &A, const TensorToken &B);
    static TensorToken mat_vec_mul(const TensorToken &A, const TensorToken &B);

    unsigned id = 0;
    unsigned Dim = 0;
    unsigned sizes[TensorCompiler::MAX_DIM] = {0};
  };
}