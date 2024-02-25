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

    void start_program();
    TensorProgram finish_program(bool print_program = false);
    void input(const TensorToken &t, std::string name);
    void output(const TensorToken &t, std::string name);
    void inout(const TensorToken &t, std::string name);

  private:
    struct Variable
    {
      unsigned Dim;
      unsigned offset;
      unsigned total_size;
      unsigned sizes[TensorProgram::MAX_DIM];

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
                                                      unsigned arg0 = 0, unsigned arg1 = 0, unsigned arg2 = 0,
                                                      unsigned arg3 = 0, unsigned arg4 = 0);
    void compactify();
    void remove_noop();
    void calculate_variable_usage_intervals();
    void calculate_variable_usage_interval_with_aliases(unsigned v_id);
    void replace_output_var(unsigned old_id, unsigned new_id);
    void reset_alias_rec(unsigned alias_id, unsigned master_id, unsigned base_offset);
    void set_alias(unsigned alias_id, unsigned master_id, unsigned from, unsigned to);
    bool optimize_unused_cycle();
    bool optimize_renaming_moves();
    bool optimize_self_applicable_commands();
    void optimize_copy_to_aliases();
    void optimize_program();
    unsigned calculate_memory_layout_naive();
    unsigned calculate_memory_layout_interval_coloring();

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
    TensorToken(const unsigned _sizes[TensorProgram::MAX_DIM]);
    explicit TensorToken(int sz_0, int sz_1 = 0, int sz_2 = 0, int sz_3 = 0,
                         int sz_4 = 0, int sz_5 = 0, int sz_6 = 0, int sz_7 = 0);
    explicit TensorToken(unsigned sz_0, unsigned sz_1 = 0, unsigned sz_2 = 0, unsigned sz_3 = 0,
                         unsigned sz_4 = 0, unsigned sz_5 = 0, unsigned sz_6 = 0, unsigned sz_7 = 0);
    TensorToken(const TensorToken &other);
    TensorToken &operator=(const TensorToken &other);

    //helper functions
    unsigned total_size() const;
    void check_dimensions_for_arithmetics(const TensorToken &other) const;
    TensorToken reshape(std::vector<unsigned> new_shape) const;
    TensorToken flatten() const;
    static void issue_command(TensorProgram::CommandType type, const TensorToken &A, const TensorToken &B, const TensorToken &C, 
                              unsigned arg0 = 0, unsigned arg1 = 0, unsigned arg2 = 0, unsigned arg3 = 0, unsigned arg4 = 0);

    //data manipulation
    TensorToken get(unsigned n) const;
    TensorToken get(std::pair<unsigned, unsigned> range) const;
    void set(unsigned n, const TensorToken &t);
    void set_flatten(std::pair<unsigned, unsigned> range, const TensorToken &t);
    void copy_to(std::pair<unsigned, unsigned> to_range, const TensorToken &t, std::pair<unsigned, unsigned> from_range);
    void set(std::pair<unsigned, unsigned> range, const TensorToken &t);
    void fill(float val);
    TensorToken add_padding(unsigned left_pad, unsigned right_pad, int Dim = 0) const;
    TensorToken flip(unsigned axis) const;
    void random(unsigned seed = 0u); // fills tensor with values, uniformly distributed in [0,1]

    static void g_2op(TensorProgram::CommandType cmd, const TensorToken &A, const TensorToken &B, const TensorToken &C,
                      unsigned steps, unsigned step_size, unsigned group_size);
    static TensorToken g_2op(TensorProgram::CommandType cmd, const TensorToken &A, const TensorToken &B,
                             unsigned steps, unsigned step_size, unsigned group_size);
    static TensorToken g_2op(TensorProgram::CommandType cmd, const TensorToken &A, const TensorToken &B,
                             unsigned start_dim = 0);
    
    //arithmetics
    TensorToken &operator+=(const TensorToken &other);
    TensorToken operator+(const TensorToken &other) const;
    friend TensorToken operator+(float A, const TensorToken &B) { return TensorToken(A) + B; }
    TensorToken &operator*=(const TensorToken &other);
    TensorToken operator*(const TensorToken &other) const;
    friend TensorToken operator*(float A, const TensorToken &B) { return TensorToken(A) * B; }
    TensorToken &operator-=(const TensorToken &other);
    TensorToken operator-(const TensorToken &other) const;
    friend TensorToken operator-(float A, const TensorToken &B) { return TensorToken(A) - B; }
    TensorToken &operator/=(const TensorToken &other);
    TensorToken operator/(const TensorToken &other) const;
    friend TensorToken operator/(float A, const TensorToken &B) { return TensorToken(A) / B; }

    static TensorToken pow(const TensorToken &A, const TensorToken &B);

    //per-element operations
    static TensorToken exp(const TensorToken &A);
    static TensorToken sin(const TensorToken &A);
    static TensorToken cos(const TensorToken &A);
    static TensorToken log(const TensorToken &A);
    static TensorToken sqrt(const TensorToken &A);

    //aggregation
    TensorToken sum(int Dims = -1) const;
    TensorToken outer_sum() const;
    TensorToken minimum(int Dims = -1) const;
    TensorToken maximum(int Dims = -1) const;

    //linear algebra operations
    TensorToken transpose(unsigned transp_dim = 0) const; //swaps (transp_dim) and (transp_dim+1)
    static TensorToken vector_outer_product(const TensorToken &A, const TensorToken &B);
    static TensorToken mat_mul_t(const TensorToken &A, const TensorToken &B);
    static TensorToken conv2D(const TensorToken &A, const TensorToken &kernel, unsigned stride = 1);
    static TensorToken conv3D(const TensorToken &A, const TensorToken &kernel, unsigned stride = 1);

    unsigned id = 0;
    unsigned Dim = 0;
    unsigned sizes[TensorProgram::MAX_DIM] = {0};
  };
}