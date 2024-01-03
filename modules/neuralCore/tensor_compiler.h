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
    TensorProgram finish_program();
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
    
    std::map<std::string, unsigned> input_vars;  //name -> var_id
    std::map<std::string, unsigned> output_vars; //name -> var_id

    std::vector<unsigned> fm, lm, fu, lu; //variable usage intervals
  };

  struct TensorToken
  {
    static TensorCompiler *tp;

    TensorToken()
    {
      Dim = 0;
      id = tp->add_var(*this);
    }
    TensorToken(float val)
    {
      Dim = 0;
      sizes[0] = 1;
      id = tp->add_var(*this);
      tp->ftt(id, val);
    }
    TensorToken(const std::vector<unsigned> &shape)
    {
      assert(shape.size() <= TensorCompiler::MAX_DIM);
      Dim = shape.size();
      for (int i=0;i<Dim;i++)
        sizes[i] = shape[i];
      id = tp->add_var(*this);
    }
    TensorToken(const unsigned _sizes[TensorCompiler::MAX_DIM]) : TensorToken(_sizes[0], _sizes[1], _sizes[2], _sizes[3])
    {
    }
    explicit TensorToken(int sz_0, int sz_1 = 0, int sz_2 = 0, int sz_3 = 0) : TensorToken((unsigned)sz_0, (unsigned)sz_1, (unsigned)sz_2, (unsigned)sz_3) {}
    explicit TensorToken(unsigned sz_0, unsigned sz_1 = 0, unsigned sz_2 = 0, unsigned sz_3 = 0)
    {
      static_assert(TensorCompiler::MAX_DIM == 4);
      Dim = sz_0 == 0 ? 0 : (sz_1 == 0 ? 1 : (sz_2 == 0 ? 2 : (sz_3 == 0 ? 3 : 4)));
      sizes[0] = sz_0;
      sizes[1] = sz_1;
      sizes[2] = sz_2;
      sizes[3] = sz_3;
      id = tp->add_var(*this);
      printf("created tensor %u - %u %u %u %u\n",id, sz_0, sz_1, sz_2, sz_3);
    }
    TensorToken(const TensorToken &other)
    {
      Dim = other.Dim;
      for (int i = 0; i < Dim; i++)
        sizes[i] = other.sizes[i];
      id = tp->add_var(*this);
      tp->add_command(TensorProgram::MOV, other.id, 0, id);
    }
    TensorToken &operator=(const TensorToken &other)
    {
      Dim = other.Dim;
      for (int i = 0; i < Dim; i++)
        sizes[i] = other.sizes[i];
      id = tp->add_var(*this);
      tp->add_command(TensorProgram::MOV, other.id, 0, id);
      return *this;
    }

    unsigned total_size() const
    {
      unsigned size = 1;
      for (int i = 0; i < Dim; i++)
        size *= sizes[i];
      return size;
    }

    TensorToken &operator+=(const TensorToken &other) 
    {
      tp->add_command(TensorProgram::ADD, id, other.id, id);
      return *this;
    }

    TensorToken operator+(const TensorToken &other) const
    {
      TensorToken res(sizes);
      tp->add_command(TensorProgram::ADD, id, other.id, res.id);
      return res;
    }

    TensorToken &operator/=(const TensorToken &other)
    {
      tp->add_command(TensorProgram::DIV, id, other.id, id);
      return *this;
    }

    TensorToken operator/(const TensorToken &other) const
    {
      TensorToken res(sizes);
      tp->add_command(TensorProgram::DIV, id, other.id, res.id);
      return res;
    }

    TensorToken exp() const
    {
      TensorToken res(sizes);
      tp->add_command(TensorProgram::EXP, id, 0, res.id);
      return res;
    }

    TensorToken sum(int Dims = -1) const
    {
      if (Dim == 0) // sum of scalar is this scalar itself
        return *this;
      if (Dims == -1)
        Dims = Dim;
      assert(Dims > 0);
      assert(Dims <= Dim);
      unsigned res_Dim = Dim - Dims; // remaining dimensions
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      for (int i = 0; i < res_Dim; i++)
        res_sizes[i] = sizes[i + Dims];

      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::SUM, id, 0, res.id);
      return res;
    }

    TensorToken transpose() const
    {
      assert(Dim > 1);
      unsigned res_Dim = Dim;
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      res_sizes[0] = sizes[1];
      res_sizes[1] = sizes[0];
      for (int i = 2; i < res_Dim; i++)
        res_sizes[i] = sizes[i];

      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::TRANSP, id, 0, res.id);
      return res;
    }

    TensorToken get(unsigned n) const
    {
      assert(Dim > 0);
      assert(n < sizes[Dim-1]);

      unsigned res_size = 1;
      unsigned res_Dim = Dim-1;
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      for (int i = 0; i < res_Dim; i++)
      {
        res_sizes[i] = sizes[i];
        res_size *= res_sizes[i];
      }
      TensorToken res(res_sizes);

      tp->add_command(TensorProgram::COPY, id, 0, res.id, n*res_size, 0, res_size);
      return res;
    }

    TensorToken get(std::pair<unsigned, unsigned> range) const
    {
      unsigned from = range.first;
      unsigned to = range.second;
      assert(Dim > 0);
      assert(from < to);
      assert(to <= sizes[Dim-1]);

      unsigned res_size = 1;
      unsigned res_Dim = Dim;
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      for (int i = 0; i < Dim-1; i++)
      {
        res_sizes[i] = sizes[i];
        res_size *= res_sizes[i];
      }
      res_sizes[Dim-1] = to-from;
      TensorToken res(res_sizes);

      tp->add_command(TensorProgram::COPY, id, 0, res.id, from*res_size, 0, (to-from)*res_size);
      return res;      
    }

    void set(unsigned n, const TensorToken &t)
    {
      assert(t.Dim == Dim-1);
      assert(n <= sizes[Dim-1]);
      unsigned sub_size = 1;
      for (int i = 0; i < Dim-1; i++)
      {
        assert(sizes[i] == t.sizes[i]);
        sub_size *= sizes[i];
      }
      tp->add_command(TensorProgram::COPY, t.id, 0, id, 0, n*sub_size, sub_size);
    }

    void set_flatten(std::pair<unsigned, unsigned> range, const TensorToken &t)
    {
      unsigned from = range.first;
      unsigned to = range.second;
      assert(t.total_size() == to-from);
      assert(from < to);
      assert(to <= total_size());

      tp->add_command(TensorProgram::COPY, t.id, 0, id, from, 0, (to-from));
    }

    void copy_to(std::pair<unsigned, unsigned> to_range, const TensorToken &t, std::pair<unsigned, unsigned> from_range)
    {
      assert(to_range.second <= total_size());
      assert(from_range.second <= t.total_size());
      assert(to_range.second-to_range.first == from_range.second-from_range.first);
      tp->add_command(TensorProgram::COPY, t.id, 0, id, from_range.first, to_range.first, to_range.second-to_range.first);
    }

    void set(std::pair<unsigned, unsigned> range, const TensorToken &t)
    {
      unsigned from = range.first;
      unsigned to = range.second;
      assert(t.Dim == Dim);
      assert(from < to);
      assert(to <= sizes[Dim-1]);
      assert(to-from <= t.sizes[Dim-1]);
      unsigned sub_size = 1;
      for (int i = 0; i < Dim-1; i++)
      {
        assert(sizes[i] == t.sizes[i]);
        sub_size *= sizes[i];
      }
      tp->add_command(TensorProgram::COPY, t.id, 0, id, 0, from*sub_size, (to-from)*sub_size);
    }

    TensorToken reshape(std::vector<unsigned> new_shape) const
    {
      unsigned new_size = 1;
      for (unsigned s : new_shape)
        new_size *= s;
      unsigned size = 1;
      for (int i = 0; i < Dim; i++)
        size *= sizes[i];

      assert(new_size == size);
      assert(new_shape.size() <= TensorCompiler::MAX_DIM);

      unsigned res_Dim = new_shape.size();
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      for (int i = 0; i < new_shape.size(); i++)
        res_sizes[i] = new_shape[i];

      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::MOV, id, 0, res.id);
      return res;    
    }
    
    TensorToken flatten() const
    {
      unsigned size = 1;
      for (int i = 0; i < Dim; i++)
        size *= sizes[i];
      return reshape({size});
    }

    static TensorToken vector_outer_product(const TensorToken &A, const TensorToken &B)
    {
      assert(A.Dim >= 1);
      assert(A.Dim < TensorCompiler::MAX_DIM);
      assert(B.Dim == 1);
      assert(A.sizes[0] == B.sizes[0]);
      
      unsigned res_Dim = A.Dim + 1;
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      res_sizes[0] = B.sizes[0];
      res_sizes[1] = A.sizes[0];
      for (int i = 2; i < res_Dim; i++)
        res_sizes[i] = A.sizes[i - 1];

      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::OUTER_P, A.id, B.id, res.id);
      return res;
    }

    static TensorToken mat_mul_t(const TensorToken &A, const TensorToken &B)
    {
      assert(A.Dim == 2);
      assert(B.Dim == 2);
      assert(A.sizes[0] == B.sizes[0]);

      unsigned res_Dim = 2;
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      res_sizes[0] = A.sizes[1];
      res_sizes[1] = B.sizes[1];
      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::MATMUL_T, A.id, B.id, res.id);
      return res;
    }

    static TensorToken mat_vec_mul(const TensorToken &A, const TensorToken &B)
    {
      assert(A.Dim == 2);
      assert(B.Dim == 1);
      assert(A.sizes[0] == B.sizes[0]);

      unsigned res_Dim = 1;
      unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
      res_sizes[0] = A.sizes[1];
      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::MATMUL_T, A.id, B.id, res.id);
      return res;
    }

    unsigned id;
    unsigned Dim;
    unsigned sizes[TensorCompiler::MAX_DIM] = {0};
  };
}