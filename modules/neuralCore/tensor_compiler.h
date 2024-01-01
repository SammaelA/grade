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
    };

    unsigned add_var(const TensorToken &t);
    void ftt(unsigned id, float val);
    void add_command(TensorProgram::CommandType type, unsigned A = 0, unsigned B = 0, unsigned C = 0, unsigned num_arg = 0);
    void optimize_program();
    unsigned calculate_memory_layout();

    std::vector<Variable> vars;
    std::vector<TensorProgram::Command> commands;
    std::vector<float> constants;
    
    std::map<std::string, unsigned> input_vars;  //name -> var_id
    std::map<std::string, unsigned> output_vars; //name -> var_id
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
      Dim = 1;
      sizes[0] = 1;
      id = tp->add_var(*this);
      tp->ftt(id, val);
    }
    TensorToken(unsigned _sizes[TensorCompiler::MAX_DIM]) : TensorToken(_sizes[0], _sizes[1], _sizes[2], _sizes[3])
    {
    }
    explicit TensorToken(int sz_0, int sz_1 = 0, int sz_2 = 0, int sz_3 = 0) : TensorToken((unsigned)sz_0, (unsigned)sz_1, (unsigned)sz_2, (unsigned)sz_3) {}
    explicit TensorToken(unsigned sz_0, unsigned sz_1 = 0, unsigned sz_2 = 0, unsigned sz_3 = 0)
    {
      assert(TensorCompiler::MAX_DIM == 4);
      Dim = sz_0 == 0 ? 0 : (sz_1 == 0 ? 1 : (sz_1 == 0 ? 2 : (sz_1 == 0 ? 3 : 4)));
      sizes[0] = sz_0;
      sizes[1] = sz_1;
      sizes[2] = sz_2;
      sizes[3] = sz_3;
      id = tp->add_var(*this);
    }
    TensorToken(const TensorToken &other)
    {
      Dim = other.Dim;
      for (int i = 0; i < Dim; i++)
        sizes[i] = other.sizes[i];
      id = tp->add_var(*this);
      tp->add_command(TensorProgram::MOV, other.id, id);
    }
    TensorToken &operator=(const TensorToken &other)
    {
      Dim = other.Dim;
      for (int i = 0; i < Dim; i++)
        sizes[i] = other.sizes[i];
      tp->add_command(TensorProgram::MOV, other.id, id);
      return *this;
    }

    TensorToken &operator+=(const TensorToken &other)
    {
      tp->add_command(TensorProgram::ADD, id, other.id, id);
      return *this;
    }

    TensorToken operator+(const TensorToken &other)
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

    TensorToken operator/(const TensorToken &other)
    {
      TensorToken res(sizes);
      tp->add_command(TensorProgram::DIV, id, other.id, res.id);
      return res;
    }

    TensorToken exp()
    {
      TensorToken res(sizes);
      tp->add_command(TensorProgram::EXP, id, res.id);
      return res;
    }

    TensorToken sum(int Dims = -1)
    {
      if (Dim == 0) // sum of scalar is this scalar itself
        return *this;
      if (Dims == -1)
        Dims = Dim;
      assert(Dims > 0);
      assert(Dims <= Dim);
      unsigned res_Dim = Dim - Dims; // remaining dimensions
      unsigned res_sizes[TensorCompiler::MAX_DIM];
      for (int i = 0; i < res_Dim; i++)
        res_sizes[i] = sizes[i + Dims];

      TensorToken res(res_sizes);
      tp->add_command(TensorProgram::SUM, id, res.id);
      return res;
    }

    unsigned id;
    unsigned Dim;
    unsigned sizes[TensorCompiler::MAX_DIM];
  };
}