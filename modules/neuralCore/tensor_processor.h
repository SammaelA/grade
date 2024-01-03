#pragma once
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <memory>

class TensorProcessorImpl;
namespace nn
{
  struct TensorProgram
  {
    //each command is [type, A, B, C, arg0, arg1, arg2]
    //A,B - input variables
    //C - output variable
    enum CommandType
    {
      NOOP,
      ADD,      // C = A+B
      MUL,      // C = A*B
      DIV,      // C = A/B
      EXP,      // C = exp(A)
      SUM,      // C = sum(A)
      MATMUL_T, // C = AxB
      MOV,      // memcpy(C,A, sizeof(float)*A.total_size)
      FTT,      // C = as_float(arg0)
      COPY,     // memcpy(C+arg1, A+arg0, sizeof(float)*arg2)
      TRANSP,   // C = transpose(A)
      OUTER_P,  // C = outer_product(A, B)
      URAND,    // C = random_float_uniform(A, B)
    };

    struct Command
    {
      CommandType type;
      unsigned args[6];
    };

    struct Variable
    {
      unsigned Dim;
      unsigned offset;
      unsigned total_size;
      unsigned sizes[4];
    };

    std::vector<Command> commands;
    std::vector<Variable> vars;

    unsigned total_memory_req;

    std::map<std::string, unsigned> input_vars;  //name -> var_id
    std::map<std::string, unsigned> output_vars; //name -> var_id
  };

  class TensorProcessor
  {
  public:
    TensorProcessor();
    void set_program(const TensorProgram &program);
    void set_input(const std::string &name, float * const data);
    void set_input(const std::map<std::string, float * const> &vars);
    void get_output(const std::map<std::string, float *> &vars);
    void get_output(const std::string &name, float *data);
    void execute();
  private:
    std::shared_ptr<TensorProcessorImpl> pImpl;
    TensorProgram program;
    std::vector<float> cpu_memory;
    std::map<std::string, bool> input_prepared;
    bool program_prepared = false;
  };
}