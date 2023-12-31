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
    enum CommandType
    {
      ADD,
      DIV,
      EXP,
      SUM,
      MATMUL_T,
      MOV,
      FTT
    };

    struct Command
    {
      CommandType type;
      unsigned args[4];
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
    void set_input(const std::map<std::string, float * const> &vars);
    void set_output(const std::map<std::string, float *> &vars);
    void execute();
  private:
    std::shared_ptr<TensorProcessorImpl> pImpl;
    TensorProgram program;
    std::vector<float> cpu_memory;
    bool input_prepared = false;
    bool program_prepared = false;
  };
}