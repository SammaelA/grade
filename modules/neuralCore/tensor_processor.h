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
      NOOP,     // do nothing
      FTT,      // C = as_float(arg0)

      MOV,      // memcpy(C,A, sizeof(float)*A.total_size)
      FILL,     // fill(C, as_float(arg0))
      COPY,     // memcpy(C+arg1, A+arg0, sizeof(float)*arg2)

      ADD,      // C = A+B
      SUB,      // C = A-B
      MUL,      // C = A*B
      DIV,      // C = A/B

      EXP,      // C = exp(A)
      POW,      // C = A^B (B is 1-element tensor)
      SIN,      // C = sin(A)
      COS,      // C = cos(A)
      LOG,      // C = log(A)

      SUM,      // C = sum(A)
      O_SUM,    // C = sum(A)

      MATMUL_T, // C = AxB
      TRANSP,   // C = transpose(A)
      OUTER_P,  // C = outer_product(A, B)
      OUTER_PS, // C = sum(outer_product(A, B)) *

      URAND,    // C = random_float_uniform(A, B)
      
      CMD_COUNT
    };

    enum CmdClass
    {
      AUXILIARY,
      MEM_MANAGEMENT,
      ARITHMETICS,
      ELEMENTWISE,
      REDUCTION,
      ALGEBRA,
      OTHER
    };

    enum CmdIsSelfApplicable
    {
      SELF_APPLICABLE_NO,
      SELF_APPLICABLE_YES
    };

    struct CmdProperties
    {
      CommandType type;
      std::string name;
      CmdClass cls;
      CmdIsSelfApplicable is_self_applicable;
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

    static std::vector<CmdProperties> cmd_properties;

    std::vector<Command> commands;
    std::vector<Variable> vars;

    unsigned total_memory_req;

    std::map<std::string, unsigned> input_vars;  //name -> var_id
    std::map<std::string, unsigned> output_vars; //name -> var_id
  };

  class TensorProcessor
  {
  public:
    static void init(std::string backend = "CPU");
    //sets given program for execution. Initializes memory etc.
    static void set_program(const TensorProgram &program);
    //transfers data to input tensor with <name>
    //if <data_size> less that tensor size, remaining part is padded with zeros
    //all inputs should be set before execution
    static void set_input(const std::string &name, float * const data, unsigned data_size);
    //transfers data from output tensor with <name> to given address
    //if <data_size> less that tensor size, only this part is tranfered
    static void get_output(const std::string &name, float *data, unsigned data_size);
    static void execute();
  private:
    TensorProcessor();
    std::shared_ptr<TensorProcessorImpl> pImpl;
    TensorProgram program;
    std::vector<float> cpu_memory;
    std::map<std::string, bool> input_prepared;
    bool program_prepared = false;
  };
}