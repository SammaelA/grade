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
    constexpr static unsigned MAX_DIM = 8;
    constexpr static unsigned CMD_ARGS = 8;
    constexpr static unsigned CONSTS_VAR_ID = 1;

    //each command is [type, A, B, C, arg0, arg1, arg2]
    //A,B - input variables
    //C - output variable
    enum CommandType
    {
      NOOP,     // do nothing

      MOV,      // memcpy(C,A, sizeof(float)*A.total_size)
      FILL,     // fill(C, as_float(arg0))
      COPY,     // memcpy(C+arg1, A+arg0, sizeof(float)*arg2)
      PAD,      // padding along the given axis
      FLIP,     // reverse order of values along the given axis
      DILATE,   // put some zero values between values from input tensor
      URAND,    // fills tensor with values, uniformly distributed in [0,1]

      ADD,      // C = A+B
      SUB,      // C = A-B
      MUL,      // C = A*B
      DIV,      // C = A/B
      GREATER,  // C = A > B
      LESS,     // C = A < B
      EQUAL,    // C = A == B (precisely)
      GE,       // C = A >= B
      LE,       // C = A <= B
      NE,       // C = A != B (precisely)
      OR,       // C = A>0 || B>0
      AND,      // C = A>0 && B>0
      WHERE,    // C = B>0 ? A : 0 (elementwise)
      MIN,      // C = min(A, B)
      MAX,      // C = min(A, B)
      POW,      // C = A^B

      EXP,      // C = exp(A)
      SQRT,     // C = sqrt(A)
      SIN,      // C = sin(A)
      COS,      // C = cos(A)
      LOG,      // C = log(A)
      NOT,      // C = A > 0 ? 0 : 1 (elementwise)

      SUM,      // C = sum(A)
      O_SUM,    // C = sum(A)
      MINIMUM,  // C = min(A)
      MAXIMUM,  // C = max(A)

      MATMUL_T, // C = Ax(B^T)
      TRANSP,   // C = transpose(A)
      OUTER_P,  // C = outer_product(A, B)
      SMAX_D,   // derivative of softmax function. It's complicated enough to have a separate command for it
      CONV_2D,  //convolution with arbitrary number of channels and filters. Borders are ignored
      MPOOL,    // C = max pooling(A) with arbitrary window size
      MPOOL_D,  // derivative of max pooling
      CONV_3D,  //convolution with arbitrary number of channels and filters. Borders are ignored
      MPOOL_3D, // C = 3D max pooling(A) with arbitrary window size
      MPOOL_3D_D,// derivative of 3D max pooling

      
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
      unsigned args[8]; //CMD_ARGS
    };

    struct Variable
    {
      unsigned Dim;
      unsigned offset;
      unsigned total_size;
      unsigned sizes[8]; //MAX_DIM
    };

    static std::vector<CmdProperties> cmd_properties;

    std::vector<Command> commands;
    std::vector<Variable> vars;
    std::vector<float> constants;

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
    static void set_input(const std::string &name, const float * const data, unsigned data_size);
    //transfers data from output tensor with <name> to given address
    //if <data_size> less that tensor size, only this part is tranfered
    static void get_output(const std::string &name, float *data, unsigned data_size);
    static void execute();
    static void print_execution_stat();
  private:
    TensorProcessor();
    std::shared_ptr<TensorProcessorImpl> pImpl;
    TensorProgram program;
    std::map<std::string, bool> input_prepared;
    bool program_prepared = false;
    std::string used_backend = "CPU";
  };
}