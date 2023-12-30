#pragma once
#include <vector>
#include <cmath>

class TensorProcessor
{
public:
  constexpr static unsigned ADD = 0;
  constexpr static unsigned DIV = 1;
  constexpr static unsigned EXP = 2;
  constexpr static unsigned SUM = 3;
  constexpr static unsigned MATMUL_T = 4;
  struct Command
  {
    unsigned type;
    unsigned args[4];
  };

  struct Variable
  {
    unsigned Dim;
    unsigned offset;
    unsigned total_size;
    unsigned sizes[4];
  };
  TensorProcessor(){};
  virtual void process(const Command *commands __attribute__((size("cmd_count"))), unsigned cmd_count,
                       const Variable *vars __attribute__((size("var_count"))), unsigned var_count,
                       float *memory __attribute__((size("data_size"))), unsigned data_size);

  virtual void CommitDeviceData() {}                                                         // will be overriden in generated class
  virtual void GetExecutionTime(const char *a_funcName, float a_out[4]) { a_out[0] = 0.0f; } // will be overriden in generated class
protected:
  virtual void kernel2D_add(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C); // C = A + B
  virtual void kernel2D_div(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C); // C = A / B
  virtual void kernel1D_exp(float *data, unsigned steps, Variable A, Variable B);                                 // B = exp(A)
  virtual void kernel2D_sum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B);             // B = sum(A)
  virtual void matmul_transposed(float *data, Variable A, Variable B, Variable C);                                // C = A * (B)^T
};