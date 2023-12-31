#pragma once
#include "tensor_processor.h"

class TensorProcessorImpl
{
protected:
  friend class nn::TensorProcessor;
  using Command = nn::TensorProgram::Command;
  using Variable = nn::TensorProgram::Variable;
  TensorProcessorImpl(){};
  virtual void process(const nn::TensorProgram &program,
                       const float *memory_in __attribute__((size("data_size"))),
                       float *memory_out __attribute__((size("data_size"))), unsigned data_size);

  virtual void CommitDeviceData() {}                                                         // will be overriden in generated class
  virtual void GetExecutionTime(const char *a_funcName, float a_out[4]) { a_out[0] = 0.0f; } // will be overriden in generated class
  virtual void kernel2D_add(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C); // C = A + B
  virtual void kernel2D_div(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C); // C = A / B
  virtual void kernel1D_exp(float *data, unsigned steps, Variable A, Variable B);                                 // B = exp(A)
  virtual void kernel1D_sum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B);             // B = sum(A)
  virtual void matmul_transposed(float *data, Variable A, Variable B, Variable C);                                // C = A * (B)^T
};