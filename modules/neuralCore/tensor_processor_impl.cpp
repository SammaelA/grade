#include "tensor_processor_impl.h"
#include <cstring>

void TensorProcessorImpl::process(const nn::TensorProgram &program,
                                  const float *memory_in, float *memory_out, unsigned data_size)
{
  memcpy(memory_out, memory_in, sizeof(float) * data_size);
  for (int i = 0; i < program.commands.size(); i++)
  {
    Variable A = program.vars[program.commands[i].args[0]];
    Variable B = program.vars[program.commands[i].args[1]];
    Variable C = program.vars[program.commands[i].args[2]];
    Variable D = program.vars[program.commands[i].args[3]];
    switch (program.commands[i].type)
    {
    case nn::TensorProgram::ADD:
      kernel2D_add(memory_out, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::DIV:
      kernel2D_div(memory_out, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::EXP:
      kernel1D_exp(memory_out, A.total_size, A, B);
      break;
    case nn::TensorProgram::SUM:
      kernel1D_sum(memory_out, B.total_size, A.total_size / B.total_size, A, B);
      break;
    case nn::TensorProgram::MATMUL_T:
      matmul_transposed(memory_out, A, B, C);
      break;
    case nn::TensorProgram::MOV:
      memcpy(memory_out + B.offset, memory_out + A.offset, sizeof(float)*A.total_size);
      break;
    default:
      break;
    }
  }
}

void TensorProcessorImpl::kernel2D_add(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A + B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] + data[B.offset + j];
}
void TensorProcessorImpl::kernel2D_div(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A / B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] / data[B.offset + j];
}
void TensorProcessorImpl::kernel1D_exp(float *data, unsigned steps, Variable A, Variable B) // B = exp(A)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = std::exp(data[A.offset + i]);
}
void TensorProcessorImpl::kernel1D_sum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = 0;
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] += data[A.offset + i * step_size + j];
  }
}
void TensorProcessorImpl::matmul_transposed(float *data, Variable A, Variable B, Variable C) // C = A * (B)^T
{
}