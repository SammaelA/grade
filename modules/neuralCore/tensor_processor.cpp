#include "tensor_processor.h"

#define GET(A, i) (data[(A).offset + i])

void TensorProcessor::process(const Command *commands, unsigned cmd_count,
                              const Variable *vars, unsigned var_count,
                              float *memory, unsigned data_size)
{
  for (int i = 0; i < cmd_count; i++)
  {
    Variable A = vars[commands[i].args[0]];
    Variable B = vars[commands[i].args[1]];
    Variable C = vars[commands[i].args[2]];
    Variable D = vars[commands[i].args[3]];
    switch (commands[i].type)
    {
    case ADD:
      kernel2D_add(memory, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case DIV:
      kernel2D_div(memory, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case EXP:
      kernel1D_exp(memory, A.total_size, A, B);
      break;
    case SUM:
      kernel2D_sum(memory, B.total_size, A.total_size / B.total_size, A, B);
      break;
    case MATMUL_T:
      matmul_transposed(memory, A, B, C);
      break;
    default:
      break;
    }
  }
}

void TensorProcessor::kernel2D_add(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A + B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] + data[B.offset + j];
}
void TensorProcessor::kernel2D_div(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A / B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] / data[B.offset + j];
}
void TensorProcessor::kernel1D_exp(float *data, unsigned steps, Variable A, Variable B) // B = exp(A)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = std::exp(data[A.offset + i]);
}
void TensorProcessor::kernel2D_sum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = 0;
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] += data[A.offset + i * step_size + j];
  }
}
void TensorProcessor::matmul_transposed(float *data, Variable A, Variable B, Variable C) // C = A * (B)^T
{
}