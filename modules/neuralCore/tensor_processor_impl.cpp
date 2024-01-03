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

    unsigned arg0 = program.commands[i].args[3];
    unsigned arg1 = program.commands[i].args[4];
    unsigned arg2 = program.commands[i].args[5];

    switch (program.commands[i].type)
    {
    case nn::TensorProgram::NOOP:
      break;
    case nn::TensorProgram::ADD:
      kernel2D_add(memory_out, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::MUL:
      kernel2D_mul(memory_out, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::DIV:
      kernel2D_div(memory_out, A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::EXP:
      kernel1D_exp(memory_out, A.total_size, A, C);
      break;
    case nn::TensorProgram::SUM:
      kernel1D_sum(memory_out, C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::MATMUL_T:
      kernel2D_matmul_transposed(memory_out, A.sizes[0], A.sizes[1], std::max(1u, C.sizes[1]), A, B, C);
      break;
    case nn::TensorProgram::MOV:
      memcpy(memory_out + C.offset, memory_out + A.offset, sizeof(float)*A.total_size);
      break;
    case nn::TensorProgram::FTT:
      kernel1D_fill(memory_out, C.total_size, C, *((float *)(&arg0)));
      break;
    case nn::TensorProgram::COPY:
      memcpy(memory_out + C.offset + arg1, memory_out + A.offset + arg0, sizeof(float)*arg2);
      break;
    case nn::TensorProgram::TRANSP:
      kernel2D_transpose(memory_out, A.total_size/(A.sizes[0]*A.sizes[1]), A.sizes[0], A.sizes[1], A, C);
      break;
    case nn::TensorProgram::OUTER_P:
      kernel2D_outer_product(memory_out, A.total_size/A.sizes[0], A.sizes[0], B.sizes[0], A, B, C);
      break;
    case nn::TensorProgram::URAND:
    {
      float from = memory_out[A.offset];
      float to = memory_out[B.offset];
      //TODO: GPU-compatible random
      for (unsigned i = 0; i < C.total_size; i++)
        memory_out[C.offset + i] = from + ((double)rand()/RAND_MAX)*(to-from);
    }
      break;
    default:
      break;
    }
  }
}

void TensorProcessorImpl::kernel1D_fill(float *data, unsigned steps, Variable A, float val)
{
  for (unsigned i = 0; i < steps; i++) 
    data[A.offset + i] = val;
}

void TensorProcessorImpl::kernel2D_add(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A + B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] + data[B.offset + j];
}
void TensorProcessorImpl::kernel2D_mul(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A * B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] * data[B.offset + j];
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

void TensorProcessorImpl::kernel2D_transpose(float *data, unsigned steps, unsigned row_len, unsigned col_len, Variable A, Variable B)
{
  for (unsigned s = 0; s < steps; s++)
    for (unsigned i = 0; i < row_len; i++)
      for (unsigned j = 0; j < col_len; j++)
        data[B.offset + s*col_len*row_len + i*col_len + j] = data[A.offset + s*col_len*row_len + j*row_len + i];
}

void TensorProcessorImpl::kernel2D_matmul_transposed(float *data, unsigned A_row_len, unsigned A_col_len, unsigned B_col_len, 
                                          Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < B_col_len; i++)
  {
    for (unsigned j = 0; j < A_col_len; j++)
    {
      data[C.offset + i*A_col_len + j] = 0;
      for (unsigned k = 0; k < A_row_len; k++)
        data[C.offset + i*A_col_len + j] += data[A.offset + j*A_row_len + k]*data[B.offset + i*A_row_len + k];
    }
  }
}

void TensorProcessorImpl::kernel2D_outer_product(float *data, unsigned steps, unsigned A_len, unsigned B_len, 
                                      Variable A, Variable B, Variable C)
{
  for (unsigned s = 0; s < steps; s++)
    for (unsigned i = 0; i < A_len; i++)
      for (unsigned j = 0; j < B_len; j++)
        data[C.offset + s*A_len*B_len + i*B_len + j] = data[A.offset + s*A_len + i]*data[B.offset + j];
}
