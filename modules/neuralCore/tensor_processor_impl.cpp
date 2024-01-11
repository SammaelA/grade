#include "tensor_processor_impl.h"
#include <cstring>
#include <chrono>
#define DEBUG 0

std::vector<float> _stat_time_cmd_num;
std::vector<float> _stat_time_cmd_id;
int _stat_execution_times;

void TensorProcessorImpl::process(const nn::TensorProgram &program)
{
  unsigned data_size = memory.size();
  #if DEBUG
  {
    printf("data [");
    for (int i=0;i<data_size;i++)
      printf("%8d ", i);
    printf("]\n");
  }
  #endif

  for (int i = 0; i < program.commands.size(); i++)
  {
    auto t1 = std::chrono::high_resolution_clock::now();

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
      kernel2D_add(memory.data(), A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::MUL:
      kernel2D_mul(memory.data(), A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::SUB:
      kernel2D_sub(memory.data(), A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::DIV:
      kernel2D_div(memory.data(), A.total_size / B.total_size, B.total_size, A, B, C);
      break;
    case nn::TensorProgram::EXP:
      kernel1D_exp(memory.data(), A.total_size, A, C);
      break;
    case nn::TensorProgram::POW:
      kernel1D_pow(memory.data(), A.total_size, A, B, C);
      break;
    case nn::TensorProgram::SIN:
      kernel1D_sin(memory.data(), A.total_size, A, C);
      break;
    case nn::TensorProgram::COS:
      kernel1D_cos(memory.data(), A.total_size, A, C);
      break;
    case nn::TensorProgram::LOG:
      kernel1D_log(memory.data(), A.total_size, A, C);
      break;
    case nn::TensorProgram::SUM:
      kernel1D_sum(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::O_SUM:
      kernel1D_osum(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::MATMUL_T:
      kernel2D_matmul_transposed(memory.data(), A.sizes[0], A.sizes[1], std::max(1u, C.sizes[1]), A, B, C);
      break;
    case nn::TensorProgram::MOV:
      kernel1D_copy(memory.data(), A.total_size, 0, 0, A, C);
      break;
    case nn::TensorProgram::FILL:
      kernel1D_fill(memory.data(), C.total_size, C, *((float *)(&arg0)));
      break;
    case nn::TensorProgram::COPY:
      kernel1D_copy(memory.data(), arg2, arg0, arg1, A, C);
      break;
    case nn::TensorProgram::TRANSP:
      kernel2D_transpose(memory.data(), A.total_size/(A.sizes[0]*A.sizes[1]), A.sizes[0], A.sizes[1], A, C);
      break;
    case nn::TensorProgram::OUTER_P:
      kernel2D_outer_product(memory.data(), A.total_size/A.sizes[0], A.sizes[0], B.sizes[0], A, B, C);
      break;
    case nn::TensorProgram::OUTER_PS:
    {
      kernel1D_fill(memory.data(), A.sizes[0]*B.sizes[0], C, 0.0f);
      kernel2D_outer_p_add(memory.data(), A.total_size/A.sizes[0], A.sizes[0], B.sizes[0], A, B, C);
      //for (unsigned s = 0; s < A.total_size/A.sizes[0]; s++)
      //{
      //  kernel2D_outer_p_add(memory.data(), s, A.sizes[0], B.sizes[0], A, B, C);
      //}
    }
      break;
    case nn::TensorProgram::URAND:
    {
      float from = memory.data()[A.offset];
      float to = memory.data()[B.offset];
      //TODO: GPU-compatible random
      for (unsigned i = 0; i < C.total_size; i++)
        memory.data()[C.offset + i] = from + ((double)rand()/RAND_MAX)*(to-from);
    }
      break;
    default:
      break;
    }
    #if DEBUG
    {
      printf("data [");
      for (int i=0;i<data_size;i++)
        printf("%8.4f ", memory.data()[i]);
      printf("]\n");
    }
    #endif

    auto t2 = std::chrono::high_resolution_clock::now();
    float total_time_us = 1e-3 * std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    _stat_time_cmd_id[program.commands[i].type] += total_time_us;
    _stat_time_cmd_num[i] += total_time_us;
  }
  _stat_execution_times++;
}

void TensorProcessorImpl::allocate_memory(unsigned size)
{
  memory.resize(size);
}

void TensorProcessorImpl::set_input(const float* in, unsigned offset, unsigned size)
{
  kernel1D_set_input(in, offset, size);
}

void TensorProcessorImpl::get_output(float* out, unsigned offset, unsigned size)
{
  kernel1D_get_output(out, offset, size);
}


void TensorProcessorImpl::kernel1D_set_input(const float* data_in, unsigned offset, unsigned size)
{
  for (unsigned i=0;i<size;i++)
    memory[offset + i] = data_in[i];
}

void TensorProcessorImpl::kernel1D_get_output(float* data_out, unsigned offset, unsigned size)
{
  for (unsigned i=0;i<size;i++)
    data_out[i] = memory[offset + i];
}

void TensorProcessorImpl::kernel1D_copy(float *data, unsigned steps, unsigned from, unsigned to, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++) 
    data[B.offset + to + i] = 1.0f*data[A.offset + from + i];
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
void TensorProcessorImpl::kernel2D_sub(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C) // C = A * B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] - data[B.offset + j];
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
void TensorProcessorImpl::kernel1D_pow(float *data, unsigned steps, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    data[C.offset + i] = std::pow(data[A.offset + i], data[B.offset]);
}
void TensorProcessorImpl::kernel1D_sin(float *data, unsigned steps, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = std::sin(data[A.offset + i]);
}
void TensorProcessorImpl::kernel1D_cos(float *data, unsigned steps, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = std::cos(data[A.offset + i]);
}
void TensorProcessorImpl::kernel1D_log(float *data, unsigned steps, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = log2(data[A.offset + i])*0.69314718056f; //slicer's libraries don't have log for some reason
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
void TensorProcessorImpl::kernel1D_osum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = 0;
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] += data[A.offset + j * steps + i];
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
        data[C.offset + s*A_len*B_len + i*B_len + j] = data[A.offset + s*A_len + i]*data[B.offset + s*B_len + j];
}

void TensorProcessorImpl::kernel2D_outer_p_add(float *data, unsigned steps, unsigned A_len, unsigned B_len, 
                                         Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < A_len; i++)
  {
    for (unsigned j = 0; j < B_len; j++)
    {
      data[C.offset + i*B_len + j] = 0;
      for (unsigned step = 0; step < steps; step++)
        data[C.offset + i*B_len + j] += data[A.offset + step*A_len + i]*data[B.offset + step*B_len + j];
    }
  }
}