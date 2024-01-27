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
    unsigned arg3 = program.commands[i].args[6];
    unsigned arg4 = program.commands[i].args[7];

    switch (program.commands[i].type)
    {
    case nn::TensorProgram::NOOP:
      break;
    case nn::TensorProgram::ADD:
      kernel2D_add(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::MUL:
      kernel2D_mul(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::SUB:
      kernel2D_sub(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::DIV:
      kernel2D_div(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::GREATER:
      kernel2D_greater(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::LESS:
      kernel2D_less(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::EQUAL:
      kernel2D_equal(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::GE:
      kernel2D_greater_equal(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::LE:
      kernel2D_less_equal(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::NE:
      kernel2D_not_equal(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::OR:
      kernel2D_logical_or(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::AND:
      kernel2D_logical_and(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
      break;
    case nn::TensorProgram::WHERE:
      kernel2D_where(memory.data(), arg0, arg1, arg2, arg3, A, B, C);
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
    case nn::TensorProgram::NOT:
      kernel1D_not(memory.data(), A.total_size, A, C);
      break;
    case nn::TensorProgram::SUM:
      kernel1D_sum(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::O_SUM:
      kernel1D_osum(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::MIN:
      kernel1D_min(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::MAX:
      kernel1D_max(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
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
    case nn::TensorProgram::SMAX_D:
      kernel1D_smax_diff(memory.data(), A.sizes[A.Dim-1], A.total_size/A.sizes[A.Dim-1], A, B, C);
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
    case nn::TensorProgram::PAD:
      kernel1D_pad(memory.data(), arg0, arg1, arg2, arg3, A, C);
      break;
    case nn::TensorProgram::CONV_2D:
    {
      int in_channels = B.sizes[2] > 0 ? B.sizes[2] : 1;
      int out_channels = B.sizes[3] > 0 ? B.sizes[3] : 1;
      int images = A.total_size/(A.sizes[0]*A.sizes[1]*in_channels);
      int steps = images * out_channels;
      int x_steps = C.sizes[0];
      int y_steps = C.sizes[1];
      int stride = arg0;

      kernel1D_fill(memory.data(), C.total_size, C, 0.0f);
      kernel3D_conv2d(memory.data(), steps,x_steps, y_steps, stride, in_channels, out_channels, A, B, C);
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
  {
    float tmp = val;
    data[A.offset + i] = tmp;
  }
}
void TensorProcessorImpl::kernel1D_pad(float *data, unsigned steps, unsigned step_size, unsigned left_pad, unsigned right_pad, 
                                       Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
  {
    unsigned B_step_size = step_size + left_pad + right_pad;
    for (unsigned j = 0; j < left_pad; j++)
      data[B.offset + i * B_step_size + j] = 0;
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i * B_step_size + left_pad + j] = data[A.offset + i * step_size + j];
    for (unsigned j = 0; j < right_pad; j++)
      data[B.offset + i * B_step_size + left_pad + step_size + j] = 0;
  }
}

void TensorProcessorImpl::kernel2D_add(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                       unsigned B_inner_step, Variable A, Variable B, Variable C) // C = A * B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] + data[B.offset + i * B_outer_step + B_inner_step * j];
}
void TensorProcessorImpl::kernel2D_mul(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                       unsigned B_inner_step, Variable A, Variable B, Variable C) // C = A * B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] * data[B.offset + i * B_outer_step + B_inner_step * j];
}
void TensorProcessorImpl::kernel2D_sub(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                       unsigned B_inner_step, Variable A, Variable B, Variable C) // C = A * B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] - data[B.offset + i * B_outer_step + B_inner_step * j];
}
void TensorProcessorImpl::kernel2D_div(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                       unsigned B_inner_step, Variable A, Variable B, Variable C) // C = A * B
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] / data[B.offset + i * B_outer_step + B_inner_step * j];
}

void TensorProcessorImpl::kernel2D_greater(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                           unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] > data[B.offset + i * B_outer_step + B_inner_step * j] ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_less(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                        unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] < data[B.offset + i * B_outer_step + B_inner_step * j] ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                         unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] == data[B.offset + i * B_outer_step + B_inner_step * j] ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_greater_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                                 unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] >= data[B.offset + i * B_outer_step + B_inner_step * j] ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_less_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                              unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] <= data[B.offset + i * B_outer_step + B_inner_step * j] ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_not_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                             unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[A.offset + i * step_size + j] != data[B.offset + i * B_outer_step + B_inner_step * j] ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_logical_or(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                              unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = (data[A.offset + i * step_size + j] > 0 || data[B.offset + i * B_outer_step + B_inner_step * j] > 0) ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_logical_and(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                               unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = (data[A.offset + i * step_size + j] > 0 && data[B.offset + i * B_outer_step + B_inner_step * j] > 0) ? 1.0f : 0.0f;
}
void TensorProcessorImpl::kernel2D_where(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step,
                                         unsigned B_inner_step, Variable A, Variable B, Variable C)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < step_size; j++)
      data[C.offset + i * step_size + j] = data[B.offset + i * B_outer_step + B_inner_step * j] > 0 ? data[A.offset + i * step_size + j] : 0.0f;
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
void TensorProcessorImpl::kernel1D_not(float *data, unsigned steps, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = data[A.offset + i] > 0 ? 0.0f : 1.0f;
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
void TensorProcessorImpl::kernel1D_min(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = data[A.offset + i * step_size + 0];
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] = data[A.offset + i * step_size + j] < data[B.offset + i] ? data[A.offset + i * step_size + j] : data[B.offset + i];
  }
}
void TensorProcessorImpl::kernel1D_max(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = data[A.offset + i * step_size + 0];
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] = data[A.offset + i * step_size + j] > data[B.offset + i] ? data[A.offset + i * step_size + j] : data[B.offset + i];
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

void TensorProcessorImpl::kernel1D_smax_diff(float *data, unsigned steps, unsigned step_size, 
                                             Variable _output, Variable dLoss_dOutput, Variable dLoss_dInput)
{
      // dLoss_dInput = dloss_dOutput * dOutput/dInput = (dOutput/dInput)^T * dloss_dOutput
      //(dOutput/dInput)_ij =  output_i*(1-output_i)   if i==j
      //                     =  -output_i*output_j      if i!=j
      // dOutput/dInput)^T = dOutput/dInput
      //
      //for (IndexType i = 0; i < input.total_size; i++)
      //{
      //  dLoss_dInput.get(i) = 0;
      //  for (IndexType j = 0; j < i; j++)
      //    dLoss_dInput.get(i) += -output.get(i) * output.get(j) * dLoss_dOutput.get(j);
      //  dLoss_dInput.get(i) += output.get(i) * (1 - output.get(i)) * dLoss_dOutput.get(i);
      //  for (IndexType j = i + 1; j < input.total_size; j++)
      //    dLoss_dInput.get(i) += -output.get(i) * output.get(j) * dLoss_dOutput.get(j);
      //}
  for (unsigned step = 0; step < steps; step++)
  {
    for (unsigned i = 0; i < step_size; i++)
    {
      unsigned o = _output.offset + step*step_size;
      data[dLoss_dInput.offset + step*step_size + i] = 0;
      for (unsigned j = 0; j < i; j++)
        data[dLoss_dInput.offset + step*step_size + i] += -data[o+i]*data[o+j]*data[dLoss_dOutput.offset + step*step_size + j];
      data[dLoss_dInput.offset + step*step_size + i] += data[o+i]*(1-data[o+i])*data[dLoss_dOutput.offset + step*step_size + i];
      for (unsigned j = i+1; j < step_size; j++)
        data[dLoss_dInput.offset + step*step_size + i] += -data[o+i]*data[o+j]*data[dLoss_dOutput.offset + step*step_size + j];
    }
  }
}

void TensorProcessorImpl::kernel3D_conv2d(float *data, int steps, int x_steps, int y_steps, int stride, int in_channels, 
                                          int out_channels, Variable A, Variable kernel, Variable res)
{
  //test nothing, assume that all sizes are valid and res is filled with zeros
  for (int step = 0; step < steps; step++)
  {
    for (int y = 0; y < y_steps; y++)
    {
      for (int x = 0; x < x_steps; x++)
      {
        int image_n = step / out_channels;
        int out_ch = step % out_channels;
        int kw = (int)(kernel.sizes[0])/2; //e.g. 1 for 3x3 and 2 for 5x5
        int kh = (int)(kernel.sizes[1])/2;
        int res_offset = (int)res.offset + (image_n * out_channels + out_ch) * x_steps * y_steps;

        for (int in_ch = 0; in_ch < in_channels; in_ch++)
        {
          int A_offset = (int)A.offset + (image_n * in_channels + in_ch) * (int)(A.sizes[0] * A.sizes[1]);
          int k_offset = (int)kernel.offset + (out_ch * in_channels + in_ch) * (2 * kh + 1) * (2 * kw + 1);
          for (int dy = -kh; dy <= kh; dy++)
            for (int dx = -kw; dx <= kw; dx++)
              data[res_offset + x_steps * y + x] += data[k_offset + (dy + kh) * (2 * kw + 1) + (dx + kw)] * data[A_offset + (stride * y + kh + dy) * ((int)A.sizes[0]) + (stride * x + kw + dx)];
          // printf("%u %u %u %f %f\n", res_offset + x_steps*y + x - res.offset, k_offset + (dy+kh)*(2*kw+1) + (dx+kw) - kernel.offset,
          // A_offset + (stride*y + kh + dy)*A.sizes[0] + (stride*x + kw + dx) - A.offset, data[k_offset + (dy+kh)*(2*kw+1) + (dx+kw)], data[A_offset + (stride*y + kh + dy)*A.sizes[0] + (stride*x + kw + dx)]);
        }
      }
    }
  }
}