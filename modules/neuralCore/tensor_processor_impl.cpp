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
  unsigned constexpr maxWorkGroupSizeY = (1<<16) - 1;
  unsigned constexpr maxWorkGroupSizeZ = (1<<16) - 1;
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
      kernel1D_add(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::MUL:
      kernel1D_mul(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::SUB:
      kernel1D_sub(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::DIV:
      kernel1D_div(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::GREATER:
      kernel1D_greater(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::LESS:
      kernel1D_less(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::EQUAL:
      kernel1D_equal(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::GE:
      kernel1D_greater_equal(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::LE:
      kernel1D_less_equal(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::NE:
      kernel1D_not_equal(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::OR:
      kernel1D_logical_or(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::AND:
      kernel1D_logical_and(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;
    case nn::TensorProgram::WHERE:
      kernel1D_where(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;      
    case nn::TensorProgram::MIN:
      kernel1D_min(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;     
    case nn::TensorProgram::MAX:
      kernel1D_max(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;     
    case nn::TensorProgram::POW:
      kernel1D_pow(memory.data(), (arg0*arg1*arg2+AGroupSize-1)/AGroupSize, arg0*arg1*arg2, arg1, arg2, arg3 == 0 ? ~0u : 0u, A, B, C);
      break;     
    case nn::TensorProgram::EXP:
      kernel1D_exp(memory.data(), A.total_size, A, C);
      break;
    case nn::TensorProgram::SQRT:
      kernel1D_sqrt(memory.data(), A.total_size, A, C);
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
    case nn::TensorProgram::MINIMUM:
      kernel1D_minimum(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::MAXIMUM:
      kernel1D_maximum(memory.data(), C.total_size, A.total_size / C.total_size, A, C);
      break;
    case nn::TensorProgram::MATMUL_T:
      #if DEBUG
      if (B.Dim == 2 && B.sizes[1] > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MATMUL_T workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", B.sizes[1]);
      #endif
      kernel2D_matmul_transposed(memory.data(), A.sizes[1], B.Dim == 2 ? B.sizes[1] : 1, A.sizes[0], A, B, C);
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
    {
      unsigned transp_dim = arg0;
      unsigned group_size = 1;
      for (unsigned d=0;d<transp_dim;d++)
        group_size *= A.sizes[d];
      unsigned steps = A.total_size/(A.sizes[transp_dim]*A.sizes[transp_dim+1]*group_size);
      unsigned row_len = A.sizes[transp_dim];
      unsigned col_len = A.sizes[transp_dim+1];
      #if DEBUG
      if (row_len > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: TRANSP workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", row_len);
      if (col_len > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: TRANSP workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", col_len);
      #endif
      kernel3D_transpose(memory.data(), steps, row_len, col_len, group_size, A, C);
    }
      break;
    case nn::TensorProgram::OUTER_P:
      #if DEBUG
      if (A.sizes[0] > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: OUTER_P workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]);
      if (B.sizes[0] > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: OUTER_P workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", B.sizes[0]);
      #endif
      kernel3D_outer_product(memory.data(), A.total_size/A.sizes[0], A.sizes[0], B.sizes[0], A, B, C);
      break;
    case nn::TensorProgram::SMAX_D:
      kernel1D_smax_diff(memory.data(), A.sizes[A.Dim-1], A.total_size/A.sizes[A.Dim-1], A, B, C);
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
      #if DEBUG
      if (x_steps > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: CONV_2D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", x_steps);
      if (y_steps > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: CONV_2D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", y_steps);
      #endif
      kernel3D_conv2d(memory.data(), steps, x_steps, y_steps, stride, in_channels, out_channels, A, B, C);
    }
      break;
    case nn::TensorProgram::FLIP:
    {
      unsigned flip_dim = arg0;
      unsigned group_size = 1;
      for (unsigned d=0;d<flip_dim;d++)
        group_size *= A.sizes[d];
      unsigned flip_size = A.sizes[flip_dim];
      unsigned steps = A.total_size/(flip_size*group_size);
      kernel1D_flip(memory.data(), steps, flip_size, group_size, A, C);
    }
      break;
    case nn::TensorProgram::MPOOL:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      kernel3D_max_pool(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, window_x, window_y, A, C);
    }
      break;
    case nn::TensorProgram::MPOOL_D:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      kernel3D_max_pool_diff(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, window_x, window_y, A, B, C);
    }
      break;
    case nn::TensorProgram::DILATE:
      kernel1D_fill(memory.data(), C.total_size, C, 0.0f);
      kernel1D_dilate(memory.data(), A.total_size, A.sizes[0], arg0, A.Dim>1?A.sizes[1]:1, arg1, A.Dim>2?A.sizes[2]:1, arg2, A, C);
      break;
    case nn::TensorProgram::URAND:
      kernel1D_urand(memory.data(), C.total_size, C, arg0 == 0 ? rand() : arg0);
      break;
    case nn::TensorProgram::CONV_3D:
    {
      int in_channels = B.sizes[3] > 0 ? B.sizes[3] : 1;
      int out_channels = B.sizes[4] > 0 ? B.sizes[4] : 1;
      int images = A.total_size/(A.sizes[0]*A.sizes[1]*A.sizes[2]*in_channels);
      int steps = images * out_channels;
      int x_steps = C.sizes[0];
      int y_steps = C.sizes[1];
      int z_steps = C.sizes[2];
      int stride = arg0;
      #if DEBUG
      if (x_steps > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: CONV_3D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", x_steps);
      if (y_steps > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: CONV_3D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", y_steps);
      #endif
      kernel3D_conv3d(memory.data(), steps, x_steps, y_steps, z_steps, stride, in_channels, out_channels, A, B, C);
    }
      break;
    case nn::TensorProgram::MPOOL_3D:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned window_z = arg2;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]*A.sizes[2]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      kernel3D_max_pool_3D(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, A.sizes[2]/window_z, window_x, window_y, window_z, A, C);
    }
      break;
    case nn::TensorProgram::MPOOL_3D_D:
    {
      unsigned window_x = arg0;
      unsigned window_y = arg1;
      unsigned window_z = arg2;
      unsigned steps = A.total_size/(A.sizes[0]*A.sizes[1]*A.sizes[2]);
      #if DEBUG
      if (A.sizes[0]/window_x > maxWorkGroupSizeY)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Y size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[0]/window_x);
      if (A.sizes[1]/window_y > maxWorkGroupSizeZ)
        fprintf(stderr, "TensorProgram: MPOOL_D workgroup Z size (%u) exceeds limit. Program won't execute correctly!\n", A.sizes[1]/window_y);
      #endif
      kernel3D_max_pool_3D_diff(memory.data(), steps, A.sizes[0]/window_x, A.sizes[1]/window_y, A.sizes[2]/window_z, window_x, window_y, window_z, A, B, C);
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
    /*
    printf("cmd %d %s, C = [", i, nn::TensorProgram::cmd_properties[program.commands[i].type].name.c_str());
    bool has_nan = false;
    long double min = 1e15;
    long double max = -1e15;
    long double a = 0.0;
    long double aa = 0.0;
    for (int k=0;k<C.total_size;k++)
    {
      min = std::min(min, (long double)memory[C.offset+k]);
      max = std::max(max, (long double)memory[C.offset+k]);
      a += memory[C.offset+k];
      aa += std::abs(memory[C.offset+k]);
      if (memory[C.offset+k] != memory[C.offset+k])
        has_nan = true;
    }
    printf("min max av abs_av has_nan %f %f %f %f %d]\n",(float)min, (float)max,(float)a/C.total_size,(float)aa/C.total_size,(int)has_nan);*/

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
void TensorProcessorImpl::kernel1D_flip(float *data, unsigned steps, unsigned flip_size, unsigned group_size, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
    for (unsigned j = 0; j < flip_size; j++)
      for (unsigned k = 0; k < group_size; k++)
        data[B.offset + i*flip_size*group_size + j*group_size + k] = data[A.offset + i*flip_size*group_size + (flip_size-j-1)*group_size + k];
}
void TensorProcessorImpl::kernel1D_dilate(float *data, unsigned steps, unsigned x_size, unsigned x_dilate, unsigned y_size, unsigned y_dilate,
                                          unsigned z_size, unsigned z_dilate, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
  {
    unsigned x = i % x_size;
    unsigned y = i/x_size%y_size;
    unsigned z = i/(x_size*y_size) % z_size;
    unsigned S = i/(x_size*y_size*z_size);

    unsigned new_x_size = x_size + (x_size-1)*x_dilate;
    unsigned new_y_size = y_size + (y_size-1)*y_dilate;
    unsigned new_z_size = z_size + (z_size-1)*z_dilate;
    unsigned B_off = S*new_x_size*new_y_size*new_z_size + z*(z_dilate+1)*new_x_size*new_y_size + y*(y_dilate+1)*new_x_size + x*(x_dilate+1);
    data[B.offset + B_off] = data[A.offset + i];
  }
}
void TensorProcessorImpl::kernel1D_urand(float *data, unsigned steps, Variable A, unsigned seed)
{
  for (unsigned i = 0; i < steps; i++)
  {
    unsigned n = seed + i;
    n = (n << 13) ^ n; 
    data[A.offset + i] = float((n * (n*n*15731+789221) + 1376312589) & 0x3fffffff) / float(0x3fffffff);
  }
}

void TensorProcessorImpl::kernel1D_add(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul, 
                                       Variable A, Variable B, Variable C) // C = A * B
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] + data[B.offset + Bi];
    }
  }
}
void TensorProcessorImpl::kernel1D_mul(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C) // C = A * B
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] * data[B.offset + Bi];
    }
  }
}
void TensorProcessorImpl::kernel1D_sub(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul, 
                                       Variable A, Variable B, Variable C) // C = A * B
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] - data[B.offset + Bi];
    }
  }
}
void TensorProcessorImpl::kernel1D_div(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul, 
                                       Variable A, Variable B, Variable C) // C = A * B
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] / data[B.offset + Bi];
    }
  }
}

void TensorProcessorImpl::kernel1D_greater(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                           Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] > data[B.offset + Bi] ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_less(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                        Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] < data[B.offset + Bi] ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_equal(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                         Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai]  == data[B.offset + Bi] ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_greater_equal(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                                 Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] >= data[B.offset + Bi] ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_less_equal(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                              Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] <= data[B.offset + Bi] ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_not_equal(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                             Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] != data[B.offset + Bi] ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_logical_or(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                              Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] > 0 || data[B.offset + Bi] > 0 ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_logical_and(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                               Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[A.offset + Ai] > 0 && data[B.offset + Bi] > 0 ? 1.0f : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_where(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                         Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = data[B.offset + Bi] > 0 ? data[A.offset + Ai] : 0.0f;
    }
  }
}
void TensorProcessorImpl::kernel1D_min(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = std::min(data[A.offset + Ai], data[B.offset + Bi]);
    }
  }
}
void TensorProcessorImpl::kernel1D_max(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = std::max(data[A.offset + Ai], data[B.offset + Bi]);
    }
  }
}
void TensorProcessorImpl::kernel1D_pow(float *data, unsigned steps, unsigned total_size, unsigned step_size, unsigned group_size, unsigned Ai_mul,
                                       Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned s1 = 0; s1 < steps; s1++)
  {
    for (unsigned s2 = 0; s2 < std::min(AGroupSize, total_size-s1*AGroupSize); s2++)
    {
    unsigned s = s1*AGroupSize + s2;
    unsigned Ai = Ai_mul & s;
    unsigned Bi = s / group_size % step_size;
    unsigned Ci = s;
    data[C.offset + Ci] = std::pow(data[A.offset + Ai], data[B.offset + Bi]);
    }
  }
}

void TensorProcessorImpl::kernel1D_exp(float *data, unsigned steps, Variable A, Variable B) // B = exp(A)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = std::exp(data[A.offset + i]);
}
void TensorProcessorImpl::kernel1D_sqrt(float *data, unsigned steps, Variable A, Variable B)
{
  for (unsigned i = 0; i < steps; i++)
    data[B.offset + i] = std::sqrt(data[A.offset + i]);
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
  #pragma omp parallel for
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = 0;
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] += data[A.offset + j * steps + i];
  }
}
void TensorProcessorImpl::kernel1D_minimum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = data[A.offset + i * step_size + 0];
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] = data[A.offset + i * step_size + j] < data[B.offset + i] ? data[A.offset + i * step_size + j] : data[B.offset + i];
  }
}
void TensorProcessorImpl::kernel1D_maximum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B) // B = sum(A)
{
  for (unsigned i = 0; i < steps; i++)
  {
    data[B.offset + i] = data[A.offset + i * step_size + 0];
    for (unsigned j = 0; j < step_size; j++)
      data[B.offset + i] = data[A.offset + i * step_size + j] > data[B.offset + i] ? data[A.offset + i * step_size + j] : data[B.offset + i];
  }
}

void TensorProcessorImpl::kernel3D_transpose(float *data, unsigned steps, unsigned row_len, unsigned col_len, unsigned group_size, Variable A, Variable B)
{
  for (unsigned s = 0; s < steps; s++)
    for (unsigned i = 0; i < row_len; i++)
      for (unsigned j = 0; j < col_len; j++)
        for (unsigned k = 0; k < group_size; k++)
          data[B.offset + (s*col_len*row_len + i*col_len + j)*group_size + k] = data[A.offset + (s*col_len*row_len + j*row_len + i)*group_size + k];
}

void TensorProcessorImpl::kernel2D_matmul_transposed(float *data, unsigned A_col_len, unsigned B_col_len, 
                                                     unsigned A_row_len, Variable A, Variable B, Variable C)
{
  #pragma omp parallel for
  for (unsigned i = 0; i < A_col_len; i++)
  {
    for (unsigned j = 0; j < B_col_len; j++)
    {
      data[C.offset + i*B_col_len + j] = 0;
      for (unsigned k = 0; k < A_row_len; k++)
        data[C.offset + i*B_col_len + j] += data[A.offset + i*A_row_len + k]*data[B.offset + j*A_row_len + k];
    }
  }
}

void TensorProcessorImpl::kernel3D_outer_product(float *data, unsigned steps, unsigned A_len, unsigned B_len, 
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
        int kw2 = (int)(kernel.sizes[0]);
        int kh2 = (int)(kernel.sizes[1]);
        int res_offset = (int)res.offset + (image_n * out_channels + out_ch) * x_steps * y_steps;

        data[res_offset + x_steps * y + x] = 0;
        for (int in_ch = 0; in_ch < in_channels; in_ch++)
        {
          int A_offset = (int)A.offset + (image_n * in_channels + in_ch) * (int)(A.sizes[0] * A.sizes[1]);
          int k_offset = (int)kernel.offset + (out_ch * in_channels + in_ch) * kh2 * kw2;
          for (int dy = -kh; dy < kh2 - kh; dy++)
            for (int dx = -kw; dx < kw2 - kw; dx++)
            {
              float t = data[k_offset + (dy + kh) * kw2 + (dx + kw)] * data[A_offset + (stride * y + kh + dy) * ((int)A.sizes[0]) + (stride * x + kw + dx)];
            data[res_offset + x_steps * y + x] += t;
            }
        }
      }
    }
  }
}

void TensorProcessorImpl::kernel3D_conv3d(float *data, int steps, int x_steps, int y_steps, int z_steps, int stride, 
                                          int in_channels, int out_channels, Variable A, Variable kernel, Variable res)
{
  //test nothing, assume that all sizes are valid and res is filled with zeros
  for (int step = 0; step < steps; step++)
  {
    for (int z = 0; z < z_steps; z++)
    {
      for (int y = 0; y < y_steps; y++)
      {
        for (int x = 0; x < x_steps; x++)
        {
          int image_n = step / out_channels;
          int out_ch = step % out_channels;
          int kw = (int)(kernel.sizes[0])/2; //e.g. 1 for 3x3 and 2 for 5x5
          int kh = (int)(kernel.sizes[1])/2;
          int kd = (int)(kernel.sizes[2])/2;
          int kw2 = (int)(kernel.sizes[0]);
          int kh2 = (int)(kernel.sizes[1]);
          int kd2 = (int)(kernel.sizes[2]);
          int res_offset = (int)res.offset + (image_n * out_channels + out_ch) * x_steps * y_steps * z_steps;

          data[res_offset + x_steps*y_steps*z + x_steps*y + x] = 0;
          for (int in_ch = 0; in_ch < in_channels; in_ch++)
          {
            int A_offset = (int)A.offset + (image_n * in_channels + in_ch) * (int)(A.sizes[0] * A.sizes[1] * A.sizes[2]);
            int k_offset = (int)kernel.offset + (out_ch * in_channels + in_ch) * kh2 * kw2 * kd2;
            for (int dz = -kd; dz < kd2 - kd; dz++)
              for (int dy = -kh; dy < kh2 - kh; dy++)
                for (int dx = -kw; dx < kw2 - kw; dx++)
                {
                  float t = data[k_offset + (dz+kd)*kh2*kw2 + (dy+kh)*kw2 + (dx+kw)] * 
                            data[A_offset + (stride*z+kd+dz)*(int)(A.sizes[0]*A.sizes[1]) + (stride*y+kh+dy)*((int)A.sizes[0]) + (stride*x+kw+dx)];
                  data[res_offset + x_steps*y_steps*z + x_steps*y + x] += t;
                }
          }
        }
      }
    }
  }
}

void TensorProcessorImpl::kernel3D_max_pool(float *data, int steps, int x_steps, int y_steps, int window_x, int window_y, 
                                            Variable A, Variable res)
{
  for (int step = 0; step < steps; step++)
  {
    for (int y = 0; y < y_steps; y++)
    {
      for (int x = 0; x < x_steps; x++)
      {
        float max_val = data[A.offset + step*x_steps*window_x*y_steps*window_y + (y*window_y+0)*x_steps + x*window_x + 0];
        for (int wy=0;wy<window_y;wy++)
        {
          for (int wx=0;wx<window_x;wx++)
          {
            float val = data[A.offset + step*x_steps*window_x*y_steps*window_y + (y*window_y+wy)*x_steps*window_x + x*window_x + wx];
            if (val > max_val)
              max_val = val;
          }
        }
        data[res.offset + step*x_steps*y_steps + y*x_steps + x] = max_val;
      }
    }
  }
}

void TensorProcessorImpl::kernel3D_max_pool_diff(float *data, int steps, int x_steps, int y_steps, int window_x, int window_y,
                                                 Variable A, Variable dLoss_dOutput, Variable dLoss_dInput)
{
  for (int step = 0; step < steps; step++)
  {
    for (int y = 0; y < y_steps; y++)
    {
      for (int x = 0; x < x_steps; x++)
      {
        float max_val = data[A.offset + step*x_steps*window_x*y_steps*window_y + (y*window_y+0)*x_steps*window_x + x*window_x + 0];
        int max_wy = 0;
        int max_wx = 0;
        for (int wy=0;wy<window_y;wy++)
        {
          for (int wx=0;wx<window_x;wx++)
          {
            data[dLoss_dInput.offset + step*x_steps*window_x*y_steps*window_y + (y*window_y+wy)*x_steps*window_x + x*window_x + wx] = 0;
            float val = data[A.offset + step*x_steps*window_x*y_steps*window_y + (y*window_y+wy)*x_steps*window_x + x*window_x + wx];
            if (val > max_val)
            {
              max_val = val;
              max_wx = wx;
              max_wy = wy;
            }
          }
        }
        data[dLoss_dInput.offset + step*x_steps*window_x*y_steps*window_y + (y*window_y+max_wy)*x_steps*window_x + x*window_x + max_wx] = 
            data[dLoss_dOutput.offset + step*x_steps*y_steps + y*x_steps + x];
      }
    }
  }
}

void TensorProcessorImpl::kernel3D_max_pool_3D(float *data, int steps, int x_steps, int y_steps, int z_steps, 
                                               int window_x, int window_y, int window_z, Variable A, Variable res)
{
  for (int step = 0; step < steps; step++)
  {
    for (int z = 0; z < z_steps; z++)
    {
      for (int y = 0; y < y_steps; y++)
      {
        for (int x = 0; x < x_steps; x++)
        {
          unsigned offset = A.offset + step*x_steps*window_x*y_steps*window_y*z_steps*window_z;
          float max_val = data[offset + (z*window_z+0)*y_steps*x_steps*window_x*window_y + (y*window_y+0)*x_steps*window_x + x*window_x + 0];
          for (int wz=0;wz<window_z;wz++)
          {
            for (int wy=0;wy<window_y;wy++)
            {
              for (int wx=0;wx<window_x;wx++)
              {
                float val = data[offset + (z*window_z+wz)*y_steps*x_steps*window_x*window_y + (y*window_y+wy)*x_steps*window_x + x*window_x + wx];
                if (val > max_val)
                  max_val = val;
              }
            }
          }
          data[res.offset + step*x_steps*y_steps*z_steps + z*y_steps*x_steps + y*x_steps + x] = max_val;
        }
      }
    }
  }
}

void TensorProcessorImpl::kernel3D_max_pool_3D_diff(float *data, int steps, int x_steps, int y_steps, int z_steps, 
                                                    int window_x, int window_y, int window_z, 
                                                    Variable A, Variable dLoss_dOutput, Variable dLoss_dInput)
{
  for (int step = 0; step < steps; step++)
  {
    for (int z = 0; z < z_steps; z++)
    {
      for (int y = 0; y < y_steps; y++)
      {
        for (int x = 0; x < x_steps; x++)
        {
          unsigned offset = A.offset + step*x_steps*window_x*y_steps*window_y*z_steps*window_z;
          float max_val = data[offset + (z*window_z+0)*y_steps*x_steps*window_x*window_y + (y*window_y+0)*x_steps*window_x + x*window_x + 0];
          int max_wz = 0;
          int max_wy = 0;
          int max_wx = 0;
          for (int wz=0;wz<window_z;wz++)
          {
            for (int wy=0;wy<window_y;wy++)
            {
              for (int wx=0;wx<window_x;wx++)
              {
                data[dLoss_dInput.offset + step*x_steps*window_x*y_steps*window_y*z_steps*window_z + 
                    (z*window_z+wz)*y_steps*window_y*x_steps*window_x +
                    (y*window_y+wy)*x_steps*window_x + x*window_x + wx] = 0;
                float val = data[offset + (z*window_z+wz)*y_steps*x_steps*window_x*window_y + (y*window_y+wy)*x_steps*window_x + x*window_x + wx];
                if (val > max_val)
                {
                  max_val = val;
                  max_wx = wx;
                  max_wy = wy;
                  max_wz = wz;
                }
              }
            }
          }
          data[dLoss_dInput.offset + step*x_steps*window_x*y_steps*window_y*z_steps*window_z + 
               (z*window_z+max_wz)*y_steps*window_y*x_steps*window_x +
               (y*window_y+max_wy)*x_steps*window_x + x*window_x + max_wx] = 
              data[dLoss_dOutput.offset + step*x_steps*y_steps*z_steps + z*y_steps*x_steps + y*x_steps + x];
        }
      }
    }
  }
}
