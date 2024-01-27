#pragma once
#include "tensor_processor.h"

extern std::vector<float> _stat_time_cmd_num;
extern std::vector<float> _stat_time_cmd_id;
extern int _stat_execution_times;

class TensorProcessorImpl
{
public:
  friend class nn::TensorProcessor;
  using Command = nn::TensorProgram::Command;
  using Variable = nn::TensorProgram::Variable;

  TensorProcessorImpl(){};
  virtual ~TensorProcessorImpl(){};

protected:
  std::vector<float> memory;

  virtual void allocate_memory(unsigned size);
  virtual void set_input(const float* in __attribute__((size("size"))), unsigned offset, unsigned size);
  virtual void get_output(float* out __attribute__((size("size"))), unsigned offset, unsigned size);
  virtual void process(const nn::TensorProgram &program);

  virtual void CommitDeviceData() {}                                                         // will be overriden in generated class
  virtual void GetExecutionTime(const char *a_funcName, float a_out[4]) { a_out[0] = 0.0f; } // will be overriden in generated class
  virtual void __attribute__((noinline)) kernel1D_fill(float *data, unsigned steps, Variable A, float val); // A = val
  virtual void __attribute__((noinline)) kernel1D_copy(float *data, unsigned steps, unsigned from, unsigned to, Variable A, Variable B); 
  virtual void __attribute__((noinline)) kernel1D_pad(float *data, unsigned steps, unsigned step_size, unsigned left_pad, unsigned right_pad, 
                                                      Variable A, Variable B); 

  virtual void __attribute__((noinline)) kernel2D_add(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A + B
  virtual void __attribute__((noinline)) kernel2D_mul(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A * B
  virtual void __attribute__((noinline)) kernel2D_sub(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A - B
  virtual void __attribute__((noinline)) kernel2D_div(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A / B
  virtual void __attribute__((noinline)) kernel2D_greater(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A / B
  virtual void __attribute__((noinline)) kernel2D_less(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A * B
  virtual void __attribute__((noinline)) kernel2D_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A - B
  virtual void __attribute__((noinline)) kernel2D_greater_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A / B
  virtual void __attribute__((noinline)) kernel2D_less_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A / B
  virtual void __attribute__((noinline)) kernel2D_not_equal(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A * B
  virtual void __attribute__((noinline)) kernel2D_logical_or(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A - B
  virtual void __attribute__((noinline)) kernel2D_logical_and(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A / B
  virtual void __attribute__((noinline)) kernel2D_where(float *data, unsigned steps, unsigned step_size, unsigned B_outer_step, 
                                                      unsigned B_inner_step, Variable A, Variable B, Variable C); // C = A / B
                                                        
  virtual void __attribute__((noinline)) kernel1D_exp(float *data, unsigned steps, Variable A, Variable B);                                 // B = exp(A)
  virtual void __attribute__((noinline)) kernel1D_pow(float *data, unsigned steps, Variable A, Variable B, Variable C);                     // C = pow(A, B)
  virtual void __attribute__((noinline)) kernel1D_sin(float *data, unsigned steps, Variable A, Variable B);                                 // B = sin(A)
  virtual void __attribute__((noinline)) kernel1D_cos(float *data, unsigned steps, Variable A, Variable B);                                 // B = cos(A)
  virtual void __attribute__((noinline)) kernel1D_log(float *data, unsigned steps, Variable A, Variable B);                                 // B = log(A)
  virtual void __attribute__((noinline)) kernel1D_not(float *data, unsigned steps, Variable A, Variable B);
   
  virtual void __attribute__((noinline)) kernel1D_sum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B);             // B = sum(A)
  virtual void __attribute__((noinline)) kernel1D_osum(float *data, unsigned steps, unsigned step_size, Variable A, Variable B);             // B = sum(A)
  virtual void __attribute__((noinline)) kernel1D_min(float *data, unsigned steps, unsigned step_size, Variable A, Variable B); 
  virtual void __attribute__((noinline)) kernel1D_max(float *data, unsigned steps, unsigned step_size, Variable A, Variable B); 

  virtual void __attribute__((noinline)) kernel2D_transpose(float *data, unsigned steps, unsigned row_len, unsigned col_len, Variable A, Variable B); // B = (A)^T
  virtual void __attribute__((noinline)) kernel2D_matmul_transposed(float *data, unsigned A_row_len, unsigned A_col_len, unsigned B_col_len, 
                                          Variable A, Variable B, Variable C);                                // C = A * (B)^T
  virtual void __attribute__((noinline)) kernel2D_outer_product(float *data, unsigned steps, unsigned A_len, unsigned B_len, 
                                         Variable A, Variable B, Variable C);
  virtual void __attribute__((noinline)) kernel1D_smax_diff(float *data, unsigned steps, unsigned step_size, Variable A, Variable B, Variable C); 
  virtual void __attribute__((noinline)) kernel3D_conv2d(float *data, int steps, int x_steps, int y_steps, int stride, int in_channels, 
                                                         int out_channels, Variable A, Variable B, Variable C); 

  virtual void __attribute__((noinline)) kernel1D_set_input(const float* data_in, unsigned offset, unsigned a_size);
  virtual void __attribute__((noinline)) kernel1D_get_output(float* data_out, unsigned offset, unsigned a_size);
};