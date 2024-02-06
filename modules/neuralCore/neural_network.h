#pragma once
#include "tensor_compiler.h"
#include <vector>
#include <cstdint>
#include <array>
#include <cmath>
#include <memory>
#include <functional>
#include <tuple>
#include <variant>

namespace nn
{
  class Layer
  {
  public:
    std::vector<unsigned> input_shape, output_shape;
    std::vector<TensorToken> weights;
    std::vector<TensorToken> dLoss_dWeights;

    virtual ~Layer() {};
    virtual void init() {};
    virtual int parameters_count() { return 0; }
    virtual TensorToken forward(const TensorToken &in) = 0;
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) = 0;
    virtual std::string get_name() = 0;
  };

  class DenseLayer : public Layer
  {
  public:
    DenseLayer(int input_size, int output_size)
    {
      input_shape.push_back(input_size);
      output_shape.push_back(output_size);
    }
    virtual void init() override;
    virtual int parameters_count() override { return (input_shape[0]+1)*output_shape[0]; };
    virtual TensorToken forward(const TensorToken &in) override;
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override;
    virtual std::string get_name() override { return "Dense"; }
  };

  class SinLayer : public Layer
  {
    float mult = 30.0f;
  public:
    SinLayer(float _mult = 30.0f){ mult = _mult; }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      return TensorToken::sin(input*mult);
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      return dLoss_dOutput * TensorToken::cos(input*mult) * mult;
    }
    virtual std::string get_name() override { return "Sin"; }
  };

  class SoftMaxLayer : public Layer
  {
  public:
    SoftMaxLayer(){ }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      TensorToken max_val = input.maximum(input.Dim-1) + 1e-15f;
      TensorToken output = TensorToken::g_2op(TensorProgram::SUB, input, max_val, 1);
      output = TensorToken::exp(output);
      TensorToken sum = output.sum(input.Dim-1) + 1e-15f;
      TensorToken res = TensorToken::g_2op(TensorProgram::DIV, output, sum, 1);
      return res;
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      TensorToken dLoss_dInput(input.sizes);
      TensorToken::issue_command(TensorProgram::SMAX_D, output, dLoss_dOutput, dLoss_dInput);
      return dLoss_dInput;
    }
    virtual std::string get_name() override { return "SoftMax"; }
  };

  class ReLULayer : public Layer
  {
  public:
    ReLULayer(){ }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      return TensorToken::g_2op(TensorProgram::WHERE, input, input); //input[i] > 0 ? input[i] : 0
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      return TensorToken::g_2op(TensorProgram::WHERE, dLoss_dOutput, input);//input[i] > 0 ? dLoss_dOutput[i] : 0
    }
    virtual std::string get_name() override { return "ReLU"; }
  };

  class TanhLayer : public Layer
  {
  public:
    TanhLayer(){ }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      TensorToken ex = TensorToken::exp(input);
      TensorToken ex2 = ex*ex;
      return (ex2 - 1.0f)/(ex2 + 1.0f); //tanh(x)
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      //output = tanh(x)
      //dtanh(x) = 1 - tanh(x)^2 
      return dLoss_dOutput*(1.0f - output*output);
    }
    virtual std::string get_name() override { return "Tanh"; }
  };

  class Conv2DLayer : public Layer
  {
  public:
    enum Padding
    {
      NO_PAD,
      SAME,
    };
    unsigned kernel_size = 3;
    unsigned stride = 1;
    Padding padding = NO_PAD;
    bool use_bias = false;
  public:
    Conv2DLayer(unsigned input_x,  unsigned input_y,  unsigned input_ch, 
                unsigned output_channels, unsigned _kernel_size = 3, unsigned _stride = 1, Padding _padding = NO_PAD, bool _use_bias = false)
    {
      assert(_kernel_size%2);
      assert(_stride == 1 || _padding == NO_PAD);
      kernel_size = _kernel_size;
      stride = _stride;
      padding = _padding;
      use_bias = _use_bias;
      input_shape = {input_x, input_y, input_ch};
      if (padding == NO_PAD)
        output_shape = {(input_x - kernel_size)/stride + 1, (input_y - kernel_size)/stride + 1, output_channels};
      else
        output_shape = {input_x, input_y, output_channels};
    }
    virtual void init() override 
    {
      weights.clear();

      weights.push_back(TensorToken(kernel_size, kernel_size, input_shape[2], output_shape[2])); // kernel
      if (use_bias)
        weights.push_back(TensorToken(output_shape[2]));                 // bias

      dLoss_dWeights.resize(weights.size());
    };
    virtual int parameters_count() override { return input_shape[2]*output_shape[2]*kernel_size*kernel_size + use_bias*output_shape[2]; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      unsigned pad = padding == NO_PAD ? 0 : (kernel_size-1)/2;
      TensorToken pad_input = (pad > 0) ? input.add_padding(pad, pad, 0).add_padding(pad, pad, 1) : input;
      TensorToken conv = TensorToken::conv2D(pad_input, weights[0], stride);
      return use_bias ? TensorToken::g_2op(TensorProgram::ADD, conv, weights[1], 2) : conv;
      //return TensorToken::g_2op(TensorProgram::ADD, TensorToken::conv2D(input, weights[0]), weights[1], input.sizes[3], input.total_size(), 0, 1);
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      unsigned OD_sizes[TensorCompiler::MAX_DIM];
      for (int i = 0; i < TensorCompiler::MAX_DIM; i++)
        OD_sizes[i] = dLoss_dOutput.sizes[i];
      OD_sizes[0] = (dLoss_dOutput.sizes[0]-1)*stride + 1;
      OD_sizes[1] = (dLoss_dOutput.sizes[1]-1)*stride + 1;
      TensorToken dLoss_dOutput_dilated = (stride == 1) ? dLoss_dOutput : TensorToken(OD_sizes);
      if (stride > 1)
        TensorToken::issue_command(TensorProgram::DILATE, dLoss_dOutput, dLoss_dOutput, dLoss_dOutput_dilated, stride-1, stride-1);
      unsigned pad = padding == NO_PAD ? 0 : (kernel_size-1)/2;
      float batch_size = (float)(input.sizes[input.Dim-1]);
      TensorToken X = input.flip(0).flip(1);
      TensorToken pad_X = (pad > 0) ? X.add_padding(pad, pad, 0).add_padding(pad, pad, 1) : X;
      dLoss_dWeights[0] = TensorToken::conv2D(pad_X.transpose(2), dLoss_dOutput_dilated.transpose(2)).transpose(2) / batch_size;
      if (use_bias)
        dLoss_dWeights[1] = dLoss_dOutput.sum(2).outer_sum() / batch_size;//dilation won't change result here

      unsigned i_pad = kernel_size - pad - 1;
      TensorToken pad_dLoss_dOutput = (i_pad > 0) ? dLoss_dOutput_dilated.add_padding(i_pad, i_pad, 0).add_padding(i_pad, i_pad, 1) : dLoss_dOutput_dilated;

      return TensorToken::conv2D(pad_dLoss_dOutput, weights[0].flip(0).flip(1).transpose(2));
    }
    virtual std::string get_name() override { return "Conv2DLayer"; }
  };

  class FlattenLayer : public Layer
  {
  public:
    FlattenLayer(unsigned input_x,  unsigned input_y,  unsigned input_ch)
    {
      input_shape = {input_x, input_y, input_ch};
      output_shape = {input_x*input_y*input_ch};
    }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      return input.reshape({(unsigned)output_shape[0], input.sizes[input.Dim-1]});
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      return dLoss_dOutput.reshape({input.sizes[0], input.sizes[1], input.sizes[2], input.sizes[3]});
    }
    virtual std::string get_name() override { return "Flatten"; }
  };

  class SigmoidLayer : public Layer
  {
  public:
    SigmoidLayer(){ }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      TensorToken one(input.sizes);
      one.fill(1.0f);
      return one/(1.0f + TensorToken::exp(-1.0f*input)); //sigmoid(x)
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      //output = sigmoid(x)
      //dsigmoid(x) = sigmoid(x)*(1-sigmoid(x))
      return dLoss_dOutput*(output*(1.0f - output));
    }
    virtual std::string get_name() override { return "Sigmoid"; }
  };

  class MaxPoolingLayer : public Layer
  {
    unsigned window_size = 2;
  public:
    MaxPoolingLayer(unsigned input_x, unsigned input_y, unsigned input_ch, unsigned _window_size = 2)
    {
      window_size = _window_size;
      input_shape = {input_x, input_y, input_ch};
      output_shape = {input_x/window_size, input_y/window_size, input_ch};
    }
    virtual void init() override { };
    virtual int parameters_count() override { return 0; };
    virtual TensorToken forward(const TensorToken &input) override
    {
      unsigned output_sizes[TensorCompiler::MAX_DIM];
      for (int i = 0; i < TensorCompiler::MAX_DIM; i++)
        output_sizes[i] = input.sizes[i];
      output_sizes[0] = output_shape[0];
      output_sizes[1] = output_shape[1];
      TensorToken output(output_sizes);
      TensorToken::issue_command(TensorProgram::MPOOL, input, input, output, window_size, window_size);
      return output;
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      TensorToken dLoss_dInput(input.sizes);
      TensorToken::issue_command(TensorProgram::MPOOL_D, input, dLoss_dOutput, dLoss_dInput, window_size, window_size);
      return dLoss_dInput;
    }
    virtual std::string get_name() override { return "MaxPooling"; }
  };

  struct OptimizerGD
  {
    explicit OptimizerGD(float _lr = 0.01) { learning_rate = _lr; }
    float learning_rate;
  };

  struct OptimizerAdam
  {
    explicit OptimizerAdam(float _lr = 0.01, float _beta_1 = 0.9, float _beta_2 = 0.999, float _eps = 1e-8) 
    { 
      learning_rate = _lr; 
      beta_1 = _beta_1;
      beta_2 = _beta_2;
      eps = _eps;
    }
    float learning_rate;
    float beta_1;
    float beta_2;
    float eps;
  };

  using Optimizer = std::variant<OptimizerGD, OptimizerAdam>;

  enum class Loss
  {
    MSE,
    CrossEntropy
  };
  enum class Initializer
  {
    Zero,
    He,
    Siren,
    GlorotNormal
  };
  enum class Metric
  {
    MSE,
    MAE,
    Accuracy,
    Precision,
    Recall,
    AUC_ROC,
    AUC_PR
  };

  class NeuralNetwork
  {
  public:
    void add_layer(std::shared_ptr<Layer> layer, Initializer initializer = Initializer::Zero);
    void set_batch_size_for_evaluate(int size);
    bool check_validity();
    void initialize();
    void initialize_with_weights(const float *weights);
    void initialize_from_file(std::string filename);
    void save_weights_to_file(std::string filename);
    void set_arch_to_file(std::string filename);
    void print_info();
    TensorProgram get_train_prog(int batch_size, Optimizer optimizer, Loss loss);
    void train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
               int batch_size, int iterations, Optimizer optimizer, Loss loss, bool verbose = false);
    void train(const float *data, const float *labels, int samples, int batch_size, int epochs, bool use_validation = false, Optimizer optimizer = OptimizerAdam(0.01f), 
               Loss loss = Loss::CrossEntropy, Metric metric = Metric::Accuracy, bool verbose = false);
    void get_evaluate_prog();
    void evaluate(std::vector<float> &input_data, std::vector<float> &output_data, int samples = -1);
    void evaluate(const float *input_data, float *output_labels, int samples);
    float calculate_metric(const float *output, const float *output_ref, int samples, Metric metric);
    NeuralNetwork(){};
    NeuralNetwork(const NeuralNetwork &other) = delete;
    NeuralNetwork &operator=(const NeuralNetwork &other) = delete;

  private:
    bool initialized = false;
    unsigned batch_size_evaluate = 256;
    std::vector<std::shared_ptr<Layer>> layers;
    std::vector<Initializer> initializers;
    std::vector<float> weights;
    unsigned total_params = 0;
    TensorProgram evaluate_prog;
  };
}