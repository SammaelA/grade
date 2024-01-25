#pragma once
#include "tensor_compiler.h"
#include <vector>
#include <cstdint>
#include <array>
#include <cmath>
#include <memory>
#include <functional>

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
      TensorToken max_val = input.max(input.Dim-1) + 1e-15f;
      TensorToken output = TensorToken::g_2op(TensorProgram::SUB, input, max_val, 
                                              max_val.total_size(), input.total_size()/max_val.total_size(), 1, 0);
      output = TensorToken::exp(output);
      TensorToken sum = output.sum(input.Dim-1);
      TensorToken res = TensorToken::g_2op(TensorProgram::DIV, output, sum, 
                                           sum.total_size(), output.total_size()/sum.total_size(), 1, 0);
      return res;
    }
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override
    {
      //return dLoss_dOutput * TensorToken::cos(input*mult) * mult;
      TensorToken dLoss_dInput(input.sizes);
      TensorToken::issue_command(TensorProgram::SMAX_D, output, dLoss_dOutput, dLoss_dInput);
      return dLoss_dInput;
    }
    virtual std::string get_name() override { return "SoftMax"; }
  };

  class NeuralNetwork
  {
  public:
    enum Opt
    {
      GD,
      Adam
    };
    enum Loss
    {
      MSE,
      CrossEntropy
    };
    enum WeightsInitializer
    {
      ZERO,
      HE,
      SIREN
    };

    void add_layer(std::shared_ptr<Layer> layer, WeightsInitializer initializer = ZERO);
    void set_batch_size_for_evaluate(int size);
    bool check_validity();
    void initialize();
    void initialize_with_weights(const float *weights);
    void initialize_from_file(std::string filename);
    void save_weights_to_file(std::string filename);
    void set_arch_to_file(std::string filename);
    void print_info();
    TensorProgram get_train_prog(int batch_size, Opt optimizer, Loss loss, float lr);
    void train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
               int batch_size, int iterations, Opt optimizer, Loss loss, float lr = 0.1f, bool verbose = false);
    void get_evaluate_prog();
    void evaluate(std::vector<float> &input_data, std::vector<float> &output_data, int samples = -1);
    //float test(const TensorView &input, const TensorView &target_output, Loss loss);
    //float calculate_loss(Loss loss, const TensorView &values, const TensorView &target_values, TensorView dLoss_dValues);
    //float calculate_score(Loss loss, const TensorView &values, const TensorView &target_values);
    NeuralNetwork(){};
    NeuralNetwork(const NeuralNetwork &other) = delete;
    NeuralNetwork &operator=(const NeuralNetwork &other) = delete;

  private:
    bool initialized = false;
    unsigned batch_size_evaluate = 256;
    std::vector<std::shared_ptr<Layer>> layers;
    std::vector<WeightsInitializer> initializers;
    std::vector<float> weights;
    unsigned total_params = 0;
    TensorProgram evaluate_prog;
  };
}