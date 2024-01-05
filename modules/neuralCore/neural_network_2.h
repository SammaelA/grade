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
  class Layer2
  {
  public:
    std::vector<unsigned> input_shape, output_shape;
    std::vector<TensorToken> weights;
    std::vector<TensorToken> dLoss_dWeights;

    virtual void init() {};
    virtual int parameters_count() { return 0; }
    virtual TensorToken forward(const TensorToken &in) = 0;
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) = 0;
  };

  class DenseLayer2 : public Layer2
  {
  public:
    DenseLayer2(int input_size, int output_size)
    {
      input_shape.push_back(input_size);
      output_shape.push_back(output_size);
    }
    virtual void init() override;
    virtual int parameters_count() { return (input_shape[0]+1)*output_shape[0]; }
    virtual TensorToken forward(const TensorToken &in) override;
    virtual TensorToken backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput) override;
  };

  class NeuralNetwork2
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

    void add_layer(std::shared_ptr<Layer2> layer, WeightsInitializer initializer = ZERO);
    bool check_validity();
    void initialize();
    void initialize_with_weights(const float *weights);
    void initialize_from_file(std::string filename);
    void save_weights_to_file(std::string filename);
    void print_info();
    TensorProgram get_train_prog(int batch_size, Opt optimizer, Loss loss, float lr);
    void train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
               int batch_size, int iterations, Opt optimizer, Loss loss, float lr = 0.1f);
    void get_evaluate_prog();
    void evaluate(std::vector<float> &input_data, std::vector<float> &output_data);
    //float test(const TensorView &input, const TensorView &target_output, Loss loss);
    //float calculate_loss(Loss loss, const TensorView &values, const TensorView &target_values, TensorView dLoss_dValues);
    //float calculate_score(Loss loss, const TensorView &values, const TensorView &target_values);
    NeuralNetwork2(){};
    NeuralNetwork2(const NeuralNetwork2 &other) = delete;
    NeuralNetwork2 &operator=(const NeuralNetwork2 &other) = delete;

  private:
    unsigned batch_size_evaluate = 256;
    std::vector<std::shared_ptr<Layer2>> layers;
    std::vector<WeightsInitializer> initializers;
    std::vector<float> weights;
    unsigned total_params = 0;
    TensorProgram evaluate_prog;
    TensorProcessor tp;
  };
}