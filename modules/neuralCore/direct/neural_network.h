#pragma once
#include "tensors.h"
#include <vector>
#include <cstdint>
#include <array>
#include <cmath>
#include <memory>
#include <functional>

namespace nnd
{
  class Layer
  {
  public:
    Shape input_shape, output_shape;

    virtual ~Layer() {};
    virtual void init(float *param_mem, float *gradient_mem, float *tmp_mem, bool initialize_random_weights) = 0;
    virtual void forward(const TensorView &input, TensorView &output) = 0;
    virtual void backward(const TensorView &input, const TensorView &output,
                          const TensorView &dLoss_dOutput, TensorView dLoss_dInput, bool first_layer) = 0;
    virtual int parameters_memory_size() = 0;
    virtual int tmp_memory_size() = 0;
  };

  class DenseLayer : public Layer
  {
  public:
    DenseLayer() = default;
    DenseLayer(int input_size, int output_size);
    virtual void init(float *param_mem, float *gradient_mem, float *tmp_mem, bool initialize_random_weights) override;
    virtual void forward(const TensorView &input, TensorView &output) override;
    virtual void backward(const TensorView &input, const TensorView &output,
                          const TensorView &dLoss_dOutput, TensorView dLoss_dInput, bool first_layer) override;
    virtual int parameters_memory_size() override;
    virtual int tmp_memory_size() override;

  private:
    TensorView A, b;
    TensorView dLoss_dA, dLoss_db;
    TensorView At, op;
  };
  class ReLULayer : public Layer
  {
  public:
    ReLULayer() = default;
    virtual void init(float *param_mem, float *gradient_mem, float *tmp_mem, bool initialize_random_weights) override
    {
    }

    virtual void forward(const TensorView &input, TensorView &output) override
    {
      for (IndexType i = 0; i < input.total_size; i++)
        output.get(i) = (input.get(i) >= 0) ? input.get(i) : 0.0f;
    }
    virtual void backward(const TensorView &input, const TensorView &output,
                          const TensorView &dLoss_dOutput, TensorView dLoss_dInput, bool first_layer) override
    {
      for (IndexType i = 0; i < input.total_size; i++)
        dLoss_dInput.get(i) = (input.get(i) >= 0) ? dLoss_dOutput.get(i) : 0.0f;
    }
    virtual int parameters_memory_size() override
    {
      return 0;
    }
    virtual int tmp_memory_size() override
    {
      return 0;
    }
  };

  class SoftMaxLayer : public Layer
  {
  public:
    SoftMaxLayer() = default;
    virtual void init(float *param_mem, float *gradient_mem, float *tmp_mem, bool initialize_random_weights) override
    {
    }

    virtual void forward(const TensorView &input, TensorView &output) override
    {
      float mx = input.get(0);
      for (IndexType i = 0; i < input.total_size; i++)
        mx = std::max(mx, input.get(i));
      float sum = 0;
      for (IndexType i = 0; i < input.total_size; i++)
      {
        float ex = exp(input.get(i) - mx + 1e-15);
        output.get(i) = ex;
        sum += ex;
      }
      float mult = 1 / sum;
      for (IndexType i = 0; i < input.total_size; i++)
        output.get(i) *= mult;
    }
    virtual void backward(const TensorView &input, const TensorView &output,
                          const TensorView &dLoss_dOutput, TensorView dLoss_dInput, bool first_layer) override
    {
      // dLoss_dInput = dloss_dOutput * dOutput/dInput = (dOutput/dInput)^T * dloss_dOutput
      //(dOutput/dInput)_ij =  output_i*(1-output_i)   if i==j
      //                     =  -output_i*output_j      if i!=j
      // dOutput/dInput)^T = dOutput/dInput
      //
      for (IndexType i = 0; i < input.total_size; i++)
      {
        dLoss_dInput.get(i) = 0;
        for (IndexType j = 0; j < i; j++)
          dLoss_dInput.get(i) += -output.get(i) * output.get(j) * dLoss_dOutput.get(j);
        dLoss_dInput.get(i) += output.get(i) * (1 - output.get(i)) * dLoss_dOutput.get(i);
        for (IndexType j = i + 1; j < input.total_size; j++)
          dLoss_dInput.get(i) += -output.get(i) * output.get(j) * dLoss_dOutput.get(j);
      }
    }
    virtual int parameters_memory_size() override
    {
      return 0;
    }
    virtual int tmp_memory_size() override
    {
      return 0;
    }
  };

  class SinLayer : public Layer
  {
  public:
    static constexpr float omega_0 = 30;
    SinLayer() = default;
    virtual void init(float *param_mem, float *gradient_mem, float *tmp_mem, bool initialize_random_weights) override
    {
    }

    virtual void forward(const TensorView &input, TensorView &output) override
    {
      for (IndexType i = 0; i < input.total_size; i++)
        output.get(i) = sinf(omega_0 * input.get(i));
    }
    virtual void backward(const TensorView &input, const TensorView &output,
                          const TensorView &dLoss_dOutput, TensorView dLoss_dInput, bool first_layer) override
    {
      for (IndexType i = 0; i < input.total_size; i++)
        dLoss_dInput.get(i) = dLoss_dOutput.get(i) * omega_0 * cosf(omega_0 * input.get(i));
    }
    virtual int parameters_memory_size() override
    {
      return 0;
    }
    virtual int tmp_memory_size() override
    {
      return 0;
    }
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

    using LossFunction = std::function<float(const TensorView &, const TensorView &, TensorView)>;

    void add_layer(std::shared_ptr<Layer> layer);
    bool check_validity();
    void initialize();
    void initialize_with_weights(const float *weights);
    void initialize_from_file(std::string filename);
    void save_weights_to_file(std::string filename);
    void print_info();
    void train(const TensorView &inputs /*[input_size, count]*/, const TensorView &outputs /*[output_size, count]*/,
               const TensorView &inputs_val, const TensorView &outputs_val,
               int batch_size, int iterations, Opt optimizer, Loss loss, float lr = 0.1f);
    void evaluate(const TensorView &input, TensorView output);
    float test(const TensorView &input, const TensorView &target_output, Loss loss);
    float calculate_loss(Loss loss, const TensorView &values, const TensorView &target_values, TensorView dLoss_dValues);
    float calculate_score(Loss loss, const TensorView &values, const TensorView &target_values);
    NeuralNetwork(){};
    NeuralNetwork(const NeuralNetwork &other) = delete;
    NeuralNetwork &operator=(const NeuralNetwork &other) = delete;

  private:
    std::vector<std::shared_ptr<Layer>> layers;

    std::vector<float> weights;
    std::vector<float> tmp_mem;

    std::vector<TensorView> layer_outputs;
  };
}