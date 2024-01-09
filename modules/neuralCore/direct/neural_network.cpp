#include "tensors.h"
#include "neural_network.h"
#include <cassert>
#include <cstdio>
#include <string>
#include <fstream>
#include <chrono>

namespace nnd
{
  void zero_initialization(TensorView t)
  {
    fill(t, 0);
  }

  void he_initialization(TensorView t, int fan_in, int fan_out)
  {
    float mn = -sqrt(6 / fan_in);
    float mx = sqrt(6 / fan_out);
    float d = mx - mn;
    for (IndexType i = 0; i < t.total_size; i++)
      t.get(i) = d * (((double)rand()) / RAND_MAX) + mn;
  }
  void SIREN_initialization(TensorView t, int fan_in, int fan_out)
  {
    float mx = (fan_in == 2) ? 0.5 : (sqrt(6.0 / fan_in) / SinLayer::omega_0);
    float mn = -mx;
    float d = mx - mn;
    for (IndexType i = 0; i < t.total_size; i++)
      t.get(i) = d * (((double)rand()) / RAND_MAX) + mn;
  }
  DenseLayer::DenseLayer(int input_size, int output_size)
  {
    input_shape.push_back(input_size);
    output_shape.push_back(output_size);
  }

  void DenseLayer::init(float *param_mem, float *gradient_mem, float *tmp_mem, bool initialize_random_weights)
  {
    A = TensorView(param_mem + 0, Shape{input_shape[0], output_shape[0]});
    b = TensorView(param_mem + A.total_size, Shape{output_shape[0]});

    dLoss_dA = TensorView(gradient_mem + 0, Shape{input_shape[0], output_shape[0]});
    dLoss_db = TensorView(gradient_mem + A.total_size, Shape{output_shape[0]});

    At = TensorView(tmp_mem + 0, Shape{output_shape[0], input_shape[0]});
    op = TensorView(tmp_mem + 0, Shape{input_shape[0], output_shape[0]}); // never used together

    if (initialize_random_weights)
    {
      SIREN_initialization(A, input_shape[0], output_shape[0]);
      SIREN_initialization(b, input_shape[0], output_shape[0]);
    }
  }

  void DenseLayer::forward(const TensorView &input, TensorView &output)
  {
    // y = A*x + b
    vec_mul(A, input, output);
    add(output, b);
  }

  void DenseLayer::backward(const TensorView &input, const TensorView &output,
                            const TensorView &dLoss_dOutput, TensorView dLoss_dInput, bool first_layer)
  {
    if (!first_layer)
    {
      // dLoss_dInput = A^T * dloss_dOutput
      transpose(A, At);
      vec_mul(At, dLoss_dOutput, dLoss_dInput);
    }

    // dLoss_dA += dLoss_dOutput ⊗ input
    // print_scheme(dLoss_dOutput.Dim, dLoss_dOutput.scheme);printf("\n");
    // print_scheme(input.Dim, input.scheme);printf("\n");
    // print_scheme(op.Dim, op.scheme);printf("\n");
    vector_outer_product(dLoss_dOutput, input, op);
    add(dLoss_dA, op);

    // dLoss_db += dLoss_dOutput;
    add(dLoss_db, dLoss_dOutput);
  }

  int DenseLayer::parameters_memory_size()
  {
    return (input_shape[0] + 1) * output_shape[0];
  }

  int DenseLayer::tmp_memory_size()
  {
    return input_shape[0] * output_shape[0];
  }
  class Optimizer
  {
  public:
    Optimizer(int _params_count)
    {
      params_count = _params_count;
    }
    virtual ~Optimizer(){};
    virtual void step(float *params_ptr, float const *grad_ptr) = 0;

  protected:
    int params_count = 0;
  };

  class OptimizerGD : public Optimizer
  {
  public:
    OptimizerGD(int _params_count, float _lr = 0.01f) : Optimizer(_params_count)
    {
      lr = _lr;
    }
    virtual void step(float *params_ptr, float const *grad_ptr) override
    {
      for (int i = 0; i < params_count; i++)
        params_ptr[i] -= lr * grad_ptr[i];
    }

  private:
    float lr;
  };

  class OptimizerAdam : public Optimizer
  {
  public:
    OptimizerAdam(int _params_count, float _lr = 0.01f, float _beta_1 = 0.9f, float _beta_2 = 0.999f, float _eps = 1e-8) : Optimizer(_params_count)
    {
      lr = _lr;
      beta_1 = _beta_1;
      beta_2 = _beta_2;
      eps = _eps;

      V = std::vector<float>(_params_count, 0);
      S = std::vector<float>(_params_count, 0);
    }
    virtual void step(float *params_ptr, float const *grad_ptr) override
    {
      for (int i = 0; i < params_count; i++)
      {
        float g = grad_ptr[i];
        V[i] = beta_1 * V[i] + (1 - beta_1) * g;
        float Vh = V[i] / (1 - pow(beta_1, iter + 1));
        S[i] = beta_2 * S[i] + (1 - beta_2) * g * g;
        float Sh = S[i] / (1 - pow(beta_2, iter + 1));
        params_ptr[i] -= lr * Vh / (sqrt(Sh) + eps);
      }
      iter++;
    }

  private:
    float lr, beta_1, beta_2, eps;
    int iter = 0;
    std::vector<float> V;
    std::vector<float> S;
  };

  float loss_MSE(const TensorView &values, const TensorView &target_values, TensorView dLoss_dValues)
  {
    int len = values.size(0);
    float loss = 0;
    if (dLoss_dValues.Dim > 0)
    {
      for (int i = 0; i < len; i++)
      {
        float l = values.get(i) - target_values.get(i);
        loss += l * l;
        dLoss_dValues.get(i) = 2 * l;
      }
    }
    else
    {
      for (int i = 0; i < len; i++)
      {
        float l = values.get(i) - target_values.get(i);
        loss += l * l;
      }
    }
    return loss;
  }

  float loss_cross_entropy(const TensorView &values, const TensorView &target_values, TensorView dLoss_dValues)
  {
    int len = values.size(0);
    float loss = 0;
    for (int i = 0; i < len; i++)
      loss += -target_values.get(i) * logf(values.get(i) + 1e-15f);
    if (dLoss_dValues.Dim > 0)
      for (int i = 0; i < len; i++)
        dLoss_dValues.get(i) = -target_values.get(i) / (values.get(i) + 1e-15f);

    return loss;
  }

  float score_MSE(const TensorView &values, const TensorView &target_values)
  {
    return loss_MSE(values, target_values, TensorView());
  }

  float score_accuracy(const TensorView &values, const TensorView &target_values)
  {
    // values are one-hot encoded;
    int count = values.size(1);
    int len = values.size(0);
    int predicted_label = 0;
    int true_label = 0;
    for (int j = 0; j < len; j++)
    {
      if (values.get(j) > values.get(predicted_label))
        predicted_label = j;
      if (target_values.get(j) > target_values.get(true_label))
        true_label = j;
    }
    return (predicted_label == true_label);
  }

  void NeuralNetwork::add_layer(std::shared_ptr<Layer> layer)
  {
    layers.push_back(layer);
  }

  bool NeuralNetwork::check_validity()
  {
    for (int i = 1; i < layers.size(); i++)
    {
      // shape-insensitive layers
      if (layers[i]->input_shape.empty() || layers[i - 1]->output_shape.empty())
      continue;

      if (layers[i]->input_shape.size() != layers[i - 1]->output_shape.size())
      {
        printf("NeuralNetwork: layers %d and %d have incompatible shapes!\n", i - 1, i);
        return false;
      }
      for (int j = 0; j < layers[i]->input_shape.size(); j++)
      {
        if (layers[i]->input_shape[j] != layers[i - 1]->output_shape[j])
        {
          printf("NeuralNetwork: layers %d and %d have incompatible sizes!\n", i - 1, i);
          return false;
        }
      }
    }
    return true;
  }

  void NeuralNetwork::initialize()
  {
    if (!check_validity())
      return;

    int total_params = 0;
    int layers_tmp_size = 1;
    int network_tmp_size = 1;
    for (int i = 0; i < layers.size(); i++)
    {
      total_params += layers[i]->parameters_memory_size();
      layers_tmp_size = std::max(layers_tmp_size, layers[i]->tmp_memory_size());
      network_tmp_size = std::max(network_tmp_size, (int)get_total_size(layers[i]->output_shape));
    }

    weights = std::vector<float>(total_params, 0);
    tmp_mem = std::vector<float>(layers_tmp_size + 2 * network_tmp_size, 0);

    // init layer for evaluation (null-pointing tensors for gradients)
    float *cur_param_ptr = weights.data();
    int li = 0;
    layer_outputs.clear();
    for (int i = 0; i < layers.size(); i++)
    {
      if (i > 0 && layers[i]->input_shape.empty())
      {
        layers[i]->input_shape = layers[i - 1]->output_shape;
        layers[i]->output_shape = layers[i - 1]->output_shape;
      }
      layers[i]->init(cur_param_ptr, nullptr, tmp_mem.data(), false);
      cur_param_ptr += layers[i]->parameters_memory_size();

      layer_outputs.push_back(TensorView(tmp_mem.data() + layers_tmp_size + li * network_tmp_size, layers[i]->output_shape));
      li = (li + 1) % 2;
    }

    printf("Neural Network succesfully created\n");
    printf("%d layers\n", (int)(layers.size()));
    for (int i = 0; i < layers.size(); i++)
      printf("Layer %d has %d parameters\n", i, layers[i]->parameters_memory_size());
    printf("%d input size\n", get_total_size(layers[0]->input_shape));
    printf("%d output size\n", get_total_size(layers.back()->output_shape));
    printf("%d weights\n", total_params);
  }

  void NeuralNetwork::initialize_with_weights(const float *init_weights)
  {
    initialize();
    weights = std::vector<float>(init_weights, init_weights + weights.size());
  }

  void NeuralNetwork::initialize_from_file(std::string filename)
  {
    initialize();
    std::ifstream in(filename, std::ios_base::binary);
    assert(in.is_open());
    in.read(reinterpret_cast<char*>(weights.data()), sizeof(float)*weights.size());
    in.close();
  }

  void NeuralNetwork::save_weights_to_file(std::string filename)
  {
    std::ofstream out(filename, std::ios_base::binary);
    assert(out.is_open());
    out.write(reinterpret_cast<char*>(weights.data()), sizeof(float)*weights.size());
    out.close();
  }

  void NeuralNetwork::train(const TensorView &inputs /*[input_size, count]*/, const TensorView &outputs /*[output_size, count]*/,
                            const TensorView &inputs_val, const TensorView &outputs_val,
                            int batch_size, int iterations, Opt opt, Loss loss_func, float lr)
  {
    // check if the input is correct
    assert(inputs.Dim == 2);
    assert(compact_dims(inputs.scheme) == 2);
    assert(outputs.Dim == 2);
    assert(compact_dims(outputs.scheme) == 2);
    assert(inputs.size(1) == outputs.size(1));

    // initialize network
    initialize();

    // initialize optimizer
    std::unique_ptr<Optimizer> optimizer;
    if (opt == Opt::GD)
      optimizer.reset(new OptimizerGD(weights.size(), lr));
    else if (opt == Opt::Adam)
      optimizer.reset(new OptimizerAdam(weights.size(), lr));

    // calculate memory requirements for training
    int io_size = 0;
    for (auto &l : layers)
      io_size += get_total_size(l->output_shape);

    // allocate memory needed for training
    std::vector<float> gradients(weights.size(), 0);
    std::vector<float> io_mem(io_size, 0);

    // create tensors to store layers' outputs
    std::vector<TensorView> layer_dLoss_dOutputs = layer_outputs;
    std::vector<TensorView> layer_saved_outputs;
    float *cur_io_ptr = io_mem.data();
    for (auto &l : layers)
    {
      layer_saved_outputs.push_back(TensorView(cur_io_ptr, l->output_shape));
      cur_io_ptr += get_total_size(l->output_shape);
    }

    // init layer for training
    float *cur_param_ptr = weights.data();
    float *cur_grad_ptr = gradients.data();
    for (auto &l : layers)
    {
      l->init(cur_param_ptr, cur_grad_ptr, tmp_mem.data(), true);
      cur_param_ptr += l->parameters_memory_size();
      cur_grad_ptr += l->parameters_memory_size();
    }

    // main training loop
    int count = inputs.size(inputs.Dim - 1);
    std::chrono::steady_clock::time_point t1, t2;
    t1 = std::chrono::steady_clock::now();

    for (int iter = 0; iter < iterations; iter++)
    {
      float loss = 0;
      std::fill_n(gradients.data(), gradients.size(), 0);

      for (int batch_id = 0; batch_id < batch_size; batch_id++)
      {
        int sample_id = rand() % count;
        TensorView input_batch = slice(inputs, sample_id);
        TensorView target_output_batch = slice(outputs, sample_id);

        // forward pass
        layers[0]->forward(input_batch, layer_saved_outputs[0]);
        for (int i = 1; i < layers.size(); i++)
          layers[i]->forward(layer_saved_outputs[i - 1], layer_saved_outputs[i]);

        loss += calculate_loss(loss_func, layer_saved_outputs.back(), target_output_batch, layer_dLoss_dOutputs.back());

        // backward pass
        for (int i = layers.size() - 1; i >= 1; i--)
          layers[i]->backward(layer_saved_outputs[i - 1], layer_saved_outputs[i], layer_dLoss_dOutputs[i], layer_dLoss_dOutputs[i - 1], false);
        layers[0]->backward(input_batch, layer_saved_outputs[0], layer_dLoss_dOutputs[0], TensorView(), true);
      }

      loss /= batch_size;
      for (auto &g : gradients)
        g /= batch_size;
      optimizer->step(weights.data(), gradients.data());

      if (iter % 100 == 0 && iter > 0)
      {
        t2 = std::chrono::steady_clock::now();
        float time_ms = 1e-6*std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        t1 = t2;
        float pass_per_msec = 100*batch_size/time_ms;
        printf("iteration %d/%d: average loss %f. %.1f Kpass/sec \n", iter, iterations, loss, pass_per_msec);
        if (inputs_val.Dim > 0)
          test(inputs_val, outputs_val, loss_func);
      }
    }

    // printf("weights [ ");
    // for (auto &w: weights)
    //   printf("%f ", w);
    // printf("]\n");
  }

  void NeuralNetwork::evaluate(const TensorView &input, TensorView output)
  {
    int elements = input.size(input.Dim - 1);
    for (int i = 0; i < elements; i++)
    {
      // forward pass
      layers[0]->forward(slice(input, i), layer_outputs[0]);
      for (int j = 1; j < layers.size(); j++)
        layers[j]->forward(layer_outputs[j - 1], layer_outputs[j]);
      copy(layer_outputs.back(), slice(output, i));
    }
  }

  float NeuralNetwork::test(const TensorView &input, const TensorView &target_output, Loss loss_func)
  {
    int elements = input.size(input.Dim - 1);
    float loss = 0;
    float score = 0;
    for (int i = 0; i < elements; i++)
    {
      // forward pass
      layers[0]->forward(slice(input, i), layer_outputs[0]);
      for (int j = 1; j < layers.size(); j++)
        layers[j]->forward(layer_outputs[j - 1], layer_outputs[j]);

      loss += calculate_loss(loss_func, layer_outputs.back(), slice(target_output, i), TensorView());
      score += calculate_score(loss_func, layer_outputs.back(), slice(target_output, i));
    }
    loss /= elements;
    score /= elements;
    printf("Validation loss %f score %f\n", loss, score);
    return loss;
  }

  float NeuralNetwork::calculate_loss(Loss loss, const TensorView &values, const TensorView &target_values, TensorView dLoss_dValues)
  {
    switch (loss)
    {
    case Loss::MSE:
      return loss_MSE(values, target_values, dLoss_dValues);
      break;
    case Loss::CrossEntropy:
      return loss_cross_entropy(values, target_values, dLoss_dValues);
      break;
    default:
      return 0;
      break;
    }
  }

  float NeuralNetwork::calculate_score(Loss loss, const TensorView &values, const TensorView &target_values)
  {
    switch (loss)
    {
    case Loss::MSE:
      return score_MSE(values, target_values);
      break;
    case Loss::CrossEntropy:
      return score_accuracy(values, target_values);
      break;
    default:
      return 0;
      break;
    }
  }
}