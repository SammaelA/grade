#include "neural_network.h"
#include <cassert>
#include <cstdio>
#include <string>
#include <fstream>
#include <chrono>
#include <cstring>
#include <random>

constexpr bool DEBUG = false;
namespace nn
{
  void DenseLayer::init()
  {
    weights.clear();

    weights.push_back(TensorToken(output_shape[0], input_shape[0])); // A
    weights.push_back(TensorToken(output_shape[0]));                 // b

    dLoss_dWeights.resize(weights.size());
  }

  TensorToken DenseLayer::forward(const TensorToken &in)
  {
    return TensorToken::mat_mul_t(in, weights[0].transpose()) + weights[1];
  }

  TensorToken DenseLayer::backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput)
  {
    float batch_size = (float)(input.sizes[1]);

    dLoss_dWeights[0] = TensorToken::mat_mul_t(input.transpose(), dLoss_dOutput.transpose())/batch_size;
    dLoss_dWeights[1] = dLoss_dOutput.outer_sum()/batch_size;

    return TensorToken::mat_mul_t(dLoss_dOutput, weights[0]);
  }

  unsigned total_size(const std::vector<unsigned> &sizes)
  {
    unsigned total_sz = 1;
    for (auto sz : sizes)
      total_sz *= sz;
    return total_sz;
  }

  void zero_initialization(float *data, int size)
  {
    std::fill_n(data, size, 0.0f);
  }

  void he_initialization(float *data, int size, int fan_in, int fan_out)
  {
    float mn = -sqrt(6.0 / fan_in);
    float mx = sqrt(6.0 / fan_in);
    float d = mx - mn;
    for (int i = 0; i < size; i++)
      data[i] = d * (((double)rand()) / RAND_MAX) + mn;
  }
  void SIREN_initialization(float *data, int size, int fan_in, int fan_out)
  {
    constexpr float omega_0 = 30.0f;
    float mx = (fan_in == 2) ? 0.5 : (sqrt(6.0 / fan_in) / omega_0);
    float mn = -mx;
    float d = mx - mn;
    for (int i = 0; i < size; i++)
      data[i] = d * (((double)rand()) / RAND_MAX) + mn;
  }

  void Glorot_normal_initialization(float *data, int size, int fan_in, int fan_out)
  {
    double stddev = sqrt(2.0f/(fan_in+fan_out));
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution distr{0.0, stddev};
    for (int i = 0; i < size; i++)
      data[i] = distr(gen);
  }

  void NeuralNetwork::add_layer(std::shared_ptr<Layer> layer, WeightsInitializer initializer)
  {
    layers.push_back(layer);
    initializers.push_back(initializer);
  }

  void NeuralNetwork::set_batch_size_for_evaluate(int size)
  {
    batch_size_evaluate = size;
    if (initialized)
      get_evaluate_prog();
  }

  bool NeuralNetwork::check_validity()
  {
    if (layers[0]->input_shape.empty() || layers[0]->output_shape.empty())
    {
      printf("NeuralNetwork: first layer must have implicit shape!\n");
      return false;
    }
    for (int i = 1; i < layers.size(); i++)
    {
      // shape-insensitive layers
      if (layers[i]->input_shape.empty() && layers[i]->output_shape.empty())
      {
        layers[i]->input_shape = layers[i-1]->output_shape;
        layers[i]->output_shape = layers[i-1]->output_shape;
      }
      else if (layers[i]->input_shape.size() != layers[i - 1]->output_shape.size())
      {
        printf("NeuralNetwork: layers %d and %d have incompatible shapes!\n", i - 1, i);
        return false;
      }
    }
    return true;
  }

  void NeuralNetwork::initialize()
  {
    if (!check_validity())
      return;
    
    for (auto &l : layers)
      total_params += l->parameters_count();
    weights.resize(total_params, 0.0f);

    int offset = 0;
    for (int i=0;i<layers.size();i++)
    {
      unsigned fan_in = total_size(layers[i]->input_shape);
      unsigned fan_out = total_size(layers[i]->output_shape);
      unsigned size = layers[i]->parameters_count();

      if (size == 0)
        continue;

      switch (initializers[i])
      {
        case ZERO:
          zero_initialization(weights.data()+offset, size);
          break;
        case HE:
          he_initialization(weights.data()+offset, size, fan_in, fan_out);
          break;
        case SIREN:
          SIREN_initialization(weights.data()+offset, size, fan_in, fan_out);
          break;
        case GLOROT_NORMAL:
          Glorot_normal_initialization(weights.data()+offset, size, fan_in, fan_out);
          break;
        default:
          break;
      }
      offset += size;
    }

    get_evaluate_prog();
    //print_info();
    initialized = true;
  }

  void NeuralNetwork::initialize_with_weights(const float *w)
  {
    initialize();
    weights = std::vector<float>(w, w + total_params);
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

  void NeuralNetwork::set_arch_to_file(std::string filename)
  {
    std::ofstream out(filename);
    assert(out.is_open());
    for (auto &l : layers)
    {
      out << l->get_name() << " ";
      if (l->input_shape.size() > 0)
      {
        out << "input shape ("<<l->input_shape[0];
        for (int i=1;i<l->input_shape.size();i++)
          out << ", "<<l->input_shape[i];
        out <<") ";
      }
      if (l->output_shape.size() > 0)
      {
        out << "output shape ("<<l->output_shape[0];
        for (int i=1;i<l->output_shape.size();i++)
          out << ", "<<l->output_shape[i];
        out <<")";
      }
      out<<"\n";
    }
    out.close();
  }

  void NeuralNetwork::print_info()
  {
    printf("Neural Network\n");
    printf("%d layers\n", (int)(layers.size()));
    for (int i = 0; i < layers.size(); i++)
      printf("Layer %d has %d parameters\n", i, layers[i]->parameters_count());
    printf("%d input size\n", total_size(layers[0]->input_shape));
    printf("%d output size\n", total_size(layers.back()->output_shape));
    printf("%d weights\n", total_params);
  }

  void NeuralNetwork::get_evaluate_prog()
  {
    TensorCompiler compiler;
    compiler.start_program();
    auto i_shape = layers[0]->input_shape;
    i_shape.push_back(batch_size_evaluate);
    auto o_shape = layers.back()->output_shape;
    o_shape.push_back(batch_size_evaluate);

    TensorToken input = TensorToken(i_shape);
    TensorToken output = TensorToken(o_shape);
    TensorToken w = TensorToken(total_params);
    unsigned offset = 0;
    for (auto &l : layers)
    {
      l->init();
      for (auto &lw : l->weights)
      {
        unsigned sz = lw.total_size();
        lw.copy_to({0, sz}, w, {offset, offset + sz});
        offset += sz;
      }
    }

    TensorToken t = input;
    for (auto &l : layers)
      t = l->forward(t);
    output = t;
    
    compiler.inout(w, "W");
    compiler.input(input, "In");
    compiler.output(output, "Out");
    evaluate_prog = compiler.finish_program();
  }

  TensorProgram NeuralNetwork::get_train_prog(int batch_size, Opt optimizer, Loss loss, float lr)
  {
    TensorCompiler compiler;
    compiler.start_program();
    auto i_shape = layers[0]->input_shape;
    i_shape.push_back(batch_size);
    auto o_shape = layers.back()->output_shape;
    o_shape.push_back(batch_size);

    TensorToken input = TensorToken(i_shape); compiler.input(input, "In");
    TensorToken target_output = TensorToken(o_shape); compiler.input(target_output, "Out");
    TensorToken w = TensorToken(total_params); compiler.inout(w, "W");
    
    unsigned offset = 0;
    for (auto &l : layers)
    {
      l->init();
      for (int i=0;i<l->weights.size();i++)
      {
        unsigned sz = l->weights[i].total_size();
        l->weights[i].copy_to({0, sz}, w, {offset, offset + sz});
        offset += sz;
      }
    }

    std::vector<TensorToken> all_outputs;
    all_outputs.push_back(layers[0]->forward(input));
        TensorToken t = input;
    for (int i=1;i<layers.size();i++)
      all_outputs.push_back(layers[i]->forward(all_outputs.back()));
  
    TensorToken output = all_outputs.back();

    //loss
    TensorToken l, dLoss_dOutput;
    if (loss == Loss::MSE)
    {
      TensorToken diff = output - target_output;
      l = (diff*diff).sum()/(float)(output.total_size());
      dLoss_dOutput = 2.0f*diff;
    }
    else if (loss == Loss::CrossEntropy)
    {
      TensorToken mo = -1.0f*target_output;
      l = (mo * TensorToken::log(output + 1e-15f)).sum()/(float)(batch_size);
      dLoss_dOutput = mo / (output + 1e-15f);
    }

    for (int i=layers.size()-1;i>0;i--)
      dLoss_dOutput = layers[i]->backward(all_outputs[i-1], all_outputs[i], dLoss_dOutput);
    layers[0]->backward(input, all_outputs[0], dLoss_dOutput);
  
    TensorToken grad = TensorToken(total_params);
    offset = 0;
    for (auto &l : layers)
    {
      for (auto &dLoss_dWeight : l->dLoss_dWeights)
      {
        unsigned sz = dLoss_dWeight.total_size();
        grad.set({offset, offset + sz}, dLoss_dWeight.flatten());
        offset += sz;
      }
    }

    //Adam optimizer
    TensorToken V = TensorToken(total_params); compiler.inout(V, "V");
    TensorToken S = TensorToken(total_params); compiler.inout(S, "S");
    TensorToken iter = TensorToken(1); compiler.input(iter, "iter");
    float beta_1 = 0.9f;
    float beta_2 = 0.999f;
    float eps = 1e-7f;

    V = beta_1*V + (1.0f - beta_1)*grad;
    TensorToken Vh = V / (1.0f - TensorToken::pow(beta_1, iter + 1.0f));
    S = beta_2*S + (1.0f - beta_2)*grad*grad;
    TensorToken Sh = S / (1.0f - TensorToken::pow(beta_2, iter + 1.0f));
    w -= lr*Vh/(TensorToken::sqrt(Sh) + eps);
    
    compiler.output(l, "loss");
    if (DEBUG)
      compiler.output(grad, "grad");
    return compiler.finish_program();
  }

  void NeuralNetwork::evaluate(std::vector<float> &input_data, std::vector<float> &output_data, int samples)
  {
    unsigned input_size = total_size(layers[0]->input_shape);
    unsigned output_size = total_size(layers.back()->output_shape);

    if (samples < 0)
      samples = input_data.size()/input_size;

    unsigned batches = (samples + batch_size_evaluate - 1)/batch_size_evaluate;
    
    TensorProcessor::set_program(evaluate_prog);
    TensorProcessor::set_input("W", weights.data(), weights.size());
    
    for (int i=0;i<batches;i++)
    {
      TensorProcessor::set_input("In", input_data.data() + i*batch_size_evaluate*input_size, samples*input_size - i*batch_size_evaluate*input_size);
      TensorProcessor::execute();
      TensorProcessor::get_output("Out", output_data.data() + i*batch_size_evaluate*output_size, samples*output_size - i*batch_size_evaluate*output_size);
    }
  }

  void NeuralNetwork::train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
                             int batch_size, int iterations, Opt optimizer, Loss loss, float lr, bool verbose)
  {
    initialize();

    TensorProgram train_prog = get_train_prog(batch_size, optimizer, loss, lr);

    unsigned input_size = total_size(layers[0]->input_shape);
    unsigned output_size = total_size(layers.back()->output_shape);
    unsigned count = inputs.size()/input_size;
    assert(inputs.size() % input_size == 0);
    assert(outputs.size() % output_size == 0);
    assert(count == outputs.size() / output_size);

    std::vector<float> V(total_params, 0);
    std::vector<float> S(total_params, 0);
    std::vector<float> in_batch(input_size*batch_size);
    std::vector<float> out_batch(output_size*batch_size);
    float iter = 0;
    TensorProcessor::set_program(train_prog);
    TensorProcessor::set_input("W", weights.data(), weights.size());
    TensorProcessor::set_input("V", V.data(), V.size());
    TensorProcessor::set_input("S", S.data(), S.size());

    float av_loss = 0;
    for (int it=0;it<iterations;it++)
    {
      for (int i=0;i<batch_size;i++)
      {
        unsigned b_id = rand()%count;
        memcpy(in_batch.data() + i*input_size, inputs.data() + b_id*input_size, sizeof(float)*input_size);
        memcpy(out_batch.data() + i*output_size, outputs.data() + b_id*output_size, sizeof(float)*output_size);
      }

      iter = it;
      TensorProcessor::set_input("In", in_batch.data(), in_batch.size());
      TensorProcessor::set_input("Out", out_batch.data(), out_batch.size());
      TensorProcessor::set_input("iter", &iter, 1);
      TensorProcessor::execute();
      float loss = -1;
      TensorProcessor::get_output("loss", &loss, 1);
      av_loss += loss;
      if (verbose && it % 100 == 0)
      {
        printf("[%d/%d] Loss = %f %f\n", it, iterations, loss, av_loss/100);
        av_loss = 0;
      }

      if (DEBUG)
      {
        std::vector<float> grad(weights.size(),0);
        TensorProcessor::get_output("grad", grad.data(), grad.size());
        TensorProcessor::get_output("W", weights.data(), weights.size());
        printf("grad = [ ");
        for (int i=0;i<weights.size();i++)
          printf("%f ", grad[i]);
        printf("]\n");

        printf("w = [ ");
        for (int i=0;i<weights.size();i++)
          printf("%f ", weights[i]);
        printf("]\n");
      }
    }
    TensorProcessor::get_output("W", weights.data(), weights.size());
  }

}