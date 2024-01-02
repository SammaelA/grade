#include "neural_network_2.h"
#include <cassert>
#include <cstdio>
#include <string>
#include <fstream>
#include <chrono>

namespace nn
{
  void DenseLayer2::init()
  {
    weights.push_back(TensorToken(input_shape[0], output_shape[0])); // A
    weights.push_back(TensorToken(output_shape[0]));                 // b

    dLoss_dWeights = weights;
  }

  TensorToken DenseLayer2::forward(const TensorToken &in)
  {
    return TensorToken::mat_vec_mul(weights[0], in) + weights[1];
  }

  TensorToken DenseLayer2::backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput)
  {
    TensorToken At = weights[0].transpose();
    TensorToken dLoss_dInput = TensorToken::mat_vec_mul(At, dLoss_dOutput);

    dLoss_dWeights[0] = TensorToken::vector_outer_product(dLoss_dOutput, input);
    dLoss_dWeights[1] = dLoss_dOutput;

    return dLoss_dInput;
  }
  void NeuralNetwork2::add_layer(std::shared_ptr<Layer2> layer)
  {
    layers.push_back(layer);
  }

  bool NeuralNetwork2::check_validity()
  {
    for (int i = 1; i < layers.size(); i++)
    {
      // shape-insensitive layers
      if (layers[i]->input_shape.empty() || layers[i - 1]->output_shape.empty())
        ;
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

  void NeuralNetwork2::initialize()
  {
    if (!check_validity())
      return;
    
    for (auto &l : layers)
      total_params += l->parameters_count();
    weights.resize(total_params);
    get_evaluate_prog();
    print_info();
  }

  void NeuralNetwork2::initialize_with_weights(const float *w)
  {
    initialize();
    weights = std::vector<float>(w, w + total_params);
  }

  void NeuralNetwork2::initialize_from_file(std::string filename)
  {
    initialize();
    std::ifstream in(filename, std::ios_base::binary);
    assert(in.is_open());
    in.read(reinterpret_cast<char*>(weights.data()), sizeof(float)*weights.size());
    in.close();
  }

  void NeuralNetwork2::save_weights_to_file(std::string filename)
  {
    std::ofstream out(filename, std::ios_base::binary);
    assert(out.is_open());
    out.write(reinterpret_cast<char*>(weights.data()), sizeof(float)*weights.size());
    out.close();
  }

  unsigned get_total_size_2(const std::vector<unsigned> &sizes)
  {
    unsigned total_sz = 1;
    for (auto sz : sizes)
      total_sz *= sz;
    return total_sz;
  }

  void NeuralNetwork2::print_info()
  {
    printf("Neural Network\n");
    printf("%d layers\n", (int)(layers.size()));
    for (int i = 0; i < layers.size(); i++)
      printf("Layer %d has %d parameters\n", i, layers[i]->parameters_count());
    printf("%d input size\n", get_total_size_2(layers[0]->input_shape));
    printf("%d output size\n", get_total_size_2(layers.back()->output_shape));
    printf("%d weights\n", total_params);
  }

  void NeuralNetwork2::get_evaluate_prog()
  {
    TensorCompiler compiler;
    compiler.start_program();

    TensorToken t = TensorToken(layers[0]->input_shape);
    TensorToken w = TensorToken(total_params);
    unsigned offset = 0;
    for (auto &l : layers)
    {
      l->init();
      for (auto &lw : l->weights)
      {
        unsigned sz = lw.total_size();
        lw.set_flatten({0, sz}, w.get({offset, offset + sz}));
        offset += sz;
      }
    }
    for (auto &l : layers)
      t = l->forward(t);
    
    compiler.input(w, "W");
    compiler.input(t, "In");
    compiler.output(t, "Out");
    evaluate_prog = compiler.finish_program();
  }

  void NeuralNetwork2::evaluate(std::vector<float> &input_data, std::vector<float> &output_data)
  {
    int cnt = input_data.size()/get_total_size_2(layers[0]->input_shape);
    for (int i=0;i<cnt;i++)
    {
      tp.set_program(evaluate_prog);
      tp.set_input({{"W", weights.data()}, {"In", input_data.data() + i*get_total_size_2(layers[0]->input_shape)}});
      tp.execute();
      tp.set_output({{"Out", output_data.data() + i*get_total_size_2(layers.back()->output_shape)}});
    }
  }

}