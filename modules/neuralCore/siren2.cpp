#include "siren2.h"
#include <cassert>

namespace nn
{
  Siren2::Siren2(Type type, int hidden_layers, int layer_size)
  {
    assert(hidden_layers > 0);
    assert(layer_size > 3);
    unsigned input_dim = (type==Type::Image) ? 2 : 3;
    point.resize(input_dim);
    distance.resize(1);

    add_layer(std::make_shared<DenseLayer2>(input_dim, layer_size), WeightsInitializer::SIREN);
    for (int i=0;i<hidden_layers-1;i++)
    {
      add_layer(std::make_shared<SinLayer2>());
      add_layer(std::make_shared<DenseLayer2>(layer_size, layer_size), WeightsInitializer::SIREN);
    }
    add_layer(std::make_shared<SinLayer2>());
    add_layer(std::make_shared<DenseLayer2>(layer_size, 1), WeightsInitializer::SIREN);
  }

  void Siren2::train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
                    int batch_size, int iterations)
  {
    NeuralNetwork2::train(inputs, outputs, batch_size, iterations, 
                          NeuralNetwork2::Opt::Adam, NeuralNetwork2::Loss::MSE, 0.0001f);
  }

  float Siren2::get(float x, float y)
  {
    point[0] = x;
    point[1] = y;

    evaluate(point, distance);

    return distance[0];
  }

  float Siren2::get(float x, float y, float z)
  {
    point[0] = x;
    point[1] = y;
    point[2] = z;
    
    evaluate(point, distance);

    return distance[0];
  }
}
