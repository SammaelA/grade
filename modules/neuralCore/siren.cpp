#include "siren.h"
#include <cassert>

namespace nn
{
  Siren::Siren(Type type, int hidden_layers, int layer_size)
  {
    assert(hidden_layers > 0);
    assert(layer_size > 3);
    unsigned input_dim = 3;
    switch (type)
    {
      case Type::Image:
        input_dim = 2;
        break;
      case Type::SDF:
        input_dim = 3;
        break;
      case Type::Gen_SDF_32124:
        input_dim = 3+3+1+3+3;
        break;
      case Type::Chair:
        input_dim = 3+6;
        break;
    }
    point.resize(input_dim);
    distance.resize(1);

    add_layer(std::make_shared<DenseLayer>(input_dim, layer_size), Initializer::Siren);
    for (int i=0;i<hidden_layers-1;i++)
    {
      add_layer(std::make_shared<SinLayer>());
      add_layer(std::make_shared<DenseLayer>(layer_size, layer_size), Initializer::Siren);
    }
    add_layer(std::make_shared<SinLayer>());
    add_layer(std::make_shared<DenseLayer>(layer_size, 1), Initializer::Siren);
  }

  void Siren::train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
                    int batch_size, int iterations, bool verbose)
  {
    NeuralNetwork::train(inputs, outputs, batch_size, iterations, OptimizerAdam(0.00002f), Loss::MSE, verbose);
  }

  float Siren::get(float x, float y)
  {
    point[0] = x;
    point[1] = y;

    evaluate(point, distance);

    return distance[0];
  }

  float Siren::get(float x, float y, float z)
  {
    point[0] = x;
    point[1] = y;
    point[2] = z;
    
    evaluate(point, distance);

    return distance[0];
  }
}
