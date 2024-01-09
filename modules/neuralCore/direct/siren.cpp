#include "siren.h"
#include <cassert>

namespace nnd
{
  Siren::Siren(Type type, int hidden_layers, int layer_size)
  {
    assert(hidden_layers > 0);
    assert(layer_size > 3);
    IndexType input_dim = (type==Type::Image) ? 2 : 3;

    one_point_tv = TensorView(one_point_tv_data.data(), Shape{input_dim, 1});
    res_tv = TensorView(&one_point_tv_res, Shape{1, 1});

    add_layer(std::make_shared<DenseLayer>(input_dim, layer_size));
    for (int i=0;i<hidden_layers-1;i++)
    {
      add_layer(std::make_shared<SinLayer>());
      add_layer(std::make_shared<DenseLayer>(layer_size, layer_size));
    }
    add_layer(std::make_shared<SinLayer>());
    add_layer(std::make_shared<DenseLayer>(layer_size, 1));
  }

  void Siren::train(const TensorView &inputs /*[input_size, count]*/, const TensorView &outputs /*[output_size, count]*/,
               int batch_size, int iterations)
  {
    NeuralNetwork::train(inputs, outputs, TensorView(), TensorView(), batch_size, iterations, 
                         NeuralNetwork::Opt::Adam, NeuralNetwork::Loss::MSE, 0.0001);
  }

  float Siren::get(float x, float y)
  {
    one_point_tv.get(0, 0) = x;
    one_point_tv.get(1, 0) = y;
    
    evaluate(one_point_tv, res_tv);

    return res_tv.get(0,0);
  }

  float Siren::get(float x, float y, float z)
  {
    one_point_tv.get(0, 0) = x;
    one_point_tv.get(1, 0) = y;
    one_point_tv.get(2, 0) = z;
    
    evaluate(one_point_tv, res_tv);

    return res_tv.get(0,0);
  }
}
