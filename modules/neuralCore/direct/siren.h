#pragma once
#include "neural_network.h"

namespace nnd
{
  class Siren : public NeuralNetwork
  {
  public:
    enum class Type
    {
      Image,
      SDF
    };

    Siren(Type type, int hidden_layers, int layer_size);
    void train(const TensorView &inputs /*[input_size, count]*/, const TensorView &outputs /*[output_size, count]*/,
               int batch_size, int iterations);
    float get(float x, float y);
    float get(float x, float y, float z);
  private:
    std::array<float, 3> one_point_tv_data;
    float one_point_tv_res;
    TensorView one_point_tv;
    TensorView res_tv;
  };
}