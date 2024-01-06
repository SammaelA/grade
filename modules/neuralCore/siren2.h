#pragma once
#include "neural_network_2.h"

namespace nn
{
  class Siren2 : public NeuralNetwork2
  {
  public:
    enum class Type
    {
      Image,
      SDF
    };

    Siren2(Type type, int hidden_layers, int layer_size);
    void train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
               int batch_size, int iterations);
    float get(float x, float y);
    float get(float x, float y, float z);
  private:
    std::vector<float> point;
    std::vector<float> distance;
  };
}