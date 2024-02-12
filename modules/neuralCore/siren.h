#pragma once
#include "neural_network.h"

namespace nn
{
  class Siren : public NeuralNetwork
  {
  public:
    enum class Type
    {
      Image,
      SDF,
      Gen_SDF_32124,
      Chair
    };

    Siren(Type type, int hidden_layers, int layer_size);
    void train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
               int batch_size, int iterations, bool verbose = false);
    float get(float x, float y);
    float get(float x, float y, float z);
  private:
    std::vector<float> point;
    std::vector<float> distance;
  };
}