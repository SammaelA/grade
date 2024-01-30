#pragma once
#include <vector>
#include <string>

namespace nn
{
  struct Dataset
  {
    std::vector<float> train_data;
    std::vector<float> train_labels;
    std::vector<float> test_data;
    std::vector<float> test_labels;

    std::vector<unsigned> element_size;
    unsigned label_size = 1;
    unsigned train_elements = 0;
    unsigned test_elements = 0;
  };

  void read_CIFAR10_dataset(std::string path, Dataset *out_dataset);
}