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

  //save and load already prepared dataset to/from binary file
  void save_dataset(std::string path, const Dataset *dataset);
  void load_dataset(std::string path, Dataset *out_dataset);
  void train_test_split(Dataset *dataset, float test_fraction = 0.1);

  //labels are trated as one-hot encoded class marks
  //this function duplicates entries of rare classes
  //and checks if all classes exist in dataset
  void rebalance_classes(Dataset *out_dataset);

  //read datasets from their specific binary formats
  void read_CIFAR10_dataset(std::string path, Dataset *out_dataset);
  void read_MNIST_dataset(std::string path, Dataset *out_dataset);
}