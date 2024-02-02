#include "dataset.h"
#include <fstream>
#include <cassert>
#include <random>
#include <algorithm>

namespace nn
{
  void CIFAR10_add_batch(std::string filename, std::vector<float> *data, std::vector<float> *labels)
  {
    unsigned images = 10000;
    unsigned image_size = 32*32*3;
    unsigned classes = 10;
    
    std::vector<unsigned char> raw_data(images*(image_size+1), 0);
    std::ifstream in(filename, std::ios_base::binary);
    assert(in.is_open());
    in.read(reinterpret_cast<char*>(raw_data.data()), raw_data.size());
    in.close();

    unsigned offset = data->size()/image_size;
    data->resize(data->size() + image_size*images);
    labels->resize(labels->size() + classes*images, 0.0f);

    //<1 x label><3072 x pixel>
    for (int i=0;i<images;i++)
    {
      unsigned label = raw_data[i*(image_size+1)];
      assert(label < classes);
      (*labels)[classes*offset + i*classes + label] = 1;
      //printf("%d label %u\n",i,label);

      for (int j=0;j<image_size;j++)
        (*data)[image_size*offset + i*image_size + j] = raw_data[i*(image_size+1) + j + 1]/255.0f;
    }
  }
  void read_CIFAR10_dataset(std::string path, Dataset *out_dataset)
  {
    //CIFAR10 dataset, files taken from https://www.cs.toronto.edu/~kriz/cifar.html 
    out_dataset->element_size = {32, 32, 3}; //32x32 RGB images
    out_dataset->label_size = 10;
    out_dataset->train_elements = 60000;

    CIFAR10_add_batch(path + "/data_batch_1.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_2.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_3.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_4.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_5.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/test_batch.bin"  , &(out_dataset->train_data), &(out_dataset->train_labels));
  } 

  void MNIST_add_batch(std::string filename, int label, std::vector<float> *data, std::vector<float> *labels)
  {
    unsigned images = 1000;
    unsigned image_size = 28*28;
    unsigned classes = 10;
    
    std::vector<unsigned char> raw_data(images*image_size, 0);
    std::ifstream in(filename, std::ios_base::binary);
    assert(in.is_open());
    in.read(reinterpret_cast<char*>(raw_data.data()), raw_data.size());
    in.close();

    unsigned offset = data->size()/image_size;
    data->resize(data->size() + image_size*images);
    labels->resize(labels->size() + classes*images, 0.0f);

    //<28 x 28 pixels>
    for (int i=0;i<images;i++)
    {
      (*labels)[classes*offset + i*classes + label] = 1;

      for (int j=0;j<image_size;j++)
        (*data)[image_size*offset + i*image_size + j] = raw_data[i*image_size + j]/255.0f;
    }
  }
  void read_MNIST_dataset(std::string path, Dataset *out_dataset)
  {
    //CIFAR10 dataset, files taken from https://www.cs.toronto.edu/~kriz/cifar.html 
    out_dataset->element_size = {28, 28, 1}; //28x28 images
    out_dataset->label_size = 10;
    out_dataset->train_elements = 10000;

    for (int i=0;i<10;i++)
      MNIST_add_batch(path + "/data"+std::to_string(i)+".bin", i, &(out_dataset->train_data), &(out_dataset->train_labels));
  }

  void train_test_split(Dataset *dataset, float test_fraction)
  {
    assert(test_fraction > 1.0f/dataset->train_elements && test_fraction < 1 - 1.0f/dataset->train_elements);
    assert(dataset->train_data.size() > 0 && dataset->test_data.empty());
    std::vector<int> indices(dataset->train_elements, 0);
    for (int i=0;i<dataset->train_elements; i++)
      indices[i] = i;
    std::shuffle(indices.begin(),indices.end(), std::default_random_engine{});

    unsigned element_size = dataset->train_data.size()/dataset->train_elements;
    unsigned test_elems = test_fraction*dataset->train_elements;
    unsigned train_elems = dataset->train_elements - test_elems;

    std::vector<float> train_data(element_size*train_elems);
    std::vector<float> train_labels(dataset->label_size*train_elems);
    std::vector<float> test_data(element_size*test_elems);
    std::vector<float> test_labels(dataset->label_size*test_elems);

    for (int i=0;i<train_elems;i++)
    {
      unsigned idx = indices[i];
      for (int j=0;j<element_size;j++)
        train_data[i*element_size + j] = dataset->train_data[idx*element_size + j];
      for (int j=0;j<dataset->label_size;j++)
        train_labels[i*dataset->label_size + j] = dataset->train_labels[idx*dataset->label_size + j];
    }
    for (int i=train_elems;i<dataset->train_elements;i++)
    {
      unsigned idx = indices[i];
      for (int j=0;j<element_size;j++)
        test_data[(i-train_elems)*element_size + j] = dataset->train_data[idx*element_size + j];
      for (int j=0;j<dataset->label_size;j++)
        test_labels[(i-train_elems)*dataset->label_size + j] = dataset->train_labels[idx*dataset->label_size + j];
    }

    dataset->train_elements = train_elems;
    dataset->train_data = train_data;
    dataset->train_labels = train_labels;

    dataset->test_elements = test_elems;
    dataset->test_data = test_data;
    dataset->test_labels = test_labels;
  }
}
