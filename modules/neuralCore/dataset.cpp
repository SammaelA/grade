#include "dataset.h"
#include <fstream>
#include <cassert>

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
      (*labels)[offset + i*classes + label] = 1;
      //printf("%d label %u\n",i,label);

      for (int j=0;j<image_size;j++)
        (*data)[offset + i*image_size + j] = raw_data[i*(image_size+1) + j + 1]/255.0f;
    }
  }
  void read_CIFAR10_dataset(std::string path, Dataset *out_dataset)
  {
    //CIFAR10 dataset, files taken from https://www.cs.toronto.edu/~kriz/cifar.html 
    out_dataset->element_size = {32, 32, 3}; //32x32 RGB images
    out_dataset->label_size = 10;
    out_dataset->train_elements = 50000;
    out_dataset->test_elements  = 10000;

    CIFAR10_add_batch(path + "/data_batch_1.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_2.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_3.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_4.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/data_batch_5.bin", &(out_dataset->train_data), &(out_dataset->train_labels));
    CIFAR10_add_batch(path + "/test_batch.bin"  , &(out_dataset->test_data ), &(out_dataset->test_labels ));
  } 
}
