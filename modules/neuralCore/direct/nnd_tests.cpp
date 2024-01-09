#include <cstdio>
#include <memory>
#include <algorithm>
#include <cstdint>
#include <cassert>
#include <vector>
#include <functional>
#include <array>
#include <type_traits>
#include <string>
#include <cstring>
#include <chrono>
#include <cmath>

#include "tensors.h"
#include "neural_network.h"
#include "siren.h"
#include "nnd_tests.h"

#include "../stb_image.h"
#include "../stb_image_write.h"


namespace nnd
{
  void write_image_rgb(std::string path, const TensorView &image_tensor)
  {
    assert(image_tensor.Dim == 3);
    assert(image_tensor.size(0) == 3);

    unsigned char *data = new unsigned char[image_tensor.total_size];
    for (int i=0;i<image_tensor.total_size;i++)
    data[i] = std::max(0, std::min(255, (int)(255*image_tensor.get(i))));
    
    stbi_write_png(path.c_str(), image_tensor.size(1), image_tensor.size(2), 3, data, 3*image_tensor.size(1));

    delete[] data;
  }

  TensorView read_image_rgb(std::string path, std::vector<float> &image_data)
  {
    int width, height, channels;
    unsigned char *imgData = stbi_load(path.c_str(), &width, &height, &channels, 0);

    image_data.resize(3*width*height, 0);
    TensorView image_tensor(image_data.data(), Shape{3, (uint)width, (uint)height});
    for (int i=0;i<height;i++)
    {
      for (int j=0;j<width;j++)
      {
        image_tensor.get(0, j, i) = imgData[channels*(i*width+j) + 0]/255.0f;
        image_tensor.get(1, j, i) = imgData[channels*(i*width+j) + 1]/255.0f;
        image_tensor.get(2, j, i) = imgData[channels*(i*width+j) + 2]/255.0f;
      }
    }

    stbi_image_free(imgData);

    return image_tensor;
  }

  void test_1_linear_regression()
  {
    std::vector<float> X, y;
    uint size = 1000;
    uint dim = 25;

    //y = 3*x + 1
    for (int i=0;i<size;i++)
    {
      float r = 1;
      for (int j=0;j<dim;j++)
      {
        float x0 = 2*((float)rand())/RAND_MAX - 1;
        X.push_back(x0);
        r += ((j%2) ? 0.5 : -0.5)*x0;
      }
      y.push_back(r);
    }

    TensorView Xv(X.data(), Shape{dim, size});
    TensorView yv(y.data(), Shape{1  , size});

    IndexType train_size = 900;
    IndexType val_size = 90;
    IndexType test_size = 10;

    TensorView X_train = slice(Xv, {0, train_size});
    TensorView y_train = slice(yv, {0, train_size});

    TensorView X_val = slice(Xv, {train_size, train_size+val_size});
    TensorView y_val = slice(yv, {train_size, train_size+val_size});

    TensorView X_test = slice(Xv, {train_size+val_size, size});
    TensorView y_test = slice(yv, {train_size+val_size, size});

    NeuralNetwork nn;
    nn.add_layer(std::make_shared<DenseLayer>(dim, 1));
    nn.train(X_train, y_train, X_val, y_val, 50, 3000, NeuralNetwork::Opt::Adam, NeuralNetwork::Loss::MSE, 0.01);
    nn.save_weights_to_file("weights.bin");

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(dim, 1));
    nn2.initialize_from_file("weights.bin");

    print(y_test);
    nn2.evaluate(X_test, y_test);
    print(y_test);
  }

  void test_2_simple_classification()
  {
    std::vector<float> X, y;
    uint size = 50000;
    uint dim = 5;

    //y = 3*x + 1
    for (int i=0;i<size;i++)
    {
      float r = 1;
      for (int j=0;j<dim;j++)
      {
        float x0 = 2*((float)rand())/RAND_MAX - 1;
        X.push_back(x0);
        r += x0;
      }
      y.push_back( sin(r) >= 0);
      y.push_back( sin(r) < 0);
    }

    TensorView Xv(X.data(), Shape{dim, size});
    TensorView yv(y.data(), Shape{2  , size});

    IndexType train_size = 0.8*size;
    IndexType val_size = 0.1*size;
    IndexType test_size = size-train_size-val_size;

    TensorView X_train = slice(Xv, {0, train_size});
    TensorView y_train = slice(yv, {0, train_size});

    TensorView X_val = slice(Xv, {train_size, train_size+val_size});
    TensorView y_val = slice(yv, {train_size, train_size+val_size});

    TensorView X_test = slice(Xv, {train_size+val_size, size});
    TensorView y_test = slice(yv, {train_size+val_size, size});

    NeuralNetwork nn;
    nn.add_layer(std::make_shared<DenseLayer>(dim, 64));
    nn.add_layer(std::make_shared<ReLULayer>());
    nn.add_layer(std::make_shared<DenseLayer>(64, 64));
    nn.add_layer(std::make_shared<ReLULayer>());
    nn.add_layer(std::make_shared<DenseLayer>(64, 2));
    nn.add_layer(std::make_shared<SoftMaxLayer>());
    nn.train(X_train, y_train, X_val, y_val, 256, 10000, NeuralNetwork::Opt::Adam, NeuralNetwork::Loss::CrossEntropy, 0.01);
  }

  void test_3_SIREN_image()
  {
    std::vector<float> image_data, pixel_data, image_data_grayscale;
    TensorView view = read_image_rgb("1a.png", image_data);
    IndexType pixel_count = view.size(1)*view.size(2);
    pixel_data.resize(2*pixel_count, 0);
    image_data_grayscale.resize(pixel_count, 0);
    for (IndexType i=0;i<view.size(2);i++)
    {
      for (IndexType j=0;j<view.size(1);j++)
      {
        pixel_data[2*(i*view.size(1) + j) + 0] = 2*(float)j/view.size(1)-1;
        pixel_data[2*(i*view.size(1) + j) + 1] = 2*(float)i/view.size(2)-1;
        image_data_grayscale[i*view.size(1) + j] = 0.2126 * view.get(0,j,i)+ 0.7152 * view.get(1,j,i) + 0.0722 * view.get(2,j,i);
        image_data_grayscale[i*view.size(1) + j] = 2*(image_data_grayscale[i*view.size(1) + j]) - 1;
      }
    }

    TensorView Xv = TensorView(pixel_data.data(), Shape{2, pixel_count}); //pixel coordinates
    TensorView yv = TensorView(image_data_grayscale.data(), Shape{1, pixel_count}); //list of pixels

    Siren siren(Siren::Type::Image, 3, 64);
    siren.train(Xv, yv, 1000, 10000);
    siren.evaluate(Xv, yv);
    yv = reshape(yv, Shape{view.size(1), view.size(2)});
    for (IndexType i=0;i<yv.total_size;i++)
      yv.get(i) = 0.5*(yv.get(i)+1);
    for (IndexType i=0;i<view.size(2);i++)
    {
      for (IndexType j=0;j<view.size(1);j++)
      {
        view.get(0,j,i) = yv.get(j,i);
        view.get(1,j,i) = yv.get(j,i);
        view.get(2,j,i) = yv.get(j,i);
      }
    }

    write_image_rgb("res.png", view);
  }

  void perform_tests()
  {
    test_1_linear_regression();
    test_2_simple_classification();
    test_3_SIREN_image();
  }
}