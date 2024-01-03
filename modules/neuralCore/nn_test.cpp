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
#include "tensor_processor.h"
#include "tensor_compiler.h"
#include "neural_network_2.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../../third_party/third_party/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../third_party/third_party/stb_image_write.h"

using namespace nn;

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

  void test_4_tensor_processor()
  {
    std::vector<float> v = {1,2,3,  4,3,12, 5};
    std::vector<TensorProgram::Variable> vars(4);
    vars[0] = {1,0,3, {3,0,0,0}};//A
    vars[1] = {1,3,3, {3,0,0,0}};//B
    vars[2] = {0,6,1, {0,0,0,0}};//c
    vars[3] = {0,6,1, {0,0,0,0}};//s
    std::vector<TensorProgram::Command> commands(4);
    commands[0] = {TensorProgram::ADD, {0, 1, 0, 0, 0, 0}}; //A = A + B
    commands[1] = {TensorProgram::ADD, {0, 2, 0, 0, 0, 0}}; //A = A + c
    commands[2] = {TensorProgram::SUM, {0, 0, 3, 0, 0, 0}}; //s = sum(A)
    commands[3] = {TensorProgram::DIV, {0, 2, 0, 0, 0, 0}}; //A = A/c

    TensorProgram p;
    p.commands = commands;
    p.vars = vars;
    p.total_memory_req = 16;
    p.input_vars = std::map<std::string, unsigned>{{"A", 0},{"B", 1}, {"c", 2}};
    p.output_vars = std::map<std::string, unsigned>{{"A", 0}};

    std::vector<float> A = {1, 2, 3};
    std::vector<float> B = {4, 3,12};
    float c = 5;

    TensorProcessor tp;
    tp.set_program(p);
    tp.set_input({{"A", A.data()},{"B", B.data()}, {"c", &c}});
    tp.execute();
    tp.get_output({{"A", A.data()}});

    printf("res = ");
    for (int i=0;i<A.size();i++)
      printf("%f ", A[i]);
    printf("\n");
  }

  void test_5_tensor_tokens()
  {
    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(3);
      TensorToken B = TensorToken(3);
      TensorToken c = TensorToken(1);

      TensorToken res = A + B + 2.0f + 3.0f;
      res /= res.sum();

      tc.input(A, "A");
      tc.input(B, "B");
      tc.output(res, "A");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1, 2, 3};
    std::vector<float> B = {4, 3,12};

    TensorProcessor tp;
    tp.set_program(p);
    tp.set_input({{"A", A.data()},{"B", B.data()}});
    tp.execute();
    tp.get_output({{"A", A.data()}});

    printf("res = ");
    for (int i=0;i<A.size();i++)
      printf("%f ", A[i]);
    printf("\n");
  }

  void test_6_tensor_operations()
  {
    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(3, 2);
      TensorToken B = TensorToken(3, 3);
      TensorToken R1 = B.get({0, 2});
      TensorToken R2 = TensorToken::mat_mul_t(A, R1);
      TensorToken R3 = B.get(2);
      TensorToken R4 = TensorToken::mat_vec_mul(A, R3);
      TensorToken R5 = A.transpose();
      TensorToken O = TensorToken::vector_outer_product(A, B.get(0));
      TensorToken R6 = O.get(0);
      TensorToken R7 = O.get(1);
      TensorToken OF = O.flatten();
      OF.set(0, 171.0f);
      TensorToken R8 = OF.get({0, 9});
      TensorToken R9 = OF.get({9, 18});

      tc.input(A, "A");
      tc.input(B, "B");
      tc.output(R1, "R1");
      tc.output(R2, "R2");
      tc.output(R3, "R3");
      tc.output(R4, "R4");
      tc.output(R5, "R5");
      tc.output(R6, "R6");
      tc.output(R7, "R7");
      tc.output(R8, "R8");
      tc.output(R9, "R9");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1, 2, 3, -1, -2, -3};
    std::vector<float> B = {1,2,3,4,5,6,7,8,9};
    std::vector<float> res(90, 0.0f);

    TensorProcessor tp;
    tp.set_program(p);
    tp.set_input({{"A", A.data()},{"B", B.data()}});
    tp.execute();
    tp.get_output({{"R1", res.data()+0*9},{"R2", res.data()+1*9},{"R3", res.data()+2*9},
                   {"R4", res.data()+3*9},{"R5", res.data()+4*9},{"R6", res.data()+5*9},
                   {"R7", res.data()+6*9},{"R8", res.data()+7*9},{"R9", res.data()+8*9}});

    for (int k=0;k<10;k++)
    {
      printf("R%d = ", k+1);
      for (int i=0;i<9;i++)
        printf("%.2f ", res[9*k + i]);
      printf("\n");
    }
    /* Reference values
    R1 = 1.00 2.00 3.00 4.00 5.00 6.00 0.00 0.00 0.00 
    R2 = 14.00 -14.00 32.00 -32.00 0.00 0.00 0.00 0.00 0.00 
    R3 = 7.00 8.00 9.00 0.00 0.00 0.00 0.00 0.00 0.00 
    R4 = 50.00 -50.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 
    R5 = 1.00 -1.00 2.00 -2.00 3.00 -3.00 0.00 0.00 0.00 
    R6 = 1.00 2.00 3.00 2.00 4.00 6.00 3.00 6.00 9.00 
    R7 = -1.00 -2.00 -3.00 -2.00 -4.00 -6.00 -3.00 -6.00 -9.00 
    R8 = 171.00 2.00 3.00 2.00 4.00 6.00 3.00 6.00 9.00 
    R9 = -1.00 -2.00 -3.00 -2.00 -4.00 -6.00 -3.00 -6.00 -9.00 
    R10 = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 
    */
  }

  void test_7_linear_regression_2()
  {
    std::vector<float> X = {1,2,3,2,4,5};
    std::vector<float> w = {1,1,-1,0.01};
    std::vector<float> r = {0,0};
    NeuralNetwork2 nn2;
    nn2.add_layer(std::make_shared<DenseLayer2>(3, 1));
    nn2.initialize_with_weights(w.data());
    nn2.evaluate(X, r);
    printf("%f %f\n", r[0], r[1]);
    //nn.train(X_train, y_train, X_val, y_val, 50, 3000, NeuralNetwork::Opt::Adam, NeuralNetwork::Loss::MSE, 0.01);
  }

  int main(int argc, char **argv)
  {
    //test_1_linear_regression();
    //test_2_simple_classification();
    //test_3_SIREN_image();
    test_4_tensor_processor();
    test_5_tensor_tokens();
    test_6_tensor_operations();
    test_7_linear_regression_2();
    //std::vector<float> data;
    //TensorView view = read_image_rgb("empty_64.png", data);
    //printf("%d %d %d %d\n", view.Dim, view.size(0), view.size(1), view.size(2));
    //write_image_rgb("test_64.png", view);

    return 0;
  }