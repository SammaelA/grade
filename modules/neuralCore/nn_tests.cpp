#include "nn_tests.h"
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

#include "tensor_processor.h"
#include "tensor_compiler.h"
#include "neural_network.h"
#include "siren.h"

#include "stb_image.h"
#include "stb_image_write.h"

namespace nn
{
  void read_image_rgb(std::string path, std::vector<float> &image_data, int &width, int &height)
  {
    int channels;
    unsigned char *imgData = stbi_load(path.c_str(), &width, &height, &channels, 0);

    image_data.resize(3*width*height, 0);
    for (int i=0;i<height;i++)
    {
      for (int j=0;j<width;j++)
      {
        image_data[3*(i*width+j) + 0] = imgData[channels*(i*width+j) + 0]/255.0f;
        image_data[3*(i*width+j) + 1] = imgData[channels*(i*width+j) + 1]/255.0f;
        image_data[3*(i*width+j) + 2] = imgData[channels*(i*width+j) + 2]/255.0f;
      }
    }

    stbi_image_free(imgData);
  }

  void write_image_rgb(std::string path, const std::vector<float> &image_data, int width, int height)
  {
    assert(3*width*height == image_data.size());
    unsigned char *data = new unsigned char[3*width*height];
    for (int i=0;i<3*width*height;i++)
      data[i] = std::max(0, std::min(255, (int)(255*image_data[i])));
    
    stbi_write_png(path.c_str(), width, height, 3, data, 3*width);

    delete[] data;
  }

void test_1_tensor_processor()
  {
    printf("TEST 1. BASIC TENSOR PROGRAMS PROCESSING\n");

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

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::set_input("c", &c, 1);
    TensorProcessor::execute();
    TensorProcessor::get_output("A", A.data(), A.size());

    printf("  1.1. %-64s", "Correct result ");
    if (abs(A[0] - 0.25f) < 1e-6 && abs(A[1] - 0.25f) < 1e-6 && abs(A[2] - 0.5f) < 1e-6)
      printf("PASSED\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", A[0], A[1], A[2], 0.25f, 0.25f, 0.5f);
  }

  void test_2_tensor_tokens()
  {
    printf("TEST 2. BASIC TENSOR PROGRAMS COMPILING\n");

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

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("A", A.data(), A.size());

    printf("  2.1. %-64s", "Correct result ");
    if (abs(A[0] - 0.25f) < 1e-6 && abs(A[1] - 0.25f) < 1e-6 && abs(A[2] - 0.5f) < 1e-6)
      printf("PASSED\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", A[0], A[1], A[2], 0.25f, 0.25f, 0.5f);
  }

  void test_3_tensor_operations()
  {
    printf("TEST 3. VARIOUS TENSOR COMMANDS\n");

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
      TensorToken O = TensorToken::vector_outer_product(A, B.get({0,2}));
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

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::execute();
    for (int i=0;i<9;i++)
      TensorProcessor::get_output("R"+std::to_string(i+1), res.data()+i*9, 9);

    std::vector<std::pair<std::string, std::vector<float>>> reference = 
    {
      {{"R1"}, {1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 0.00, 0.00, 0.00,}},
      {{"R2"}, {14.00, -14.00, 32.00, -32.00, 0.00, 0.00, 0.00, 0.00, 0.00,}}, 
      {{"R3"}, {7.00, 8.00, 9.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,}}, 
      {{"R4"}, {50.00, -50.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,}}, 
      {{"R5"}, {1.00, -1.00, 2.00, -2.00, 3.00, -3.00, 0.00, 0.00, 0.00,}}, 
      {{"R6"}, {1.00, 2.00, 3.00, 2.00, 4.00, 6.00, 3.00, 6.00, 9.00,}},
      {{"R7"}, {-4.00, -5.00, -6.00, -8.00, -10.00, -12.00, -12.00, -15.00, -18.00,}}, 
      {{"R8"}, {171.00, 2.00, 3.00, 2.00, 4.00, 6.00, 3.00, 6.00, 9.00,}}, 
      {{"R9"}, {-4.00, -5.00, -6.00, -8.00, -10.00, -12.00, -12.00, -15.00, -18.00,}}
    };

    for (int k=0;k<9;k++)
    {
      printf("  3.%d. %-64s", k+1,(reference[k].first+" correct ").c_str());
      float diff = 0.0f;
      for (int i=0;i<9;i++)
        diff += abs(reference[k].second[i] - res[9*k + i]);
      if (diff < 1e-6)
        printf("PASSED\n");
      else
        printf("FAILED\n");
    }
  }

  void test_4_linear_regression()
  {
    printf("TEST 4. LINEAR REGRESSION EVALUATION\n");
    std::vector<float> X = {1,2,3, 1,3,2, 1,4,1};
    std::vector<float> w = {1,1,-1,0.01};
    std::vector<float> r = {0,0,0};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(3, 1));
    nn2.initialize_with_weights(w.data());
    nn2.evaluate(X, r);
    printf("  4.1. %-64s", "Correct result ");
    if (abs(r[0] - 0.01f) < 1e-6 && abs(r[1] - 2.01f) < 1e-6 && abs(r[2] - 4.01f) < 1e-6)
      printf("PASSED\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", r[0], r[1], r[2], 0.01f, 2.01f, 4.01f);
  }

  void test_5_aliases()
  {
    printf("TEST 5. ALIASES\n");
    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken X = TensorToken(2, 2, 3);
      TensorToken A = TensorToken(2, 2);
      TensorToken r = TensorToken(2, 2, 3);
      for (int i=0;i<3;i++)
      {
        TensorToken Xm = X.get(i);
        TensorToken rm = TensorToken(2, 2);
        for (int j=0;j<2;j++)
        {
          rm.set(j, TensorToken::mat_vec_mul(A, Xm.get(j)));
        }
        r.set(i, rm);
      }
      tc.input(X, "X");
      tc.input(A, "A");
      tc.output(r, "r");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> X = {1,2,3,4, 5,6,7,8, 9,10,11,12};
    std::vector<float> A = {1,1,2,2};
    std::vector<float> r = {0,0,0,0, 0,0,0,0, 0,0,0,0};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("X", X.data(), X.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("r", r.data(), r.size());

    std::vector<float> ref = {3, 6, 7, 14, 11, 22, 15, 30, 19, 38, 23, 46};
    float diff = 0.0f;
    for (int i=0;i<ref.size();i++)
      diff += std::abs(r[i] - ref[i]);
    
    printf("  5.1. %-64s", "Correct result ");
    if (diff < 1e-6)
      printf("PASSED\n");
    else
      printf("FAILED\n");
  }

  void test_6_linear_regression_train()
  {
    printf("TEST 6. LINEAR REGRESSION\n");
    int sz = 1000;
    int dim = 32;
    std::vector<float> X(dim*sz,0);
    std::vector<float> res(sz,0);
    for (int i=0;i<sz;i++)
    {
      float s = 0;
      for (int j=0;j<dim;j++)
      {
        X[dim*i + j] = 2*((double)rand())/RAND_MAX - 1;
        s += (j%2 ? 1 : -1)*X[dim*i + j];
      }
      res[i] = s + 0.01;
    }

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(dim, 64), NeuralNetwork::HE);
    nn2.add_layer(std::make_shared<DenseLayer>(64, 1), NeuralNetwork::HE);
    nn2.train(X, res, 512, 500, NeuralNetwork::Adam, NeuralNetwork::MSE, 0.01f);

    std::vector<float> y(sz,0);
    nn2.evaluate(X, y);
    float diff = 0.0f;
    for (int i=0;i<sz;i++)
      diff += abs(y[i]-res[i]);
    diff /= sz;
    printf("  6.1. %-64s", "Perfect optimization ");
    if (diff < 1e-6)
      printf("PASSED\n");
    else
      printf("FAILED\n");
  }

  void test_7_SIREN_image()
  {
    printf("TEST 7. SIREN 2D\n");
    std::vector<float> image_data, pixel_data, image_data_grayscale;
    int width=0, height=0;
    read_image_rgb("1a.png", image_data, width, height);
    assert(width*height > 0);
    int pixel_count = image_data.size()/3;
    pixel_data.resize(2*pixel_count, 0);
    image_data_grayscale.resize(pixel_count, 0);
    for (int i=0;i<height;i++)
    {
      for (int j=0;j<width;j++)
      {
        pixel_data[2*(i*width + j) + 0] = 2*(float)j/width-1;
        pixel_data[2*(i*width + j) + 1] = 2*(float)i/height-1;
        image_data_grayscale[i*width + j] = 0.2126 * image_data[3*(i*width + j)+0] + 
                                            0.7152 * image_data[3*(i*width + j)+1] + 
                                            0.0722 * image_data[3*(i*width + j)+2];
        image_data_grayscale[i*width + j] = 2*(image_data_grayscale[i*width + j]) - 1;
      }
    }

    NeuralNetwork nn2;
    nn2.set_batch_size_for_evaluate(2048);
    nn2.add_layer(std::make_shared<DenseLayer>( 2, 64), NeuralNetwork::SIREN);
    nn2.add_layer(std::make_shared<SinLayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 64), NeuralNetwork::SIREN);
    nn2.add_layer(std::make_shared<SinLayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 64), NeuralNetwork::SIREN);
    nn2.add_layer(std::make_shared<SinLayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64,  1), NeuralNetwork::SIREN);
    nn2.train(pixel_data, image_data_grayscale, 1000, 2500, NeuralNetwork::Adam, NeuralNetwork::MSE, 0.0001f);

    std::vector<float> image_data_res(pixel_count, 0);
    nn2.evaluate(pixel_data, image_data_res);

    for (float &val : image_data_res)
      val = 0.5*(val+1);
    for (int i=0;i<height;i++)
    {
      for (int j=0;j<width;j++)
      {
        image_data[3*(i*width + j)+0] = image_data_res[i*width+j];
        image_data[3*(i*width + j)+1] = image_data_res[i*width+j];
        image_data[3*(i*width + j)+2] = image_data_res[i*width+j];
      }
    }

    double diff = 0.0;
    for (int i=0;i<image_data_res.size();i++)
      diff += abs(image_data_res[i] - image_data_grayscale[i]);
    diff /= image_data_res.size();
    
    printf("  7.1. %-64s", "Decent image optimization ");
    if (diff < 0.25f)
      printf("PASSED\n");
    else
      printf("FAILED %f >= %f\n", diff, 0.25f);
  }

  void test_8_SIREN_SDF()
  {
    printf("TEST 8. SIREN 3D\n");
    auto circle_sdf = [](float x, float y, float z) -> float
    {
      return std::sqrt(x*x+y*y+z*z) - 0.75;
    };

    int samples = 25000;
    std::vector<float> points(3*samples);
    std::vector<float> distances(samples);

    for (int i=0;i<samples;i++)
    {
      float x = 2*((double)rand())/RAND_MAX - 1;
      float y = 2*((double)rand())/RAND_MAX - 1;
      float z = 2*((double)rand())/RAND_MAX - 1;

      points[3*i+0] = x;
      points[3*i+1] = y;
      points[3*i+2] = z;
      distances[i] = circle_sdf(x,y,z);
    }

    Siren siren(Siren::Type::SDF, 3, 64);
    siren.train(points, distances, 512, 1000);

    std::vector<float> predicted_distances(samples);
    siren.evaluate(points, predicted_distances);

    double diff = 0.0;
    for (int i=0;i<predicted_distances.size();i++)
      diff += abs(predicted_distances[i] - distances[i]);
    diff /= predicted_distances.size();
    
    printf("  8.1. %-64s", "Good simple SDF optimization ");
    if (diff < 0.01f)
      printf("PASSED\n");
    else
      printf("FAILED %f >= %f\n", diff, 0.01f);
  }

  void perform_tests()
  {
    printf("NEURAL CORE CPU TESTS\n");
    test_1_tensor_processor();
    test_2_tensor_tokens();
    test_3_tensor_operations();
    test_4_linear_regression();
    test_5_aliases();
    test_6_linear_regression_train();
    test_7_SIREN_image();
    test_8_SIREN_SDF();

    printf("NEURAL CORE GPU TESTS\n");
    TensorProcessor::init("GPU");
    test_1_tensor_processor();
    test_2_tensor_tokens();
    test_3_tensor_operations();
    test_4_linear_regression();
    test_5_aliases();
    test_6_linear_regression_train();
    test_7_SIREN_image();
    test_8_SIREN_SDF();
  }
}