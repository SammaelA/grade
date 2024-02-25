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
#include "dataset.h"

#include "stb_image.h"
#include "stb_image_write.h"

namespace nn
{
  std::string base_path = "";

  void read_image_rgb(std::string path, std::vector<float> &image_data, int &width, int &height)
  {
    int channels;
    unsigned char *imgData = stbi_load((base_path+path).c_str(), &width, &height, &channels, 0);

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
    
    stbi_write_png((base_path+path).c_str(), width, height, 3, data, 3*width);

    delete[] data;
  }

void tp_test_1_tensor_processor()
  {
    printf("TEST 1. BASIC TENSOR PROGRAMS PROCESSING\n");

    std::vector<float> v = {1,2,3,  4,3,12, 5};
    std::vector<TensorProgram::Variable> vars(4);
    vars[0] = {1,0,3, {3,0,0,0}};//A
    vars[1] = {1,3,3, {3,0,0,0}};//B
    vars[2] = {0,6,1, {0,0,0,0}};//c
    vars[3] = {0,6,1, {0,0,0,0}};//s
    std::vector<TensorProgram::Command> commands(4);
    commands[0] = {TensorProgram::ADD, {0, 1, 0, 1, 3, 1, 0, 0}}; //A = A + B
    commands[1] = {TensorProgram::ADD, {0, 2, 0, 3, 1, 1, 0, 0}}; //A = A + c
    commands[2] = {TensorProgram::SUM, {0, 0, 3, 0, 0, 0, 0, 0}}; //s = sum(A)
    commands[3] = {TensorProgram::DIV, {0, 2, 0, 3, 1, 1, 0, 0}}; //A = A/s

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
      printf("passed\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", A[0], A[1], A[2], 0.25f, 0.25f, 0.5f);
  }

  void tp_test_2_tensor_tokens()
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
      printf("passed\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", A[0], A[1], A[2], 0.25f, 0.25f, 0.5f);
  }

  void tp_test_3_tensor_operations()
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
      TensorToken R4 = TensorToken::mat_mul_t(A, R3);
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
      {{"R2"}, {14.00, 32.00, -14.00, -32.00, 0.00, 0.00, 0.00, 0.00, 0.00,}}, 
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
        printf("passed\n");
      else
        printf("FAILED\n");
    }
  }

  void nn_test_1_linear_regression()
  {
    printf("TEST 1. LINEAR REGRESSION EVALUATION\n");
    std::vector<float> X = {1,2,3, 1,3,2, 1,4,1};
    std::vector<float> w = {1,1,-1,0.01};
    std::vector<float> r = {0,0,0};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(3, 1));
    nn2.initialize_with_weights(w.data());
    nn2.evaluate(X, r);
    printf("  1.1. %-64s", "Correct result ");
    if (abs(r[0] - 0.01f) < 1e-6 && abs(r[1] - 2.01f) < 1e-6 && abs(r[2] - 4.01f) < 1e-6)
      printf("passed\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", r[0], r[1], r[2], 0.01f, 2.01f, 4.01f);
  }

  void tp_test_4_aliases()
  {
    printf("TEST 4. ALIASES\n");
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
          rm.set(j, TensorToken::mat_mul_t(A, Xm.get(j)));
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
    
    printf("  4.1. %-64s", "Correct result ");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED\n");
  }

  void nn_test_2_linear_regression_train()
  {
    printf("TEST 2. LINEAR REGRESSION\n");
    int sz = 1024;
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
    nn2.add_layer(std::make_shared<DenseLayer>(dim, 64), Initializer::He);
    nn2.add_layer(std::make_shared<DenseLayer>(64, 1), Initializer::He);
    nn2.train(X.data(), res.data(), sz, 512, 512, true, OptimizerAdam(0.01f), Loss::MSE, Metric::MAE);

    std::vector<float> y(sz,0);
    nn2.evaluate(X, y);
    float diff = 0.0f;
    for (int i=0;i<sz;i++)
      diff += abs(y[i]-res[i]);
    diff /= sz;
    printf("  2.1. %-64s", "Perfect optimization ");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED %f > %f\n", diff, 1e-6);
  }

  void nn_test_3_SIREN_image()
  {
    printf("TEST 3. SIREN 2D\n");
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
    nn2.add_layer(std::make_shared<DenseLayer>( 2, 64), Initializer::Siren);
    nn2.add_layer(std::make_shared<SinLayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 64), Initializer::Siren);
    nn2.add_layer(std::make_shared<SinLayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 64), Initializer::Siren);
    nn2.add_layer(std::make_shared<SinLayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64,  1), Initializer::Siren);
    nn2.train(pixel_data, image_data_grayscale, 1000, 2500, OptimizerAdam(0.0001f), Loss::MSE);

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
    
    printf("  3.1. %-64s", "Decent image optimization ");
    if (diff < 0.25f)
      printf("passed\n");
    else
      printf("FAILED %f >= %f\n", diff, 0.25f);
  }

  void nn_test_4_SIREN_SDF()
  {
    printf("TEST 4. SIREN 3D\n");
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
    
    printf("  4.1. %-64s", "Good simple SDF optimization ");
    if (diff < 0.01f)
      printf("passed\n");
    else
      printf("FAILED %f >= %f\n", diff, 0.01f);
  }

  void tp_test_5_softmax()
  {
    printf("TEST 5. SOFTMAX FUNCTION\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken input = TensorToken(3, 2);

      TensorToken max_val = input.maximum(input.Dim-1) + 1e-15f;
      TensorToken output = TensorToken::g_2op(TensorProgram::SUB, input, max_val, 1);
      output = TensorToken::exp(output);
      TensorToken sum = output.sum(input.Dim-1);
      TensorToken res = TensorToken::g_2op(TensorProgram::DIV, output, sum, 1);

      tc.input(input, "A");
      tc.output(res, "A");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1, 2, 2, -1, -3, 2};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("A", A.data(), A.size());
    //printf("%f %f %f %f %f %f\n", A[0], A[1], A[2], A[3], A[4], A[5]);

    printf("  5.1. %-64s", "Correct result ");
    if (abs(A[0] - 0.155362f) < 1e-6 && abs(A[1] - 0.422319f) < 1e-6 && abs(A[5] - 0.946499f) < 1e-6)
      printf("passed\n");
    else
      printf("FAILED %f %f %f != %f %f %f\n", A[0], A[1], A[2], 0.25f, 0.25f, 0.5f);
    
  }

  void nn_test_5_simple_classifier()
  {
    printf("TEST 5. SIMPLE CLASSIFICATION\n");
    int sz = 1000;
    int dim = 2;
    std::vector<float> X(dim*sz,0);
    std::vector<float> res(2*sz,0);
    for (int i=0;i<sz;i++)
    {
      float s = 0;
      for (int j=0;j<dim;j++)
      {
        X[dim*i + j] = 2*((double)rand())/RAND_MAX - 1;
        s += X[dim*i + j];
      }
      res[2*i+0] = s > 0;
      res[2*i+1] = s <= 0;
    }

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(dim, 64), Initializer::He);
    nn2.add_layer(std::make_shared<DenseLayer>(64, 2), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(X, res, 512, 500, OptimizerAdam(0.01f), Loss::CrossEntropy);

    std::vector<float> y(2*sz,0);
    nn2.evaluate(X, y);
    float diff = 0.0f;
    for (int i=0;i<2*sz;i++)
    {
     // printf("(%f %f)\n", y[i], res[i]);
      diff += (y[i]>0.5) != (res[i]>0.5);
    }
    float error_rate = diff/sz;
    printf("  5.1. %-64s", "Error rate <3% ");
    if (error_rate < 0.03f)
      printf("passed\n");
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void nn_test_16_ReLU_classifier()
  {
    printf("TEST 16. CLASSIFICATION WITH RELU\n");
    int sz = 25000;
    int dim = 10;
    std::vector<float> X(dim*sz,0);
    std::vector<float> res(2*sz,0);
    for (int i=0;i<sz;i++)
    {
      float s = 0;
      for (int j=0;j<dim;j++)
      {
        X[dim*i + j] = 2*((double)rand())/RAND_MAX - 1;
        s += X[dim*i + j];
      }
      res[2*i+0] = s > 0;
      res[2*i+1] = s <= 0;
    }

    std::vector<float> X_test(dim*sz,0);
    std::vector<float> res_test(2*sz,0);
    for (int i=0;i<sz;i++)
    {
      float s = 0;
      for (int j=0;j<dim;j++)
      {
        X_test[dim*i + j] = 2*((double)rand())/RAND_MAX - 1;
        s += X_test[dim*i + j];
      }
      res_test[2*i+0] = s > 0;
      res_test[2*i+1] = s <= 0;
    }

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(dim, 64), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 64), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 2), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(X, res, 256, 5000, OptimizerAdam(0.01f), Loss::CrossEntropy);

    std::vector<float> y(2*sz,0);
    nn2.evaluate(X_test, y);
    float diff = 0.0f;
    for (int i=0;i<2*sz;i++)
      diff += (y[i]>0.5) != (res_test[i]>0.5);

    float error_rate = diff/sz;
    printf(" 16.1. %-64s", "Error rate <3% ");
    if (error_rate < 0.03f)
      printf("passed\n");
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void tp_test_6_logic_operations()
  {
    printf("TEST 6. LOGIC OPERATIONS\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(3, 2);
      TensorToken B = TensorToken(3);
      TensorToken R1 = TensorToken::g_2op(TensorProgram::GE, A, 0.0f);
      TensorToken R2 = TensorToken::g_2op(TensorProgram::WHERE, A, R1);
      TensorToken R3 = TensorToken::g_2op(TensorProgram::GREATER, A, B);
      TensorToken R4 = TensorToken::g_2op(TensorProgram::LESS, A, B);
      TensorToken R5 = TensorToken::g_2op(TensorProgram::GE, A, B);
      TensorToken R6 = TensorToken::g_2op(TensorProgram::LE, A, B);
      TensorToken R7 = TensorToken::g_2op(TensorProgram::EQUAL, A, B);
      TensorToken R8 = TensorToken::g_2op(TensorProgram::NE, A, B);
      TensorToken R9 = TensorToken::g_2op(TensorProgram::AND, R5, R6);
      TensorToken R10= TensorToken::g_2op(TensorProgram::OR, R3, R4);

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
      tc.output(R10,"R10");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1,2,0, -1,-3,2};
    std::vector<float> B = {0,2,3};
    std::vector<float> res(6*10, 0.0f);

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::execute();
    for (int i=0;i<10;i++)
      TensorProcessor::get_output("R"+std::to_string(i+1), res.data()+i*6, 6);
  
    std::vector<std::pair<std::string, std::vector<float>>> reference = 
    {
      {{"R1"}, {1,1,1,0,0,1}},
      {{"R2"}, {1,2,0,0,0,2}}, 
      {{"R3"}, {1,0,0,0,0,0}}, 
      {{"R4"}, {0,0,1,1,1,1}}, 
      {{"R5"}, {1,1,0,0,0,0}}, 
      {{"R6"}, {0,1,1,1,1,1}},
      {{"R7"}, {0,1,0,0,0,0}}, 
      {{"R8"}, {1,0,1,1,1,1}}, 
      {{"R9"}, {0,1,0,0,0,0}},
      {{"R10"},{1,0,1,1,1,1}}
    };

    for (int k=0;k<10;k++)
    {
      printf("  6.%2d. %-63s", k+1,(reference[k].first+" correct ").c_str());
      float diff = 0.0f;
      for (int i=0;i<6;i++)
        diff += abs(reference[k].second[i] - res[6*k + i]);
      if (diff < 1e-6)
        printf("passed\n");
      else
        printf("FAILED\n");
    }
  }

  void tp_test_7_padding()
  {
    printf("TEST 7. PADDING\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(3, 2);
      TensorToken PadX = A.add_padding(1, 2, 0);
      TensorToken PadY = A.add_padding(2, 1, 1);
      TensorToken PadXY = A.add_padding(2, 0, 0).add_padding(2, 0, 1);
      tc.input(A, "A");
      tc.output(PadX, "PadX");
      tc.output(PadY, "PadY");
      tc.output(PadXY, "PadXY");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1,2,3,4,5,6};
    std::vector<float> PadX(2*6, 0.0f), PadX_ref = {0,1,2,3,0,0,
                                                    0,4,5,6,0,0};
    std::vector<float> PadY(5*3, 0.0f), PadY_ref = {0,0,0,
                                                    0,0,0,
                                                    1,2,3,
                                                    4,5,6,
                                                    0,0,0};
    std::vector<float> PadXY(4*5, 0.0f), PadXY_ref = {0,0,0,0,0,
                                                      0,0,0,0,0,
                                                      0,0,1,2,3,
                                                      0,0,4,5,6,};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("PadX", PadX.data(), PadX.size());
    TensorProcessor::get_output("PadY", PadY.data(), PadY.size());
    TensorProcessor::get_output("PadXY", PadXY.data(), PadXY.size());

    {
    float diff = 0.0f;
    for (int i=0;i<PadX.size();i++)
      diff += abs(PadX[i] - PadX_ref[i]);
    
    printf("  7.1. %-64s","X padding");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<PadY.size();i++)
      diff += abs(PadY[i] - PadY_ref[i]);
    
    printf("  7.2. %-64s","Y padding");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<PadXY.size();i++)
      diff += abs(PadXY[i] - PadXY_ref[i]);
    
    printf("  7.3. %-64s","X and Y padding");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
  }

  void tp_test_8_conv2D()
  {
    printf("TEST 8. 2D CONVOLUTION\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(5,5);
      TensorToken B = TensorToken(5,5,2,3);
      TensorToken B_kernel = TensorToken(3,3,2,2);
      TensorToken dx_kernel = TensorToken(3,1);
      TensorToken dy_kernel = TensorToken(1,3);
      TensorToken R1 = TensorToken::conv2D(A, dx_kernel);
      TensorToken R2 = TensorToken::conv2D(A, dy_kernel);
      TensorToken sum_kernel = TensorToken(3,3);
      TensorToken R3 = TensorToken::conv2D(A, sum_kernel, 2);
      TensorToken R4 = TensorToken::conv2D(B, B_kernel);
      tc.input(A, "A");
      tc.input(B, "B");
      tc.input(dx_kernel, "dx_kernel");
      tc.input(dy_kernel, "dy_kernel");
      tc.input(sum_kernel, "sum_kernel");
      tc.input(B_kernel, "B_kernel");
      tc.output(R1, "R1");
      tc.output(R2, "R2");
      tc.output(R3, "R3");
      tc.output(R4, "R4");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {0 +0,0 +1,0 +4,0 +9,0 +16,
                            1 +0,1 +1,1 +4,1 +9,1 +16,
                            4 +0,4 +1,4 +4,4 +9,4 +16,
                            9 +0,9 +1,9 +4,9 +9,9 +16,
                            16+0,16+1,16+4,16+9,16+16,};
    std::vector<float> B = {0, 1, 4, 9, 16,
                            0, 1, 4, 9, 16,
                            0, 1, 4, 9, 16,
                            0, 1, 4, 9, 16,
                            0, 1, 4, 9, 16,
                            
                            0, 0, 0, 0, 0,
                            1, 1, 1, 1, 1,
                            4, 4, 4, 4, 4,
                            9, 9, 9, 9, 9,
                            16,16,16,16,16,
                            
                            
                            -0, -1, -4, -9, -16,
                            -0, -1, -4, -9, -16,
                            -0, -1, -4, -9, -16,
                            -0, -1, -4, -9, -16,
                            -0, -1, -4, -9, -16,
                            
                            -0, -0, -0, -0, -0,
                            -1, -1, -1, -1, -1,
                            -4, -4, -4, -4, -4,
                            -9, -9, -9, -9, -9,
                            -16,-16,-16,-16,-16,


                            1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1,
                            
                            1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1,
                            1,1,1,1,1};
    
    std::vector<float> B_kernel = {0,0,0,-1,0,1,0,0,0,  0,0,0,-1,0,1,0,0,0,
                                   0,-1,0,0,0,0,0,1,0,  0,-1,0,0,0,0,0,1,0};

    std::vector<float> k1 = {-1,0,1};
    std::vector<float> k2 = {1,1,1,
                             1,1,1,
                             1,1,1};
    std::vector<float> R1(15, 0.0f), R1_ref = {4, 8, 12, 4, 8, 12, 4, 8, 12, 4, 8, 12, 4, 8, 12,};
    std::vector<float> R2(15, 0.0f), R2_ref = {4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 12, 12, 12, 12, 12,};
    std::vector<float> R3(4, 0.0f), R3_ref = {30, 102, 102, 174,};
    std::vector<float> R4(3*3*2*3, 0.0f), R4_ref = {4,8,12, 4,8,12, 4,8,12, 
                                                    4,4,4, 8,8,8, 12,12,12, 
                                                    -4,-8,-12, -4,-8,-12, -4,-8,-12, 
                                                    -4,-4,-4, -8,-8,-8, -12,-12,-12, 
                                                    0,0,0, 0,0,0, 0,0,0, 
                                                    0,0,0, 0,0,0, 0,0,0, };

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::set_input("dx_kernel", k1.data(), k1.size());
    TensorProcessor::set_input("dy_kernel", k1.data(), k1.size());
    TensorProcessor::set_input("sum_kernel", k2.data(), k2.size());
    TensorProcessor::set_input("B_kernel", B_kernel.data(), B_kernel.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("R1", R1.data(), R1.size());
    TensorProcessor::get_output("R2", R2.data(), R2.size());
    TensorProcessor::get_output("R3", R3.data(), R3.size());
    TensorProcessor::get_output("R4", R4.data(), R4.size());

    /*
    for (auto &v : R1)
      printf("%f, ", v);
    printf("\n");
    for (auto &v : R2)
      printf("%f, ", v);
    printf("\n");
    for (auto &v : R3)
      printf("%f, ", v);
    printf("\n");    
    for (auto &v : R4)
      printf("%f, ", v);
    printf("\n");*/


    {
    float diff = 0.0f;
    for (int i=0;i<R1.size();i++)
      diff += abs(R1[i] - R1_ref[i]);
    
    printf("  8.1. %-64s","kernel X");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<R2.size();i++)
      diff += abs(R2[i] - R2_ref[i]);
    
    printf("  8.2. %-64s","kernel Y");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<R3.size();i++)
      diff += abs(R3[i] - R3_ref[i]);
    
    printf("  8.3. %-64s","kernel with stride > 1");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<R4.size();i++)
      diff += abs(R4[i] - R4_ref[i]);
    
    printf("  8.4. %-64s","multi-layered kernel");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
  }

  void tp_test_9_conv2D_blur()
  {
    printf("TEST 9. BLUR WITH 2D CONVOLUTION\n");

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

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken image = TensorToken(width, height);
      TensorToken blur_kernel = TensorToken(11,11);
      TensorToken blurred_image = TensorToken::conv2D(image.add_padding(5,5,0).add_padding(5,5,1), blur_kernel);
      tc.input(image, "image");
      tc.input(blur_kernel, "blur_kernel");
      tc.output(blurred_image, "blurred_image");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> blur_kernel(11*11, 1.0/(11*11));
    std::vector<float> image_data_res(pixel_count, 0);

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("image", image_data_grayscale.data(), image_data_grayscale.size());
    TensorProcessor::set_input("blur_kernel", blur_kernel.data(), blur_kernel.size());
    auto t_prev = std::chrono::steady_clock::now();
    TensorProcessor::execute();
    auto t_now = std::chrono::steady_clock::now();
    float ms = 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_prev).count();
    TensorProcessor::get_output("blurred_image", image_data_res.data(), image_data_res.size());


    for (int i=0;i<height;i++)
    {
      for (int j=0;j<width;j++)
      {
        image_data[3*(i*width + j)+0] = image_data_res[i*width+j];
        image_data[3*(i*width + j)+1] = image_data_res[i*width+j];
        image_data[3*(i*width + j)+2] = image_data_res[i*width+j];
      }
    }
    write_image_rgb("1a_blurred.png", image_data, width, height);
    
    for (int i=0;i<height;i++)
    {
      for (int j=0;j<width;j++)
      {
        image_data[3*(i*width + j)+0] = image_data_grayscale[i*width+j];
        image_data[3*(i*width + j)+1] = image_data_grayscale[i*width+j];
        image_data[3*(i*width + j)+2] = image_data_grayscale[i*width+j];
      }
    }
    write_image_rgb("1a_gray.png", image_data, width, height);
    
    printf("  9.1. blur took %4.1f ms                                               ", ms);
    printf("passed\n");
  }

  void nn_test_6_conv2D_forward()
  {
    printf("TEST 6. CONV_2D FORWARD PASS\n");
    std::vector<float> X = {2 ,5 ,10, 
                            5 ,8 ,13,
                            10,13,18};
    std::vector<float> w = {0,-1,0,0,0,0,0,1,0, 0,0,0,-1,0,1,0,0,0, 0,0};
    std::vector<float> r(18, 0.0f);
    std::vector<float> r_ref = {5, 8, 13, 8, 8, 8, -5, -8, -13, 5, 8, -5, 8, 8, -8, 13, 8, -13,};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(3,3,1, 2, 3, 1, Conv2DLayer::SAME, true));
    nn2.add_layer(std::make_shared<FlattenLayer>(3,3,2));
    nn2.initialize_with_weights(w.data());
    nn2.evaluate(X, r);

    float diff = 0.0f;
    for (int i=0;i<r_ref.size();i++)
      diff += std::abs(r[i] - r_ref[i]);
    
    printf("  6.1. %-64s", "Correct result ");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED\n");

    //for (auto &v : r)
    //  printf("%f, ", v);
    //printf("\n");
  }

  void nn_test_7_conv2D_backward()
  {
    printf("TEST 7. CONV_2D BACKWARD PASS\n");
    std::vector<float> X = {2 ,5 ,10, 
                            5 ,8 ,13,
                            10,13,18};
    std::vector<float> w = {0,-1,0,0,0,0,0,1,0, 0,0,0,-1,0,1,0,0,0};
    std::vector<float> r(18, 0.0f);
    std::vector<float> r_ref = {5, 8, 13, 8, 8, 8, -5, -8, -13, 5, 8, -5, 8, 8, -8, 13, 8, -13,};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(3,3,1, 4, 3, 1, Conv2DLayer::SAME), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<Conv2DLayer>(3,3,4, 2, 3, 1, Conv2DLayer::SAME), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<FlattenLayer>(3,3,2));
    nn2.add_layer(std::make_shared<DenseLayer>(18, 18), Initializer::He);
    nn2.train(X, r_ref, 1, 500, OptimizerAdam(0.001f), Loss::MSE);
    nn2.evaluate(X, r);
    //for (auto &v : r)
    //  printf("%f, ", v);
    //printf("\n");
    
    float diff = 0.0f;
    for (int i=0;i<r_ref.size();i++)
      diff += std::abs(r[i] - r_ref[i]);
    
    printf("  7.1. %-64s", "Correct result ");
    if (diff < 1.0f)
      printf("passed\n");
    else
      printf("FAILED diff %f > %f\n", diff, 1.0f);
  }

  void nn_test_8_conv2D_no_padding()
  {
    printf("TEST 8. CONV_2D NO PADDING\n");
    std::vector<float> X = {1,2,3,4,5,6,7,8,9};
    std::vector<float> w = {0,-1,0,0,0,0,0,1,0, 0,0,0,-1,0,1,0,0,0};
    std::vector<float> r(1, 0.0f);
    std::vector<float> r_ref = {45, 15};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(3,3,1, 2), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<FlattenLayer>(1,1,2));
    nn2.add_layer(std::make_shared<DenseLayer>(2, 2), Initializer::He);
    nn2.train(X, r_ref, 1, 500, OptimizerAdam(0.01f), Loss::MSE);
    nn2.evaluate(X, r);
    //for (auto &v : r)
    //  printf("%f, ", v);
    //printf("\n");
    
    float diff = 0.0f;
    for (int i=0;i<r_ref.size();i++)
      diff += std::abs(r[i] - r_ref[i]);
    
    printf("  8.1. %-64s", "Correct result ");
    if (diff < 1e-4)
      printf("passed\n");
    else
      printf("FAILED\n");
  }

  void tp_test_11_max_pooling()
  {
    printf("TEST 11. MAX POOLING\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(4, 4);
      TensorToken B = TensorToken(2, 2);
      TensorToken::issue_command(TensorProgram::MPOOL, A, A, B, 2, 2);
      TensorToken dOut = TensorToken(2, 2);
      TensorToken dIn = TensorToken(4, 4);
      TensorToken::issue_command(TensorProgram::MPOOL_D, A, dOut, dIn, 2, 2);
      tc.input(A, "A");
      tc.input(dOut, "dOut");
      tc.output(B, "B");
      tc.output(dIn, "dIn");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1,0, -2,-1,
                            0,0, -3,-4,

                            3,4,  6, 7,
                            7,0,  8, 9};
    std::vector<float> dOut = {1,2,3,4};
    std::vector<float> B(2*2, 0.0f), B_ref = {1,-1,7,9};
    std::vector<float> dIn(4*4, 0.0f), dIn_ref = {1,0, 0,2,
                                                  0,0, 0,0,

                                                  0,0, 0,0,
                                                  3,0, 0,4};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("dOut", dOut.data(), dOut.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("B", B.data(), B.size());
    TensorProcessor::get_output("dIn", dIn.data(), dIn.size());
    //for (auto &v : B)
    //  printf("%f, ", v);
    //printf("\n");
    //for (auto &v : dIn)
    //  printf("%f, ", v);
    //printf("\n");
    {
    float diff = 0.0f;
    for (int i=0;i<B.size();i++)
      diff += abs(B[i] - B_ref[i]);
    
    printf(" 11.1. %-64s","forward pass");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<dIn.size();i++)
      diff += abs(dIn[i] - dIn_ref[i]);
    
    printf(" 11.2. %-64s","backward pass");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
  }

  void nn_test_9_synthetic_images_classifier()
  {
    printf("TEST 9. SYNTHETIC IMAGES CLASSIFIER\n");
    unsigned image_size = 16;
    unsigned image_count = 5000;

    std::vector<float> X(image_count*image_size*image_size, 0.0f);
    std::vector<float> y(2*image_count, 0.0f);

    for (int i=0;i<image_count;i++)
    {
      int type = rand()%2;
      int line = rand()%image_size;
      y[2*i + type] = 1;

      if (type == 0) //horizontal line
      {
        for (int j=0;j<image_size;j++)
          X[i*image_size*image_size + line*image_size + j] = 1.0f;
      }
      else //vertical line
      {
        for (int j=0;j<image_size;j++)
          X[i*image_size*image_size + j*image_size + line] = 1.0f;
      }
    }

    int train_sz = 0.9*image_count;
    std::vector<float> X_train = std::vector<float>(X.begin(), X.begin() + train_sz*image_size*image_size);
    std::vector<float> X_test  = std::vector<float>(X.begin() + train_sz*image_size*image_size, X.end());
    std::vector<float> y_train = std::vector<float>(y.begin(), y.begin() + train_sz*2);
    std::vector<float> y_test  = std::vector<float>(y.begin() + train_sz*2, y.end());

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(16,16,1, 8), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<MaxPoolingLayer>(14, 14, 8));
    nn2.add_layer(std::make_shared<Conv2DLayer>(7,7,8, 16), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<FlattenLayer>(5,5,16));
    nn2.add_layer(std::make_shared<DenseLayer>(5*5*16, 64), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 2), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(X_train, y_train, 64, 1000, OptimizerAdam(0.0002f), Loss::CrossEntropy);

    std::vector<float> y_res(y_test.size(),0);
    nn2.evaluate(X_test, y_res);
    float diff = 0.0f;
    for (int i=0;i<y_test.size();i++)
      diff += (y_res[i]>0.5) != (y_test[i]>0.5);

    float error_rate = diff/(y_test.size());
    printf("  9.1. %-64s", "Error rate <10% ");
    if (error_rate < 0.1f)
      printf("passed\n");
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void nn_test_10_MNIST()
  {
    printf("TEST 10. FULLY CONNECTED NN TRAINING ON MNIST DATASET\n");
    Dataset dataset;
    read_MNIST_dataset(base_path + std::string("../../resources/MNIST-dataset"), &dataset);
    train_test_split(&dataset, 0.1);

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<FlattenLayer>(28,28,1));
    nn2.add_layer(std::make_shared<DenseLayer>(28*28*1, 200), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(200, 200), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(200, 10), Initializer::He);
    //nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    /*
    nn2.add_layer(std::make_shared<Conv2DLayer>(32,32,3, 32, 3, Conv2DLayer::SAME), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<MaxPoolingLayer>(32, 32, 32));
    nn2.add_layer(std::make_shared<Conv2DLayer>(16,16,32, 64), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<MaxPoolingLayer>(14, 14, 64));
    nn2.add_layer(std::make_shared<FlattenLayer>(7,7,64));
    nn2.add_layer(std::make_shared<DenseLayer>(7*7*64, 200), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(200, 10), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());*/
    nn2.train(dataset.train_data, dataset.train_labels, 128, 5000, OptimizerAdam(0.0001f), Loss::CrossEntropy);

    std::vector<float> y_res(dataset.test_labels.size(),0);
    nn2.evaluate(dataset.test_data, y_res);
    float acc = 0.0f;
    float cnt = 0.0f;
    for (int i=0;i<dataset.test_labels.size();i+=10)
    {
      int max_pos = 0;
      int ref_max_pos = 0;
      for (int j=0;j<10;j++)
      {
        if (y_res[i+j] > y_res[i+max_pos])
          max_pos = j;
        if (dataset.test_labels[i+j] > dataset.test_labels[i+ref_max_pos])
          ref_max_pos = j;
      }

      acc += (max_pos == ref_max_pos);
      cnt++;
    }
    float error_rate = 1 - acc/cnt;
    printf(" 10.1. %-64s", "Error rate <5% ");
    if (error_rate < 0.1f)
      printf("passed\n");
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void tp_test_12_arithmetics_benchmark()
  {
    printf("TEST 12. ARITHMETICS BENCHMARK\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(64000, 100);
      TensorToken B = TensorToken(64000);
      for (int i=0;i<100;i++)
      A = A * B;
      tc.inout(A, "A");
      tc.input(B, "B");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A(100*64000, 1.5f);
    std::vector<float> B(64000, 1.0f);

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    auto t_prev = std::chrono::steady_clock::now();
    TensorProcessor::execute();
    auto t_now = std::chrono::steady_clock::now();
    float ms = 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_prev).count();
    printf(" 12.1. operation took %6.1f ms                                        ", ms);
    printf("passed\n");
  }

  void nn_test_11_CIFAR10()
  {
    printf("TEST 11. CONVOLUTIONAL NN TRAINING ON CIFAR10 DATASET\n");
    printf("TEMPORARY DISABLED. test takes too long\n");
    return;
    Dataset dataset;
    read_CIFAR10_dataset("../../resources/cifar-10-dataset", &dataset);
    train_test_split(&dataset, 0.1);

    NeuralNetwork nn2;
    //nn2.add_layer(std::make_shared<Conv2DLayer>(32,32,3, 32, 5), Initializer::GlorotNormal);
    //nn2.add_layer(std::make_shared<ReLULayer>());
    //nn2.add_layer(std::make_shared<MaxPoolingLayer>(28, 28, 32));
    //nn2.add_layer(std::make_shared<Conv2DLayer>(14,14,32, 32), Initializer::GlorotNormal);
    //nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<Conv2DLayer>(32,32,3, 8, 5), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<MaxPoolingLayer>(28, 28, 8));
    nn2.add_layer(std::make_shared<FlattenLayer>(14,14,8));
    nn2.add_layer(std::make_shared<DenseLayer>(14*14*8, 512), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(512, 512), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(512, 10), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(dataset.train_data, dataset.train_labels, 128, 25000, OptimizerAdam(0.001f), Loss::CrossEntropy);

    std::vector<float> y_res(dataset.test_labels.size(),0);
    nn2.evaluate(dataset.test_data, y_res);
    float acc = 0.0f;
    float cnt = 0.0f;
    for (int i=0;i<dataset.test_labels.size();i+=10)
    {
      int max_pos = 0;
      int ref_max_pos = 0;
      for (int j=0;j<10;j++)
      {
        if (y_res[i+j] > y_res[i+max_pos])
          max_pos = j;
        if (dataset.test_labels[i+j] > dataset.test_labels[i+ref_max_pos])
          ref_max_pos = j;
      }

      acc += (max_pos == ref_max_pos);
      cnt++;
    }
    float error_rate = 1 - acc/cnt;
    printf(" 11.1. %-64s", "Error rate <60% ");
    if (error_rate < 0.6f)
      printf("passed\n");
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void tp_test_13_dilation()
  {
    printf("TEST 13. DILATION\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(3, 2);
      TensorToken B = TensorToken(3, 3);
      TensorToken A_res = TensorToken(7, 2);
      TensorToken B_res = TensorToken(5, 7);
      TensorToken::issue_command(TensorProgram::DILATE, A, A, A_res, 2);
      TensorToken::issue_command(TensorProgram::DILATE, B, B, B_res, 1, 2);

      tc.input(A, "A");
      tc.input(B, "B");
      tc.output(A_res, "A_res");
      tc.output(B_res, "B_res");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1,2,3,4,5,6}, A_res(14), A_ref = {1,0,0,2,0,0,3, 4,0,0,5,0,0,6};
    std::vector<float> B = {1,2,3,4,5,6,7,8,9}, B_res(35), B_ref = {1,0,2,0,3,
                                                                    0,0,0,0,0,
                                                                    0,0,0,0,0,
                                                                    4,0,5,0,6,
                                                                    0,0,0,0,0,
                                                                    0,0,0,0,0,
                                                                    7,0,8,0,9};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("A_res", A_res.data(), A_res.size());
    TensorProcessor::get_output("B_res", B_res.data(), B_res.size());
    //for (auto &v : A_res)
    //  printf("%f, ", v);
    //printf("\n");
    //for (auto &v : B_res)
    //  printf("%f, ", v); q
    //printf("\n");
    {
    float diff = 0.0f;
    for (int i=0;i<A_res.size();i++)
      diff += abs(A_res[i] - A_ref[i]);
    
    printf(" 13.1. %-64s","1D dilation");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<B_res.size();i++)
      diff += abs(B_res[i] - B_ref[i]);
    
    printf(" 13.2. %-64s","2D dilation");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }    
  }

  void nn_test_12_conv2D_stride()
  {
    printf("TEST 12. CONV_2D STRIDE\n");
    std::vector<float> X = { 1, 1,0,-1,-1,
                             1, 1,0,-1,-1,
                             0, 0,0, 0, 0,
                            -1,-1,0, 1, 1,
                            -1,-1,0, 1, 1,};
    std::vector<float> r(8, 0.0f);
    std::vector<float> r_ref = {1,-1,-1,1, -1,1,1,-1};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(5,5,1, 2, 3, 2), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<FlattenLayer>(2,2,2));
    nn2.add_layer(std::make_shared<DenseLayer>(8, 8), Initializer::He);
    nn2.train(X, r_ref, 1, 1000, OptimizerAdam(0.001f), Loss::MSE);
    nn2.evaluate(X, r);
    //for (auto &v : r)
    //  printf("%f, ", v);
    //printf("\n");
    
    float diff = 0.0f;
    for (int i=0;i<r_ref.size();i++)
      diff += std::abs(r[i] - r_ref[i]);
    
    printf(" 12.1. %-64s", "Correct result ");
    if (diff < 1.0f)
      printf("passed\n");
    else
      printf("FAILED diff %f > %f\n", diff, 1.0f);
  }

  void nn_test_13_binary_classification_metrics()
  {
    printf("TEST 13. BINARY CLASSIFICATION METRICS\n");
    unsigned image_size = 16;
    unsigned image_count = 5000;

    std::vector<float> X(image_count*image_size*image_size, 0.0f);
    std::vector<float> y(2*image_count, 0.0f);

    for (int i=0;i<image_count;i++)
    {
      int type = rand()%2;
      int line = rand()%image_size;
      y[2*i + type] = 1;

      if (type == 0) //horizontal line
      {
        for (int j=0;j<image_size;j++)
          X[i*image_size*image_size + line*image_size + j] = 1.0f;
      }
      else //vertical line
      {
        for (int j=0;j<image_size;j++)
          X[i*image_size*image_size + j*image_size + line] = 1.0f;
      }
    }

    int train_sz = 0.9*image_count;
    std::vector<float> X_train = std::vector<float>(X.begin(), X.begin() + train_sz*image_size*image_size);
    std::vector<float> X_test  = std::vector<float>(X.begin() + train_sz*image_size*image_size, X.end());
    std::vector<float> y_train = std::vector<float>(y.begin(), y.begin() + train_sz*2);
    std::vector<float> y_test  = std::vector<float>(y.begin() + train_sz*2, y.end());

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(16,16,1, 8), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<MaxPoolingLayer>(14, 14, 8));
    nn2.add_layer(std::make_shared<Conv2DLayer>(7,7,8, 16), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<FlattenLayer>(5,5,16));
    nn2.add_layer(std::make_shared<DenseLayer>(5*5*16, 64), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 2), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(X_train, y_train, 64, 1000, OptimizerAdam(0.0002f), Loss::CrossEntropy);

    std::vector<float> y_res(y_test.size(),0);
    nn2.evaluate(X_test, y_res);

    float acc = nn2.calculate_metric(y_res.data(), y_test.data(), y_test.size()/2, Metric::Accuracy);
    float pr = nn2.calculate_metric(y_res.data(), y_test.data(), y_test.size()/2, Metric::Precision);
    float recall = nn2.calculate_metric(y_res.data(), y_test.data(), y_test.size()/2, Metric::Recall);
    float auc_roc = nn2.calculate_metric(y_res.data(), y_test.data(), y_test.size()/2, Metric::AUC_ROC);
    float auc_pr = nn2.calculate_metric(y_res.data(), y_test.data(), y_test.size()/2, Metric::AUC_PR);

    printf(" 13.1. %-64s", "High accuracy");
    if (acc > 0.95)
      printf("passed\n");
    else
      printf("FAILED, accuracy %f\n", acc);

    printf(" 13.2. %-64s", "High precision");
    if (pr > 0.95)
      printf("passed\n");
    else
      printf("FAILED, precision %f\n", pr);
  
    printf(" 13.3. %-64s", "High recall");
    if (recall > 0.95)
      printf("passed\n");
    else
      printf("FAILED, recall %f\n", recall);
      
    printf(" 13.4. %-64s", "High AUC ROC");
    if (auc_roc > 0.95)
      printf("passed\n");
    else
      printf("FAILED, AUC ROC %f\n", auc_roc);
      
    printf(" 13.5. %-64s", "High AUC PR");
    if (auc_pr > 0.95)
      printf("passed\n");
    else
      printf("FAILED, AUC PR %f\n", auc_pr);
  }

  void nn_test_14_batch_normalization()
  {
    printf("TEST 14. BATCH NORMALIZATION\n");
    unsigned image_size = 16;
    unsigned image_count = 5000;

    std::vector<float> X(image_count*image_size*image_size, 0.0f);
    std::vector<float> y(2*image_count, 0.0f);

    for (int i=0;i<image_count;i++)
    {
      int type = rand()%2;
      int line = rand()%image_size;
      y[2*i + type] = 1;

      if (type == 0) //horizontal line
      {
        for (int j=0;j<image_size;j++)
          X[i*image_size*image_size + line*image_size + j] = 1.0f;
      }
      else //vertical line
      {
        for (int j=0;j<image_size;j++)
          X[i*image_size*image_size + j*image_size + line] = 1.0f;
      }
    }

    int train_sz = 0.9*image_count;
    std::vector<float> X_train = std::vector<float>(X.begin(), X.begin() + train_sz*image_size*image_size);
    std::vector<float> X_test  = std::vector<float>(X.begin() + train_sz*image_size*image_size, X.end());
    std::vector<float> y_train = std::vector<float>(y.begin(), y.begin() + train_sz*2);
    std::vector<float> y_test  = std::vector<float>(y.begin() + train_sz*2, y.end());

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv2DLayer>(16,16,1, 8), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<BatchNormLayer>(), Initializer::BatchNorm);
    nn2.add_layer(std::make_shared<MaxPoolingLayer>(14, 14, 8));
    nn2.add_layer(std::make_shared<Conv2DLayer>(7,7,8, 16), Initializer::GlorotNormal);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<BatchNormLayer>(), Initializer::BatchNorm);
    nn2.add_layer(std::make_shared<FlattenLayer>(5,5,16));
    nn2.add_layer(std::make_shared<DenseLayer>(5*5*16, 64), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(64, 2), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(X_train.data(), y_train.data(), train_sz, 64, 2, true, OptimizerAdam(0.001f), Loss::CrossEntropy, 
              Metric::Accuracy);


    std::vector<float> y_res(y_test.size(),0);
    nn2.evaluate(X_test, y_res);
    float acc = nn2.calculate_metric(y_res.data(), y_test.data(), y_test.size()/2, Metric::Accuracy);
    printf(" 14.1. %-64s", "Accuracy > 95% ");
    if (acc > 0.95f)
      printf("passed\n");
    else
      printf("FAILED, accuracy %f\n", acc);
  }

  void tp_test_10_random()
  {
    srand(time(NULL));
    printf("TEST 10. RANDOM\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(3, 3);
      A.random(1u); //should be the same each iteration
      TensorToken B = TensorToken(100000);
      B.random(); //should change each iteration
      tc.output(A, "A_res");
      tc.output(B, "B_res");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A_res(9), A_ref = {0.22637305, 0.70636714, 0.25718519, 0.41423926, 0.71294135, 0.68870342, 0.87693751, 0.81305557, 0.03246965};
    std::vector<float> B_res(10000);

    TensorProcessor::set_program(p);
    TensorProcessor::execute();
    TensorProcessor::get_output("A_res", A_res.data(), A_res.size());
    TensorProcessor::get_output("B_res", B_res.data(), B_res.size());
    {
    float diff = 0.0f;
    for (int i=0;i<A_res.size();i++)
      diff += abs(A_res[i] - A_ref[i]);
    
    printf(" 10.1. %-64s","Same seed gives same result");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED\n");
    }
    {
    float mean = 0.0f;
    for (int i=0;i<B_res.size();i++)
      mean += B_res[i];
    mean /= B_res.size();
    float variance = 0.0;
    for (int i=0;i<B_res.size();i++)
      variance += (B_res[i] - mean)*(B_res[i] - mean);
    variance = variance/B_res.size();
    printf(" 10.2. %-64s","Correct mean and variance");
    if (abs(mean-0.5) < 0.01 && abs(variance-1.0/12) < 0.01)
      printf("passed\n");
    else
      printf("FAILED mean %f variance %f\n", mean, variance);
    }     
  }

  void nn_test_15_dropout()
  {
    printf("TEST 15. DROPOUT\n");
    Dataset dataset;
    read_MNIST_dataset(base_path + "../../resources/MNIST-dataset", &dataset);
    train_test_split(&dataset, 0.1);

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<FlattenLayer>(28,28,1));
    nn2.add_layer(std::make_shared<DenseLayer>(28*28*1, 200), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DropoutLayer>(0.2));
    nn2.add_layer(std::make_shared<DenseLayer>(200, 200), Initializer::He);
    nn2.add_layer(std::make_shared<ReLULayer>());
    nn2.add_layer(std::make_shared<DenseLayer>(200, 10), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());

    nn2.train(dataset.train_data, dataset.train_labels, 128, 10000, OptimizerAdam(0.001f), Loss::CrossEntropy);

    std::vector<float> y_res(dataset.test_labels.size(),0);
    nn2.evaluate(dataset.test_data, y_res);
    float acc = 0.0f;
    float cnt = 0.0f;
    for (int i=0;i<dataset.test_labels.size();i+=10)
    {
      int max_pos = 0;
      int ref_max_pos = 0;
      for (int j=0;j<10;j++)
      {
        if (y_res[i+j] > y_res[i+max_pos])
          max_pos = j;
        if (dataset.test_labels[i+j] > dataset.test_labels[i+ref_max_pos])
          ref_max_pos = j;
      }

      acc += (max_pos == ref_max_pos);
      cnt++;
    }
    float error_rate = 1 - acc/cnt;
    printf(" 15.1. %-64s", "Error rate <10% ");
    if (error_rate < 0.1f)
      printf("passed %f\n",error_rate);
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void nn_test_17_Leaky_ReLU_classifier()
  {
    printf("TEST 17. CLASSIFICATION WITH Leaky RELU\n");
    int sz = 25000;
    int dim = 10;
    std::vector<float> X(dim*sz,0);
    std::vector<float> res(2*sz,0);
    for (int i=0;i<sz;i++)
    {
      float s = 0;
      for (int j=0;j<dim;j++)
      {
        X[dim*i + j] = 2*((double)rand())/RAND_MAX - 1;
        s += X[dim*i + j];
      }
      res[2*i+0] = s > 0;
      res[2*i+1] = s <= 0;
    }

    std::vector<float> X_test(dim*sz,0);
    std::vector<float> res_test(2*sz,0);
    for (int i=0;i<sz;i++)
    {
      float s = 0;
      for (int j=0;j<dim;j++)
      {
        X_test[dim*i + j] = 2*((double)rand())/RAND_MAX - 1;
        s += X_test[dim*i + j];
      }
      res_test[2*i+0] = s > 0;
      res_test[2*i+1] = s <= 0;
    }

    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<DenseLayer>(dim, 64), Initializer::He);
    nn2.add_layer(std::make_shared<LeakyReLULayer>(0.01f));
    nn2.add_layer(std::make_shared<DenseLayer>(64, 64), Initializer::He);
    nn2.add_layer(std::make_shared<LeakyReLULayer>(0.01f));
    nn2.add_layer(std::make_shared<DenseLayer>(64, 2), Initializer::He);
    nn2.add_layer(std::make_shared<SoftMaxLayer>());
    nn2.train(X, res, 256, 5000, OptimizerAdam(0.01f), Loss::CrossEntropy);

    std::vector<float> y(2*sz,0);
    nn2.evaluate(X_test, y);
    float diff = 0.0f;
    for (int i=0;i<2*sz;i++)
      diff += (y[i]>0.5) != (res_test[i]>0.5);

    float error_rate = diff/sz;
    printf(" 16.1. %-64s", "Error rate <3% ");
    if (error_rate < 0.03f)
      printf("passed\n");
    else
      printf("FAILED, error rate %f\n", error_rate);
  }

  void tp_test_14_conv3D()
  {
    printf("TEST 14. 3D convolution\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(5,5,5);
      TensorToken B = TensorToken(3,3,3);
      TensorToken res_1 = TensorToken::conv3D(A, B);
      TensorToken res_2 = TensorToken::conv3D(A, B, 2);

      tc.input(A, "A");
      tc.input(B, "B");
      tc.output(res_1, "res_1");
      tc.output(res_2, "res_2");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = { 1, 0, 2, 0, 3,
                             0, 0, 0, 0, 0,
                             4, 0, 5, 0, 6,
                             0, 0, 0, 0, 0,
                             7, 0, 8, 0, 9,

                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,

                            10, 0,11, 0,12,
                             0, 0, 0, 0, 0,
                            13, 0,14, 0,15,
                             0, 0, 0, 0, 0,
                            16, 0,17, 0,18,

                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,

                            19, 0,20, 0,21,
                             0, 0, 0, 0, 0,
                            22, 0,23, 0,24,
                             0, 0, 0, 0, 0,
                            25, 0,26, 0,27,};

    //A.resize(125*2);
    //for (int i=0;i<125;i++)
    //  A[i+125] = -A[i];

    std::vector<float> B = {1,0,1,
                            0,1,0,
                            1,0,1,
                            
                            0,1,0,
                            1,0,1,
                            0,1,0,

                            1,0,1,
                            0,1,0,
                            1,0,1,};

    std::vector<float> res_1(27, 0.0f);
    std::vector<float> res_1_ref = {60, 0, 68, 0, 19, 0, 84, 0, 92, 0, 25, 0, 27, 0, 29, 0, 31, 0, 132, 0, 140, 0, 37, 0, 156, 0, 164};
    std::vector<float> res_2(8, 0.0f);
    std::vector<float> res_2_ref = {60, 68, 84, 92, 132, 140, 156, 164};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("B", B.data(), B.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("res_1", res_1.data(), res_1.size());
    TensorProcessor::get_output("res_2", res_2.data(), res_2.size());
    //for (auto &v : res_1)
    //  printf("%f, ", v);
    //printf("\n");
    //for (auto &v : res_2)
    //  printf("%f, ", v);
    //printf("\n");

    {
    float diff = 0.0f;
    for (int i=0;i<res_1.size();i++)
      diff += abs(res_1[i] - res_1_ref[i]);
    
    printf(" 14.1. %-64s","stride = 1");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<res_2.size();i++)
      diff += abs(res_2[i] - res_2_ref[i]);
    
    printf(" 14.2. %-64s","stride = 2");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }   
  }

  void nn_test_18_conv3D_backward()
  {
    printf("TEST 18. CONV_3D BACKWARD PASS\n");
    std::vector<float> X = { 1, 0, 2, 0, 3,
                             0, 0, 0, 0, 0,
                             4, 0, 5, 0, 6,
                             0, 0, 0, 0, 0,
                             7, 0, 8, 0, 9,

                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,

                            10, 0,11, 0,12,
                             0, 0, 0, 0, 0,
                            13, 0,14, 0,15,
                             0, 0, 0, 0, 0,
                            16, 0,17, 0,18,

                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,

                            19, 0,20, 0,21,
                             0, 0, 0, 0, 0,
                            22, 0,23, 0,24,
                             0, 0, 0, 0, 0,
                            25, 0,26, 0,27,};
    std::vector<float> r(27, 0.0f);
    std::vector<float> r_ref = {60, 0, 68, 0, 19, 0, 84, 0, 92, 0, 25, 0, 27, 0, 29, 0, 31, 0, 132, 0, 140, 0, 37, 0, 156, 0, 164};
    NeuralNetwork nn2;
    nn2.add_layer(std::make_shared<Conv3DLayer>(5,5,5,1,1), Initializer::GlorotNormal);
    nn2.train(X, r_ref, 1, 5000, OptimizerAdam(0.005f), Loss::MSE, false);
    nn2.evaluate(X, r);
    //for (auto &v : r)
    //  printf("%f, ", v);
    //printf("\n");
  
    float diff = 0.0f;
    for (int i=0;i<r_ref.size();i++)
      diff += std::abs(r[i] - r_ref[i]);
    
    printf(" 18.1. %-64s", "Correct result ");
    if (diff < 1.0f)
      printf("passed\n");
    else
      printf("FAILED diff %f > %f\n", diff, 1.0f);
  }

  void tp_test_15_max_pooling_3D()
  {
    printf("TEST 15. MAX POOLING 3D\n");

    TensorCompiler tc;
    {
      tc.start_program();
      TensorToken A = TensorToken(4, 4, 2);
      TensorToken B = TensorToken(2, 2, 1);
      TensorToken::issue_command(TensorProgram::MPOOL_3D, A, A, B, 2, 2, 2);
      TensorToken dOut = TensorToken(2, 2, 1);
      TensorToken dIn = TensorToken(4, 4, 2);
      TensorToken::issue_command(TensorProgram::MPOOL_3D_D, A, dOut, dIn, 2, 2, 2);
      tc.input(A, "A");
      tc.input(dOut, "dOut");
      tc.output(B, "B");
      tc.output(dIn, "dIn");
    }
    TensorProgram p = tc.finish_program();

    std::vector<float> A = {1,0, -5,-6,
                            0,0, -7,-8,

                            3,4,  6, 7,
                            7,0,  8, 9,
                              
                            1,0, -2,-1,
                            0,0, -3,-4,

                            3,4,  6, 7,
                            0,0,  9, 9};
    std::vector<float> dOut = {1,2,3,4};
    std::vector<float> B(2*2, 0.0f), B_ref = {1,-1,7,9};
    std::vector<float> dIn(4*4*2, 0.0f), dIn_ref = {1,0, 0,0,
                                                    0,0, 0,0,

                                                    0,0, 0,0,
                                                    3,0, 0,4,
                                                    
                                                    0,0, 0,2,
                                                    0,0, 0,0,

                                                    0,0, 0,0,
                                                    0,0, 0,0};

    TensorProcessor::set_program(p);
    TensorProcessor::set_input("A", A.data(), A.size());
    TensorProcessor::set_input("dOut", dOut.data(), dOut.size());
    TensorProcessor::execute();
    TensorProcessor::get_output("B", B.data(), B.size());
    TensorProcessor::get_output("dIn", dIn.data(), dIn.size());
    //for (auto &v : B)
    //  printf("%f, ", v);
    //printf("\n");
    //for (auto &v : dIn)
    //  printf("%f, ", v);
    //printf("\n");
    {
    float diff = 0.0f;
    for (int i=0;i<B.size();i++)
      diff += abs(B[i] - B_ref[i]);
    
    printf(" 15.1. %-64s","forward pass");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
    {
    float diff = 0.0f;
    for (int i=0;i<dIn.size();i++)
      diff += abs(dIn[i] - dIn_ref[i]);
    
    printf(" 15.2. %-64s","backward pass");
    if (diff < 1e-6)
      printf("passed\n");
    else
      printf("FAILED diff %f >= %f\n", diff, 1e-6f);
    }
  }

  void perform_tests_tensor_processor(const std::vector<int> &test_ids)
  {
    srand(time(NULL));
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      tp_test_1_tensor_processor,  
      tp_test_2_tensor_tokens,
      tp_test_3_tensor_operations,
      tp_test_4_aliases,
      tp_test_5_softmax,
      tp_test_6_logic_operations,
      tp_test_7_padding,
      tp_test_8_conv2D,
      tp_test_9_conv2D_blur,
      tp_test_10_random,
      tp_test_11_max_pooling,
      tp_test_12_arithmetics_benchmark,
      tp_test_13_dilation,
      tp_test_14_conv3D,
      tp_test_15_max_pooling_3D,
    };

    if (tests.empty())
    {
      tests.resize(test_functions.size());
      for (int i=0;i<test_functions.size();i++)
        tests[i] = i+1;
    }

    TensorProcessor::init("CPU");

    for (int i=0;i<80;i++)
      printf("#");
    printf("\nTENSOR PROCESSOR TESTS\n");
    for (int i=0;i<80;i++)
      printf("#");
    printf("\n");
    
    for (int i : tests)
    {
      assert(i > 0 && i <= test_functions.size());
      test_functions[i-1]();
    }
  }

  void perform_tests_tensor_processor_GPU(const std::vector<int> &test_ids)
  {
#if defined(USE_GPU)
    srand(time(NULL));
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      tp_test_1_tensor_processor,  
      tp_test_2_tensor_tokens,
      tp_test_3_tensor_operations,
      tp_test_4_aliases,
      tp_test_5_softmax,
      tp_test_6_logic_operations,
      tp_test_7_padding,
      tp_test_8_conv2D,
      tp_test_9_conv2D_blur,
      tp_test_10_random,
      tp_test_11_max_pooling,
      tp_test_12_arithmetics_benchmark,
      tp_test_13_dilation,
      tp_test_14_conv3D,
      tp_test_15_max_pooling_3D,
    };

    if (tests.empty())
    {
      tests.resize(test_functions.size());
      for (int i=0;i<test_functions.size();i++)
        tests[i] = i+1;
    }

    TensorProcessor::init("GPU");

    for (int i=0;i<80;i++)
      printf("#");
    printf("\nGPU TENSOR PROCESSOR TESTS\n");
    for (int i=0;i<80;i++)
      printf("#");
    printf("\n");
    
    for (int i : tests)
    {
      assert(i > 0 && i <= test_functions.size());
      test_functions[i-1]();
    }
#else
    for (int i=0;i<80;i++)
      printf("#");
    printf("\nGPU TENSOR PROCESSOR TESTS: GPU is not supported. Enable VULKAN module\n");
    for (int i=0;i<80;i++)
      printf("#");
    printf("\n");
#endif
  }

  void perform_tests_neural_networks(const std::vector<int> &test_ids)
  {
    srand(time(NULL));
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      nn_test_1_linear_regression,
      nn_test_2_linear_regression_train,
      nn_test_3_SIREN_image,
      nn_test_4_SIREN_SDF,
      nn_test_5_simple_classifier,
      nn_test_6_conv2D_forward,
      nn_test_7_conv2D_backward,
      nn_test_8_conv2D_no_padding,
      nn_test_9_synthetic_images_classifier,
      nn_test_10_MNIST,
      nn_test_11_CIFAR10,
      nn_test_12_conv2D_stride,
      nn_test_13_binary_classification_metrics,
      nn_test_14_batch_normalization,
      nn_test_15_dropout,
      nn_test_16_ReLU_classifier,
      nn_test_17_Leaky_ReLU_classifier,
      nn_test_18_conv3D_backward,
    };

    if (tests.empty())
    {
      tests.resize(test_functions.size());
      for (int i=0;i<test_functions.size();i++)
        tests[i] = i+1;
    }
    
    TensorProcessor::init("GPU");

    for (int i=0;i<80;i++)
      printf("#");
    printf("\nNEURAL NETWORKS TESTS\n");
    for (int i=0;i<80;i++)
      printf("#");
    printf("\n");
    
    for (int i : tests)
    {
      assert(i > 0 && i <= test_functions.size());
      test_functions[i-1]();
    }
  }

  void perform_tests()
  {
    perform_tests_tensor_processor({});
    perform_tests_tensor_processor_GPU({});
    perform_tests_neural_networks({});
  }
}