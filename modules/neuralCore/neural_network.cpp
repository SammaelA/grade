#include "neural_network.h"
#include <cassert>
#include <cstdio>
#include <string>
#include <fstream>
#include <chrono>
#include <cstring>
#include <random>
#include <algorithm>

constexpr bool DEBUG = false;
namespace nn
{
  void DenseLayer::init()
  {
    weights.clear();

    weights.push_back(TensorToken(output_shape[0], input_shape[0])); // A
    weights.push_back(TensorToken(output_shape[0]));                 // b

    dLoss_dWeights.resize(weights.size());
  }

  TensorToken DenseLayer::forward(const TensorToken &in)
  {
    return TensorToken::mat_mul_t(in, weights[0].transpose()) + weights[1];
  }

  TensorToken DenseLayer::backward(const TensorToken &input, const TensorToken &output, const TensorToken &dLoss_dOutput)
  {
    float batch_size = (float)(input.sizes[1]);

    dLoss_dWeights[0] = TensorToken::mat_mul_t(input.transpose(), dLoss_dOutput.transpose())/batch_size;
    dLoss_dWeights[1] = dLoss_dOutput.outer_sum()/batch_size;

    return TensorToken::mat_mul_t(dLoss_dOutput, weights[0]);
  }

  unsigned total_size(const std::vector<unsigned> &sizes)
  {
    unsigned total_sz = 1;
    for (auto sz : sizes)
      total_sz *= sz;
    return total_sz;
  }

  std::string time_pretty_str(float t_ms)
  {
    char str[16];
    float t_s = t_ms/1000;
    float t_m = t_s/60;
    float t_h = t_m/60;
    if (t_ms < 1000)
      snprintf(str, 16, "%.1f ms", t_ms);
    else if (t_s < 60)
      snprintf(str, 16, "%.1f s", t_s);
    else if (t_m < 60)
      snprintf(str, 16, "%.1f m", t_m);
    else
      snprintf(str, 16, "%.1f h", t_h);
    return std::string(str);
  }

  void zero_initialization(float *data, int size)
  {
    std::fill_n(data, size, 0.0f);
  }

  void he_initialization(float *data, int size, int fan_in, int fan_out)
  {
    float mn = -sqrt(6.0 / fan_in);
    float mx = sqrt(6.0 / fan_in);
    float d = mx - mn;
    for (int i = 0; i < size; i++)
      data[i] = d * (((double)rand()) / RAND_MAX) + mn;
  }
  void SIREN_initialization(float *data, int size, int fan_in, int fan_out)
  {
    constexpr float omega_0 = 30.0f;
    float mx = (fan_in == 2) ? 0.5 : (sqrt(6.0 / fan_in) / omega_0);
    float mn = -mx;
    float d = mx - mn;
    for (int i = 0; i < size; i++)
      data[i] = d * (((double)rand()) / RAND_MAX) + mn;
  }

  void Glorot_normal_initialization(float *data, int size, int fan_in, int fan_out)
  {
    double stddev = sqrt(2.0f/(fan_in+fan_out));
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> distr{0.0, stddev};
    for (int i = 0; i < size; i++)
      data[i] = distr(gen);
  }

  void BatchNorm_initialization(float *data, int size)
  {
    assert(size % 2 == 0);
    for (int i = 0; i < size/2; i++)
      data[i] = 1;
    for (int i = size/2; i < size; i++)
      data[i] = 0;
  }

  void NeuralNetwork::add_layer(std::shared_ptr<Layer> layer, Initializer initializer)
  {
    layers.push_back(layer);
    initializers.push_back(initializer);
  }

  void NeuralNetwork::set_batch_size_for_evaluate(int size)
  {
    batch_size_evaluate = size;
    if (initialized)
      get_evaluate_prog();
  }

  bool NeuralNetwork::check_validity()
  {
    if (layers[0]->input_shape.empty() || layers[0]->output_shape.empty())
    {
      printf("NeuralNetwork: first layer must have implicit shape!\n");
      return false;
    }
    for (int i = 1; i < layers.size(); i++)
    {
      // shape-insensitive layers
      if (layers[i]->input_shape.empty() && layers[i]->output_shape.empty())
      {
        layers[i]->input_shape = layers[i-1]->output_shape;
        layers[i]->output_shape = layers[i-1]->output_shape;
      }
      else if (layers[i]->input_shape.size() != layers[i - 1]->output_shape.size())
      {
        printf("NeuralNetwork: layers %d and %d have incompatible shapes!\n", i - 1, i);
        return false;
      }
    }
    return true;
  }

  void NeuralNetwork::initialize()
  {
    if (!check_validity())
      return;
    
    for (auto &l : layers)
      total_params += l->parameters_count();
    weights.resize(total_params, 0.0f);

    int offset = 0;
    for (int i=0;i<layers.size();i++)
    {
      unsigned fan_in = total_size(layers[i]->input_shape);
      unsigned fan_out = total_size(layers[i]->output_shape);
      unsigned size = layers[i]->parameters_count();

      if (size == 0)
        continue;

      switch (initializers[i])
      {
        case Initializer::Zero:
          zero_initialization(weights.data()+offset, size);
          break;
        case Initializer::He:
          he_initialization(weights.data()+offset, size, fan_in, fan_out);
          break;
        case Initializer::Siren:
          SIREN_initialization(weights.data()+offset, size, fan_in, fan_out);
          break;
        case Initializer::GlorotNormal:
          Glorot_normal_initialization(weights.data()+offset, size, fan_in, fan_out);
          break;
        case Initializer::BatchNorm:
          assert(dynamic_cast<BatchNormLayer*>(layers[i].get()));
          BatchNorm_initialization(weights.data()+offset, size);
          break;
        default:
          break;
      }
      offset += size;
    }

    get_evaluate_prog();
    //print_info();
    initialized = true;
  }

  void NeuralNetwork::initialize_with_weights(const float *w)
  {
    initialize();
    weights = std::vector<float>(w, w + total_params);
  }

  void NeuralNetwork::initialize_from_file(std::string filename)
  {
    initialize();
    std::ifstream in(filename, std::ios_base::binary);
    assert(in.is_open());
    in.read(reinterpret_cast<char*>(weights.data()), sizeof(float)*weights.size());
    in.close();
  }

  void NeuralNetwork::save_weights_to_file(std::string filename)
  {
    std::ofstream out(filename, std::ios_base::binary);
    assert(out.is_open());
    out.write(reinterpret_cast<char*>(weights.data()), sizeof(float)*weights.size());
    out.close();
  }

  void NeuralNetwork::set_arch_to_file(std::string filename)
  {
    std::ofstream out(filename);
    assert(out.is_open());
    for (auto &l : layers)
    {
      out << l->get_name() << " ";
      if (l->input_shape.size() > 0)
      {
        out << "input shape ("<<l->input_shape[0];
        for (int i=1;i<l->input_shape.size();i++)
          out << ", "<<l->input_shape[i];
        out <<") ";
      }
      if (l->output_shape.size() > 0)
      {
        out << "output shape ("<<l->output_shape[0];
        for (int i=1;i<l->output_shape.size();i++)
          out << ", "<<l->output_shape[i];
        out <<")";
      }
      out<<"\n";
    }
    out.close();
  }

  void NeuralNetwork::print_info()
  {
    printf("Neural Network\n");
    printf("%d layers\n", (int)(layers.size()));
    for (int i = 0; i < layers.size(); i++)
      printf("Layer %d has %d parameters\n", i, layers[i]->parameters_count());
    printf("%d input size\n", total_size(layers[0]->input_shape));
    printf("%d output size\n", total_size(layers.back()->output_shape));
    printf("%d weights\n", total_params);
  }

  void NeuralNetwork::get_evaluate_prog()
  {
    for (auto &l : layers)
      l->training_mode = true;

    TensorCompiler compiler;
    compiler.start_program();
    auto i_shape = layers[0]->input_shape;
    i_shape.push_back(batch_size_evaluate);
    auto o_shape = layers.back()->output_shape;
    o_shape.push_back(batch_size_evaluate);

    TensorToken input = TensorToken(i_shape);
    TensorToken output = TensorToken(o_shape);
    TensorToken w = TensorToken(total_params);
    unsigned offset = 0;
    for (auto &l : layers)
    {
      l->init();
      for (auto &lw : l->weights)
      {
        unsigned sz = lw.total_size();
        lw.copy_to({0, sz}, w, {offset, offset + sz});
        offset += sz;
      }
    }

    TensorToken t = input;
    for (auto &l : layers)
      t = l->forward(t);
    output = t;
    
    compiler.inout(w, "W");
    compiler.input(input, "In");
    compiler.output(output, "Out");
    evaluate_prog = compiler.finish_program();
  }

  TensorProgram NeuralNetwork::get_train_prog(int batch_size, Optimizer optimizer, Loss loss)
  {
    for (auto &l : layers)
      l->training_mode = true;

    TensorCompiler compiler;
    compiler.start_program();
    auto i_shape = layers[0]->input_shape;
    i_shape.push_back(batch_size);
    auto o_shape = layers.back()->output_shape;
    o_shape.push_back(batch_size);

    TensorToken input = TensorToken(i_shape); compiler.input(input, "In");
    TensorToken target_output = TensorToken(o_shape); compiler.input(target_output, "Out");
    TensorToken w = TensorToken(total_params); compiler.inout(w, "W");
    
    unsigned offset = 0;
    for (auto &l : layers)
    {
      l->init();
      for (int i=0;i<l->weights.size();i++)
      {
        unsigned sz = l->weights[i].total_size();
        l->weights[i].copy_to({0, sz}, w, {offset, offset + sz});
        offset += sz;
      }
    }

    std::vector<TensorToken> all_outputs;
    all_outputs.push_back(layers[0]->forward(input));
        TensorToken t = input;
    for (int i=1;i<layers.size();i++)
      all_outputs.push_back(layers[i]->forward(all_outputs.back()));
  
    TensorToken output = all_outputs.back();

    //loss
    TensorToken l, dLoss_dOutput;
    if (loss == Loss::MSE)
    {
      TensorToken diff = output - target_output;
      l = (diff*diff).sum()/(float)(output.total_size());
      dLoss_dOutput = 2.0f*diff;
    }
    else if (loss == Loss::CrossEntropy)
    {
      TensorToken mo = -1.0f*target_output;
      l = (mo * TensorToken::log(output + 1e-15f)).sum()/(float)(batch_size);
      dLoss_dOutput = mo / (output + 1e-15f);
    }

    for (int i=layers.size()-1;i>0;i--)
      dLoss_dOutput = layers[i]->backward(all_outputs[i-1], all_outputs[i], dLoss_dOutput);
    layers[0]->backward(input, all_outputs[0], dLoss_dOutput);
  
    TensorToken grad = TensorToken(total_params);
    offset = 0;
    for (auto &l : layers)
    {
      for (auto &dLoss_dWeight : l->dLoss_dWeights)
      {
        unsigned sz = dLoss_dWeight.total_size();
        grad.set({offset, offset + sz}, dLoss_dWeight.flatten());
        offset += sz;
      }
    }

    //Adam optimizer
    TensorToken V = TensorToken(total_params); compiler.inout(V, "V");
    TensorToken S = TensorToken(total_params); compiler.inout(S, "S");
    TensorToken iter = TensorToken(1); compiler.input(iter, "iter");
    TensorToken learning_rate = TensorToken(1); compiler.input(learning_rate, "learning_rate");

    if (std::holds_alternative<OptimizerGD>(optimizer))
    {
      OptimizerGD opt = std::get<OptimizerGD>(optimizer);

      w -= learning_rate*grad;
    }
    else if (std::holds_alternative<OptimizerAdam>(optimizer))
    {
      OptimizerAdam opt = std::get<OptimizerAdam>(optimizer);

      V = opt.beta_1*V + (1.0f - opt.beta_1)*grad;
      TensorToken Vh = V / (1.0f - TensorToken::pow(opt.beta_1, iter + 1.0f));
      S = opt.beta_2*S + (1.0f - opt.beta_2)*grad*grad;
      TensorToken Sh = S / (1.0f - TensorToken::pow(opt.beta_2, iter + 1.0f));
      w -= learning_rate*Vh/(TensorToken::sqrt(Sh) + opt.eps);
    }
    else if (std::holds_alternative<OptimizerRMSProp>(optimizer)) 
    {
      OptimizerRMSProp opt = std::get<OptimizerRMSProp>(optimizer);

      S = opt.beta*S + (1.0f - opt.beta)*grad*grad;
      w -= learning_rate*grad/(TensorToken::sqrt(S) + opt.eps);
    }
    else if (std::holds_alternative<OptimizerMomentum>(optimizer)) 
    {
      OptimizerMomentum opt = std::get<OptimizerMomentum>(optimizer);

      V = opt.momentum*S + (1.0f - opt.momentum)*grad;
      w -= learning_rate*V;      
    }
    
    compiler.output(l, "loss");
    if (DEBUG)
      compiler.output(grad, "grad");
    return compiler.finish_program();
  }

  void NeuralNetwork::evaluate(std::vector<float> &input_data, std::vector<float> &output_data, int samples)
  {
    unsigned input_size = total_size(layers[0]->input_shape);
    if (samples < 0)
      samples = input_data.size()/input_size;
    evaluate(input_data.data(), output_data.data(), samples);
  }

  void NeuralNetwork::evaluate(const float *input_data, float *output_labels, int samples)
  {
    unsigned input_size = total_size(layers[0]->input_shape);
    unsigned output_size = total_size(layers.back()->output_shape);

    unsigned batches = (samples + batch_size_evaluate - 1)/batch_size_evaluate;
    
    TensorProcessor::set_program(evaluate_prog);
    TensorProcessor::set_input("W", weights.data(), weights.size());
    
    for (int i=0;i<batches;i++)
    {
      TensorProcessor::set_input("In", input_data + i*batch_size_evaluate*input_size, samples*input_size - i*batch_size_evaluate*input_size);
      TensorProcessor::execute();
      TensorProcessor::get_output("Out", output_labels + i*batch_size_evaluate*output_size, samples*output_size - i*batch_size_evaluate*output_size);
    }
  }

  void NeuralNetwork::train(const std::vector<float> &inputs /*[input_size, count]*/, const std::vector<float> &outputs /*[output_size, count]*/,
                             int batch_size, int iterations, Optimizer optimizer, Loss loss, bool verbose)
  {
    unsigned input_size = total_size(layers[0]->input_shape);
    unsigned count = inputs.size()/input_size;
    train(inputs.data(), outputs.data(), count, batch_size, ceil(batch_size*iterations/(float)count), false, optimizer, loss, Metric::Accuracy, verbose);
  }

  void NeuralNetwork::train(const float *data, const float *labels, int samples, int batch_size, int epochs, bool use_validation, Optimizer optimizer, 
                            Loss loss, Metric metric, bool verbose)
  {
    initialize();

    TensorProgram train_prog = get_train_prog(batch_size, optimizer, loss);

    unsigned input_size = total_size(layers[0]->input_shape);
    unsigned output_size = total_size(layers.back()->output_shape);
    float validation_frac = use_validation ? 0.1f : 0.0f;
    unsigned valid_count = samples*validation_frac;
    unsigned count = samples - valid_count;
    unsigned iters_per_epoch = std::max(1u, count/batch_size);
    unsigned iterations = epochs * iters_per_epoch;
    unsigned iters_per_validation = std::max(100u, iters_per_epoch);

    float start_learning_rate = 0;
    float end_learning_rate = 0;
    if (std::holds_alternative<OptimizerGD>(optimizer))
    {
      OptimizerGD opt = std::get<OptimizerGD>(optimizer);
      start_learning_rate = opt.learning_rate;
      end_learning_rate = opt.learning_rate;
    }
    else if (std::holds_alternative<OptimizerAdam>(optimizer))
    {
      OptimizerAdam opt = std::get<OptimizerAdam>(optimizer);
      start_learning_rate = opt.learning_rate;
      end_learning_rate = opt.minimum_learning_rate;
    }
    else if (std::holds_alternative<OptimizerRMSProp>(optimizer)) 
    {
      OptimizerRMSProp opt = std::get<OptimizerRMSProp>(optimizer);
      start_learning_rate = opt.learning_rate;
      end_learning_rate = opt.minimum_learning_rate;
    }
    else if (std::holds_alternative<OptimizerMomentum>(optimizer)) 
    {
      OptimizerMomentum opt = std::get<OptimizerMomentum>(optimizer);
      start_learning_rate = opt.learning_rate;
      end_learning_rate = opt.learning_rate;    
    }
    else
    {
      printf("NeuralNetwork: unknown optimizer!!!\n");
      assert(false);
    }

    std::vector<float> V(total_params, 0);
    std::vector<float> S(total_params, 0);
    std::vector<float> in_batch(input_size*batch_size);
    std::vector<float> out_batch(output_size*batch_size);
    std::vector<float> validation_labels(output_size*valid_count);

    std::vector<float> best_weights = weights;
    float best_metric = (metric == Metric::MSE || metric == Metric::MAE) ? 1e9 : -1e9;

    TensorProcessor::set_program(train_prog);
    TensorProcessor::set_input("W", weights.data(), weights.size());
    TensorProcessor::set_input("V", V.data(), V.size());
    TensorProcessor::set_input("S", S.data(), S.size());

    if (verbose)
      printf("started training %u iterations %d epochs\n", iterations, epochs);
    float av_loss = 0;
    auto t_prev = std::chrono::steady_clock::now();

    for (int it=0;it<iterations;it++)
    {
      for (int i=0;i<batch_size;i++)
      {
        unsigned b_id = rand()%count;
        memcpy(in_batch.data() + i*input_size, data + b_id*input_size, sizeof(float)*input_size);
        memcpy(out_batch.data() + i*output_size, labels + b_id*output_size, sizeof(float)*output_size);
      }

      float iter = it;
      float r = it/(float)iterations;
      float learning_rate = (1-r)*start_learning_rate + r*end_learning_rate;
      TensorProcessor::set_input("In", in_batch.data(), in_batch.size());
      TensorProcessor::set_input("Out", out_batch.data(), out_batch.size());
      TensorProcessor::set_input("iter", &iter, 1);
      TensorProcessor::set_input("learning_rate", &learning_rate, 1);
      TensorProcessor::execute();
      float loss = -1;
      TensorProcessor::get_output("loss", &loss, 1);
      av_loss += loss;
      if (it > 0 && it % iters_per_validation == 0)
      {
        if (use_validation)
        {
          TensorProcessor::get_output("W", weights.data(), weights.size());
          TensorProcessor::get_output("V", V.data(), V.size());
          TensorProcessor::get_output("S", S.data(), S.size());
          evaluate(data + count*input_size, validation_labels.data(), valid_count);
          TensorProcessor::set_program(train_prog);
          TensorProcessor::set_input("W", weights.data(), weights.size());
          TensorProcessor::set_input("V", V.data(), V.size());
          TensorProcessor::set_input("S", S.data(), S.size());

          float m = calculate_metric(validation_labels.data(), labels + count*output_size, valid_count, metric);
          if (( (metric == Metric::MSE || metric == Metric::MAE) && m <= best_metric) ||
              (!(metric == Metric::MSE || metric == Metric::MAE) && m >= best_metric))
          {
            best_metric = m;
            memcpy(best_weights.data(), weights.data(), sizeof(float)*weights.size());
          }

          if (verbose)
            printf("[%d/%d] Loss = %f Metric = %f ", it/iters_per_epoch, iterations/iters_per_epoch, av_loss/iters_per_validation, m);
        }
        else if (verbose)
          printf("[%d/%d] Loss = %f ", it/iters_per_epoch, iterations/iters_per_epoch, av_loss/iters_per_validation);
        
        auto t = std::chrono::steady_clock::now();               
        double ms = 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t - t_prev).count()/iters_per_validation;
        t_prev = t;
        if (verbose)
          printf("%s/epoch, ETA: %s\n", time_pretty_str(ms*iters_per_epoch).c_str(), time_pretty_str(ms*(iterations-it)).c_str());

        av_loss = 0;
      }

      if (DEBUG)
      {
        std::vector<float> grad(weights.size(),0);
        TensorProcessor::get_output("grad", grad.data(), grad.size());
        TensorProcessor::get_output("W", weights.data(), weights.size());
        printf("grad = [ ");
        for (int i=0;i<weights.size();i++)
          printf("%f ", grad[i]);
        printf("]\n");

        printf("w = [ ");
        for (int i=0;i<weights.size();i++)
          printf("%f ", weights[i]);
        printf("]\n");
      }
    }
    TensorProcessor::get_output("W", weights.data(), weights.size());
  }

  void get_confusion_matrix(const float *output, const float *output_ref, int samples, float threshold, float out_confusion_matrix[4])
  {
    for (int i=0;i<4;i++)
      out_confusion_matrix[i]=0;
    for (int i=0;i<samples;i++)
    {
      out_confusion_matrix[0] += output_ref[2*i+0] >= threshold && output[2*i+0] >= threshold; //true  positive
      out_confusion_matrix[1] += output_ref[2*i+0] < threshold  && output[2*i+0] >= threshold; //false positive
      out_confusion_matrix[2] += output_ref[2*i+0] >= threshold && output[2*i+0] < threshold;  //false negative
      out_confusion_matrix[3] += output_ref[2*i+0] < threshold  && output[2*i+0] < threshold;  //true  negative
    }
    for (int i=0;i<4;i++)
      out_confusion_matrix[i] /= samples;
  }

  float NeuralNetwork::calculate_metric(const float *output, const float *output_ref, int samples, Metric metric)
  {
    unsigned output_size = total_size(layers.back()->output_shape);
    if (metric == Metric::MSE || metric == Metric::MAE)
    {
      //regression metrics
      double res = 0;
      if (metric == Metric::MSE)
      {
        #pragma omp parallel for reduction(+:res)
        for (int i=0;i<output_size*samples;i++)
          res += (output[i] - output_ref[i])*(output[i] - output_ref[i]);
      }
      else if (metric == Metric::MAE)
      {
        #pragma omp parallel for reduction(+:res)
        for (int i=0;i<output_size*samples;i++)
          res += std::abs(output[i] - output_ref[i]);        
      }
      return res/(output_size*samples);
    }
    else if (metric == Metric::Accuracy)
    {
      float thr = 0;
      int right_answers = 0;

      #pragma omp parallel for reduction(+:right_answers)
      for (int i=0;i<samples;i++)
      {
        int ref_class = 0;
        for (int j=0;j<output_size;j++)
          if (output_ref[i*output_size + j] > output_ref[i*output_size + ref_class])
            ref_class = j;
        
        int pred_class = 0;
        for (int j=0;j<output_size;j++)
          if (output[i*output_size + j] > output[i*output_size + pred_class])
            pred_class = j;
        
        if (ref_class == pred_class)
          right_answers += 1;
      }
      return right_answers/(float)samples;
    }
    else
    {
      //metrics for binary classification
      assert(output_size == 2);
      float confusion_matrix[4];
      get_confusion_matrix(output, output_ref, samples, 0.5, confusion_matrix);
      //printf("conf %f %f %f %f\n", confusion_matrix[0], confusion_matrix[1], confusion_matrix[2], confusion_matrix[3]);
      if (metric == Metric::Precision)
        return confusion_matrix[0]/(confusion_matrix[0]+confusion_matrix[1]);
      else if (metric == Metric::Recall)
        return confusion_matrix[0]/(confusion_matrix[0]+confusion_matrix[2]);
      else if (metric == Metric::AUC_ROC)
      {
        constexpr unsigned steps = 1000;
        std::vector<std::pair<float, float>> tpr_fpr(steps+2);
        tpr_fpr[0] = {0,0};
        for (int i=0;i<steps;i++)
        {
          get_confusion_matrix(output, output_ref, samples, i/(float)steps, confusion_matrix);
          tpr_fpr[i+1] = std::pair<float, float>(confusion_matrix[0]/(confusion_matrix[0]+confusion_matrix[2] + 1e-9), 
                                                 confusion_matrix[1]/(confusion_matrix[1]+confusion_matrix[3] + 1e-9));
        }
        tpr_fpr[steps+1] = {1,1};

        std::sort(tpr_fpr.begin(), tpr_fpr.end(), [](std::pair<float, float> a, std::pair<float, float> b)-> bool { 
          return a.second==b.second ? a.first<b.first : a.second<b.second;});

        float auc_roc = 0;
        for (int i=0;i<=steps;i++)
        {
          auc_roc += 0.5*(tpr_fpr[i+1].first + tpr_fpr[i].first)*(tpr_fpr[i+1].second-tpr_fpr[i].second);
          //printf("%f %f %f\n", tpr_fpr[i].first, tpr_fpr[i].second, auc_roc);
        }
        
        return auc_roc;
      }
      else if (metric == Metric::AUC_PR)
      {
        constexpr unsigned steps = 1000;
        std::vector<std::pair<float, float>> rec_pr(steps+2);
        rec_pr[0] = {0,1};
        for (int i=0;i<steps;i++)
        {
          get_confusion_matrix(output, output_ref, samples, i/(float)steps, confusion_matrix);
          rec_pr[i+1] = std::pair<float, float>(confusion_matrix[0]/(confusion_matrix[0]+confusion_matrix[2] + 1e-9), 
                                                 confusion_matrix[0]/(confusion_matrix[0]+confusion_matrix[1] + 1e-9));
        }
        rec_pr[steps+1] = {1,0};

        std::sort(rec_pr.begin(), rec_pr.end(), [](std::pair<float, float> a, std::pair<float, float> b)-> bool { return a.second<b.second;});

        float auc_pr = 0;
        for (int i=0;i<=steps;i++)
        {
          auc_pr += 0.5*(rec_pr[i+1].first + rec_pr[i].first)*(rec_pr[i+1].second-rec_pr[i].second);
          //printf("%f %f %f\n", rec_pr[i].first, rec_pr[i].second, auc_pr);
        }
        
        return auc_pr;
      }
    }
    return 0;
  }
}