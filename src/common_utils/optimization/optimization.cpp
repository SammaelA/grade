#include "optimization.h"
#include "common_utils/distribution.h"

namespace opt
{
  float get_average(const std::vector<float> &X)
  {
    float av = 0;
    for (auto &v : X)
      av += v;
    av /= X.size();
    return av;
  }

  float get_dev(const std::vector<float> &X)
  {
    float av = get_average(X);
    float dev = 0;
    for (auto &v : X)
      dev += (v - av)*(v - av);
    dev = sqrt(dev)/(X.size()-1);
    return dev;
  }
  
  void check_gradients(opt_func_with_grad &F, const std::vector<float> &min_X, const std::vector<float> &max_X,
                      int samples, float h, bool show_detailed_info)
  {
    for (int i=0;i<min_X.size();i++)
    {
      std::vector<float> grad_num;
      std::vector<float> grad_gt;
      std::vector<float> errors;
      float grad_sum = 1e-9;
      for (int j=0;j<samples;j++)
      {
        std::vector<float> X;
        for (int k=0;k<min_X.size();k++)
          X.push_back(urand(min_X[k] + h, max_X[k] - h));
        
        auto res = F(X);

        X[i] -= h; //F(X - h)
        float y_1 = F(X).first;
        X[i] += 2*h;//F(X + h)
        float y_2 = F(X).first;
        float dF_dxj = (y_2 - y_1)/(2*h);
        grad_num.push_back(dF_dxj);
        grad_gt.push_back(res.second[i]);
        errors.push_back(abs(dF_dxj - res.second[i]));
        grad_sum += 0.5*abs(dF_dxj) + 0.5*abs(res.second[i]);
      }

      float e_sum = 0;
      float e_max = errors[0];
      for (float &e : errors)
      {
        if (e > e_max)
          e_max = e;
        e_sum += e;
      }
      e_sum /= samples;
      grad_sum /= samples;

      int misdirections = 0;
      int zero_grads_gt = 0;
      int zero_grads_num = 0;

      for (int k=0; k<samples; k++)
      {
        int s1 = grad_gt[k] > 1e-6 ? 1 : (grad_gt[k] < -1e-6 ? -1 : 0);
        int s2 = grad_num[k] > 1e-6 ? 1 : (grad_num[k] < -1e-6 ? -1 : 0);
        misdirections += (s1 != s2);
        zero_grads_gt += abs(grad_gt[k]) < 1e-6;
        zero_grads_num += abs(grad_num[k]) < 1e-6;
      }

      debug("Param %d\n", i);
      debug("GT derivative Average: %.4f, Std dev: %.4f\n", get_average(grad_gt), get_dev(grad_gt));
      debug("Num derivative Average: %.4f, Std dev: %.4f\n", get_average(grad_num), get_dev(grad_num));
      debug("Average error: %.5f (%.2f%), Max: %.5f (%.2f%)\n", e_sum, 100*(e_sum/grad_sum), e_max, 100*(e_max/grad_sum));
      debug("%d/%d wrong derivative signs\n", misdirections, samples);
      debug("%d/%d and %d/%d zero derivatives\n", zero_grads_gt, samples, zero_grads_num, samples);
    }
  }
}