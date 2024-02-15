#pragma once
#include <vector>
#include <string>

namespace nn
{
  void perform_tests();
  void perform_tests_tensor_processor(const std::vector<int> &test_ids);
  void perform_tests_tensor_processor_GPU(const std::vector<int> &test_ids);
  void perform_tests_neural_networks(const std::vector<int> &test_ids);

  extern std::string base_path;
}