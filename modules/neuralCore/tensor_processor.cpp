#include "tensor_processor_impl.h"
#include <cassert>
#include <cstring>

#define USE_GPU 0
#if USE_GPU
#include "tensor_processor_impl_gpu.h"
std::shared_ptr<TensorProcessorImpl> CreateTensorProcessorImpl_GPU(vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated);
#endif

namespace nn
{
  TensorProcessor::TensorProcessor()
  {
    #if USE_GPU
    pImpl = CreateTensorProcessorImpl_GPU(vk_utils::globalContextGet(), 256);
    #else
    pImpl = std::shared_ptr<TensorProcessorImpl>(new TensorProcessorImpl());
    #endif
  }

  void TensorProcessor::set_program(const TensorProgram &_program)
  {
    program = _program;
    cpu_memory = std::vector<float>(program.total_memory_req, 0);
    program_prepared = true;
  }

  void TensorProcessor::set_input(const std::map<std::string, float * const> &vars)
  {
    assert(vars.size() == program.input_vars.size());

    for (auto &p : program.input_vars)
    {
      float *data = vars.at(p.first);
      unsigned offset = program.vars[p.second].offset;
      unsigned size = program.vars[p.second].total_size;
      memcpy(cpu_memory.data() + offset, data, sizeof(float) * size);
    }
    
    input_prepared = true;
  }

  void TensorProcessor::set_output(const std::map<std::string, float *> &vars)
  {
    assert(vars.size() == program.output_vars.size());
    for (auto &p : program.output_vars)
    {
      float *data = vars.at(p.first);
      unsigned offset = program.vars[p.second].offset;
      unsigned size = program.vars[p.second].total_size;
      memcpy(data, cpu_memory.data() + offset, sizeof(float) * size);
    }
  }

  void TensorProcessor::execute()
  {
    assert((pImpl.get() != nullptr) && input_prepared && program_prepared);
    pImpl->process(program, cpu_memory.data(), cpu_memory.data(), cpu_memory.size());
  }
}