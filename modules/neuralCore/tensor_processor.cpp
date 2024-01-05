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
  std::unique_ptr<TensorProcessor> proc;
  void TensorProcessor::init(std::string backend)
  {
    if (!proc)
      proc.reset(new TensorProcessor());
    #if USE_GPU
    if (backend == "GPU")
      proc->pImpl = CreateTensorProcessorImpl_GPU(vk_utils::globalContextGet(), 256);
    else
    #endif
      proc->pImpl = std::shared_ptr<TensorProcessorImpl>(new TensorProcessorImpl());
  }
  TensorProcessor::TensorProcessor()
  {

  }

  void TensorProcessor::set_program(const TensorProgram &_program)
  {
    if (!proc)
      TensorProcessor::init();
    proc->program = _program;
    proc->cpu_memory = std::vector<float>(proc->program.total_memory_req, 0);
    proc->program_prepared = true;

    proc->input_prepared = {};
    for (auto &p : proc->program.input_vars)
      proc->input_prepared[p.first] = false;
  }

  void TensorProcessor::set_input(const std::string &name, float * const data, unsigned data_size)
  {
    if (proc->input_prepared.find(name) == proc->input_prepared.end())
    {
      printf("TensorProcessor::trying to set input variable %s that does not exist!\n",name.c_str());
    }
    else
    {
      unsigned var_id = proc->program.input_vars.at(name);
      unsigned offset = proc->program.vars[var_id].offset;
      unsigned size = std::min(data_size, proc->program.vars[var_id].total_size);
      memcpy(proc->cpu_memory.data() + offset, data, sizeof(float) * size);
      proc->input_prepared[name] = true;
    }
  }

  void TensorProcessor::get_output(const std::string &name, float * data, unsigned data_size)
  {
    if (proc->program.output_vars.find(name) == proc->program.output_vars.end())
    {
      printf("TensorProcessor::trying to get output variable %s that does not exist!\n",name.c_str());
    }
    else
    {
      unsigned var_id = proc->program.output_vars.at(name);
      unsigned offset = proc->program.vars[var_id].offset;
      unsigned size = std::min(data_size, proc->program.vars[var_id].total_size);
      memcpy(data, proc->cpu_memory.data() + offset, sizeof(float) * size);
    }
  }

  void TensorProcessor::execute()
  {
    assert((proc->pImpl.get() != nullptr) && proc->program_prepared);
    for (auto &p : proc->input_prepared)
    {
      if (!p.second)
      {
        printf("TensorProcessor::input variable %s is not set! Execution failed!\n", p.first.c_str());
      }
    }
    proc->pImpl->process(proc->program, proc->cpu_memory.data(), proc->cpu_memory.data(), proc->cpu_memory.size());
  }
}