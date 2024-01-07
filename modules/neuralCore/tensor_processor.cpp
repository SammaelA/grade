#include "tensor_processor_impl.h"
#include <cassert>
#include <cstring>

#if defined(USE_GPU)
#include "tensor_processor_impl_gpu.h"
std::shared_ptr<TensorProcessorImpl> CreateTensorProcessorImpl_GPU(vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated);
#endif

namespace nn
{

  std::vector<TensorProgram::CmdProperties> TensorProgram::cmd_properties =
  {
    {NOOP     , "NOOP"     , AUXILIARY     , SELF_APPLICABLE_NO },
    {FTT      , "FTT"      , AUXILIARY     , SELF_APPLICABLE_NO },

    {MOV      , "MOV"      , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {FILL     , "FILL"     , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {COPY     , "COPY"     , MEM_MANAGEMENT, SELF_APPLICABLE_NO },

    {ADD      , "ADD"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {SUB      , "SUB"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {MUL      , "MUL"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {DIV      , "DIV"      , ARITHMETICS   , SELF_APPLICABLE_YES},

    {EXP      , "EXP"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {POW      , "POW"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {SIN      , "SIN"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {COS      , "COS"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {LOG      , "LOG"      , ELEMENTWISE   , SELF_APPLICABLE_YES},

    {SUM      , "SUM"      , REDUCTION     , SELF_APPLICABLE_NO },
    {O_SUM    , "O_SUM"    , REDUCTION     , SELF_APPLICABLE_NO },

    {MATMUL_T , "MATMUL_T" , ALGEBRA       , SELF_APPLICABLE_NO },
    {TRANSP   , "TRANSP"   , ALGEBRA       , SELF_APPLICABLE_NO },
    {OUTER_P  , "OUTER_P"  , ALGEBRA       , SELF_APPLICABLE_NO },
    {OUTER_PS , "OUTER_PS" , ALGEBRA       , SELF_APPLICABLE_NO },

    {URAND    , "URAND"    , OTHER         , SELF_APPLICABLE_NO },
  };

  std::unique_ptr<TensorProcessor> proc;
  void TensorProcessor::init(std::string backend)
  {
    assert(TensorProgram::cmd_properties.size() == TensorProgram::CMD_COUNT);
    for (int i=0;i<TensorProgram::cmd_properties.size();i++)
      assert(TensorProgram::cmd_properties[i].type == i);
    if (!proc)
      proc.reset(new TensorProcessor());
    #if defined(USE_GPU)
    if (backend == "GPU")
      proc->pImpl = CreateTensorProcessorImpl_GPU(vk_utils::globalContextGet(false, 0u), 256);
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
    proc->pImpl->allocate_memory(proc->program.total_memory_req);
    proc->pImpl->CommitDeviceData();
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
      proc->pImpl->set_input(data, offset, size);
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
      proc->pImpl->get_output(data, offset, size);
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
    proc->pImpl->process(proc->program);
  }
}