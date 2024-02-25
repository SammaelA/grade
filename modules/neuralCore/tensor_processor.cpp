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

    {MOV      , "MOV"      , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {FILL     , "FILL"     , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {COPY     , "COPY"     , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {PAD      , "PAD"      , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {FLIP     , "FLIP"     , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {DILATE   , "DILATE"   , MEM_MANAGEMENT, SELF_APPLICABLE_NO },
    {URAND    , "URAND"    , MEM_MANAGEMENT, SELF_APPLICABLE_NO },

    {ADD      , "ADD"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {SUB      , "SUB"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {MUL      , "MUL"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {DIV      , "DIV"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {GREATER  , "GREATER"  , ARITHMETICS   , SELF_APPLICABLE_YES},
    {LESS     , "LESS"     , ARITHMETICS   , SELF_APPLICABLE_YES},
    {EQUAL    , "EQUAL"    , ARITHMETICS   , SELF_APPLICABLE_YES},
    {GE       , "GE"       , ARITHMETICS   , SELF_APPLICABLE_YES},
    {LE       , "LE"       , ARITHMETICS   , SELF_APPLICABLE_YES},
    {NE       , "NE"       , ARITHMETICS   , SELF_APPLICABLE_YES},
    {OR       , "OR"       , ARITHMETICS   , SELF_APPLICABLE_YES},
    {AND      , "AND"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {WHERE    , "WHERE"    , ARITHMETICS   , SELF_APPLICABLE_YES},
    {MIN      , "MIN"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {MAX      , "MAX"      , ARITHMETICS   , SELF_APPLICABLE_YES},
    {POW      , "POW"      , ARITHMETICS   , SELF_APPLICABLE_YES},

    {EXP      , "EXP"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {SQRT     , "SQRT"     , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {SIN      , "SIN"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {COS      , "COS"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {LOG      , "LOG"      , ELEMENTWISE   , SELF_APPLICABLE_YES},
    {NOT      , "NOT"      , ELEMENTWISE   , SELF_APPLICABLE_YES},

    {SUM      , "SUM"      , REDUCTION     , SELF_APPLICABLE_NO },
    {O_SUM    , "O_SUM"    , REDUCTION     , SELF_APPLICABLE_NO },
    {MINIMUM  , "MINIMUM"  , REDUCTION     , SELF_APPLICABLE_NO },
    {MAXIMUM  , "MAXIMUM"  , REDUCTION     , SELF_APPLICABLE_NO },

    {MATMUL_T , "MATMUL_T" , ALGEBRA       , SELF_APPLICABLE_NO },
    {TRANSP   , "TRANSP"   , ALGEBRA       , SELF_APPLICABLE_NO },
    {OUTER_P  , "OUTER_P"  , ALGEBRA       , SELF_APPLICABLE_NO },
    {SMAX_D   , "SMAX_D"   , ALGEBRA       , SELF_APPLICABLE_NO },
    {CONV_2D  , "CONV_2D"  , ALGEBRA       , SELF_APPLICABLE_NO },
    {MPOOL    ,  "MPOOL"   , ALGEBRA       , SELF_APPLICABLE_NO },
    {MPOOL_D  , "MPOOL_D"  , ALGEBRA       , SELF_APPLICABLE_NO },
    {CONV_3D  , "CONV_3D"  , ALGEBRA       , SELF_APPLICABLE_NO },
    {MPOOL_3D ,  "MPOOL_3D", ALGEBRA       , SELF_APPLICABLE_NO },
    {MPOOL_3D_D, "MPOOL_3D_D", ALGEBRA       , SELF_APPLICABLE_NO },
  };

  std::unique_ptr<TensorProcessor> proc;
  void TensorProcessor::init(std::string backend)
  {
    assert(TensorProgram::cmd_properties.size() == TensorProgram::CMD_COUNT);
    for (int i=0;i<TensorProgram::cmd_properties.size();i++)
      assert(TensorProgram::cmd_properties[i].type == i);
    if (!proc)
    {
      proc.reset(new TensorProcessor());
      proc->used_backend = backend;
    }
    #if defined(USE_GPU)
    if (proc->used_backend == "GPU")
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
    TensorProcessor::init();
    proc->program = _program;
    proc->pImpl->allocate_memory(proc->program.total_memory_req);
    proc->pImpl->CommitDeviceData();
    proc->program_prepared = true;

    proc->input_prepared = {};
    for (auto &p : proc->program.input_vars)
      proc->input_prepared[p.first] = false;

    if (proc->program.constants.size() > 0)
      proc->pImpl->set_input(proc->program.constants.data(), 
                             proc->program.vars[TensorProgram::CONSTS_VAR_ID].offset,
                             proc->program.vars[TensorProgram::CONSTS_VAR_ID].total_size);

    _stat_execution_times = 0;
    _stat_time_cmd_id  = std::vector<float>(nn::TensorProgram::cmd_properties.size(), 0.0f);
    _stat_time_cmd_num = std::vector<float>(_program.commands.size(), 0.0f);
  
  }

  void TensorProcessor::set_input(const std::string &name, const float * const data, unsigned data_size)
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

  void TensorProcessor::print_execution_stat()
  {
    float total_time = 0.0;
    for (auto &v : _stat_time_cmd_id)
      total_time += v;
    total_time *= 1e-6;//us to seconds
    printf("TensorProcessor: program execution statistics\n");
    printf("Executions: %d\n", _stat_execution_times);
    printf("Took %.2f s (%.4f s/exec)\n", total_time, total_time/_stat_execution_times);
    printf("Time spent by command type:\n");
    for (int i=0;i<TensorProgram::CMD_COUNT;i++)
      printf("%8s: %.3f s (%4.1f%%)\n", TensorProgram::cmd_properties[i].name.c_str(), 1e-6*_stat_time_cmd_id[i],
                                          100*1e-6*_stat_time_cmd_id[i]/total_time);
    printf("Time spent by command number:\n");
    for (int i=0;i<proc->program.commands.size();i++)
      printf("%3d[%8s]: %.3f s (%4.1f%%)\n", i, TensorProgram::cmd_properties[proc->program.commands[i].type].name.c_str(), 
                                                  1e-6*_stat_time_cmd_num[i],
                                                  100*1e-6*_stat_time_cmd_num[i]/total_time);
    printf("\n");
  }
}