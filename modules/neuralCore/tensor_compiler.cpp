#include "tensor_compiler.h"
#include <string>
#include <cstring>

namespace nn
{
  TensorCompiler *TensorToken::tp = nullptr;

  unsigned TensorCompiler::add_var(const TensorToken &t)
  {
    vars.emplace_back();

    vars.back().Dim = t.Dim;
    for (int i = 0; i < t.Dim; i++)
      vars.back().sizes[i] = t.sizes[i];
    vars.back().total_size = 1;
    for (int i = 0; i < t.Dim; i++)
      vars.back().total_size *= t.sizes[i];

    return vars.size() - 1;
  }
  void TensorCompiler::ftt(unsigned id, float val)
  {
    commands.push_back({TensorProgram::FTT, {0, 0, id, *((unsigned *)(&val))}});
  }

  void TensorCompiler::input(const TensorToken &t, std::string name)
  {
    assert(input_vars.find(name) == input_vars.end());
    input_vars[name] = t.id;
    vars[t.id].is_input = true;
  }

  void TensorCompiler::output(const TensorToken &t, std::string name)
  {
    assert(output_vars.find(name) == output_vars.end());
    output_vars[name] = t.id;
    vars[t.id].is_output = true;
  }

  void TensorCompiler::start_program()
  {
    TensorToken::tp = this;
    printf("started recording tensor program\n");
    add_command(TensorProgram::NOOP);
    vars.emplace_back();
  }

  void TensorCompiler::compactify()
  {
    unsigned used_vars = 1; //var[0] is fixed
    std::vector<unsigned> var_remap(vars.size(), 0);
    for (unsigned i=1;i<vars.size();i++)
    {
      if (vars[i].is_input || vars[i].is_output)
        var_remap[i] = used_vars++;
    }

    for (unsigned i=1;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      unsigned A = commands[i].args[0];
      unsigned B = commands[i].args[1];
      unsigned C = commands[i].args[2];
      if (A > 0 && var_remap[A] == 0)
        var_remap[A] = used_vars++;
      if (B > 0 && var_remap[B] == 0)
        var_remap[B] = used_vars++;
      if (C > 0 && var_remap[C] == 0)
        var_remap[C] = used_vars++;
    }

    std::vector<TensorProgram::Command> new_commands;
    new_commands.push_back(commands[0]);
    for (unsigned i=1;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      new_commands.push_back({commands[i].type, {var_remap[commands[i].args[0]], var_remap[commands[i].args[1]], var_remap[commands[i].args[2]], 
                                                 commands[i].args[3], commands[i].args[4], commands[i].args[5]}});
    }
    commands = new_commands;

    std::vector<Variable> new_vars(used_vars);
    new_vars[0] = vars[0];
    for (unsigned i=1;i<vars.size();i++)
    {
      if (var_remap[i] > 0)
      {
        new_vars[var_remap[i]] = vars[i];
      }
    }
    vars = new_vars;

    for (auto &p : input_vars)
      p.second = var_remap[p.second];
    for (auto &p : output_vars)
      p.second = var_remap[p.second];
  }

  bool TensorCompiler::optimize_mov_cycle()
  {
    std::vector<unsigned> fm(vars.size(), commands.size());
    std::vector<unsigned> fu(vars.size(), commands.size());
    std::vector<unsigned> lm(vars.size(), 0);
    std::vector<unsigned> lu(vars.size(), 0);

    for (unsigned i=0;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      unsigned A = commands[i].args[0];
      unsigned B = commands[i].args[1];
      unsigned C = commands[i].args[2];

      fu[A] = std::min(fu[A], i);
      fu[B] = std::min(fu[B], i);
      fm[C] = std::min(fm[C], i);

      lu[A] = std::max(lu[A], i);
      lu[B] = std::max(lu[B], i);
      lm[C] = std::max(lm[C], i);
    }

    for (unsigned i=0;i<vars.size();i++)
    {
      if (vars[i].is_input)
        fm[i] = 0;
      if (vars[i].is_output)
        lu[i] = commands.size();
      //printf("%u - %u %u - %u %u\n",i,fu[i],lu[i],fm[i],lm[i]);
    }

    bool optimized = false;
    for (unsigned i=0;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      unsigned A = commands[i].args[0];
      unsigned B = commands[i].args[1];
      unsigned C = commands[i].args[2];

      if ((A > 0 && fm[A] >= i) || (B > 0 && fm[B] >= i) || (C > 0 && lu[C] <= i))
      {
        optimized = true;
        commands[i].type = TensorProgram::NOOP;
      }
    }
    return optimized;
  }

  void TensorCompiler::optimize_program()
  {
    //initial optimization
    //removes redundant MOV operations
    while (optimize_mov_cycle());
    compactify();
  }

  unsigned TensorCompiler::calculate_memory_layout()
  {
    unsigned total_memory = 0;
    for (auto &var : vars)
    {
      var.offset = total_memory;
      total_memory += var.total_size;
    }

    return total_memory;
  }

  TensorProgram TensorCompiler::finish_program()
  {
    std::vector<std::string> names = {"NOOP", "ADD", "MUL", "DIV", "EXP", "SUM", "MATMUL_T",
                                      "MOV", "FTT", "COPY", "TRANSP", "OUTER_P",
                                      "", "", "", "", "",
                                      "", "", "", "", ""};


    optimize_program();
    unsigned total_memory_req = calculate_memory_layout();

    TensorProgram pr;
    pr.commands = commands;
    pr.output_vars = output_vars;
    pr.input_vars = input_vars;
    pr.total_memory_req = total_memory_req;
    pr.vars = std::vector<TensorProgram::Variable>(vars.size());
    for (int i=0; i<vars.size(); i++)
    {
      pr.vars[i].Dim = vars[i].Dim;
      pr.vars[i].offset = vars[i].offset;
      memcpy(pr.vars[i].sizes, vars[i].sizes, sizeof(pr.vars[i].sizes));
      pr.vars[i].total_size = vars[i].total_size;
    }

    printf("finished recording tensor program\n");
    printf("requires %d bytes of memory\n", (int)(sizeof(float)*total_memory_req));

    printf("%d variables\n", (int)vars.size());
    for (unsigned vid = 0; vid < vars.size(); vid++)
    {
      auto &var = vars[vid];
      printf("V%2u: %u %5u [%2u %2u %2u %2u] %s %s\n", vid, var.Dim, var.total_size, var.sizes[0], var.sizes[1], var.sizes[2], var.sizes[3],
             var.is_input ? "  input" : "", var.is_output ? "output" : "");
    }

    printf("%d commands\n", (int)commands.size());
    unsigned cid = 0;
    for (auto &cmd : commands)
    {
      printf("Cmd %2u: %-8s %2u %2u %2u - %3u %3u %3u\n", cid, names[cmd.type].c_str(), cmd.args[0], cmd.args[1], cmd.args[2],
             cmd.args[3], cmd.args[4], cmd.args[5]);
      cid++;
    }

    return pr;
  }

  void TensorCompiler::add_command(TensorProgram::CommandType type, unsigned A, unsigned B, unsigned C, unsigned arg0, unsigned arg1, unsigned arg2)
  {
    commands.push_back({type, {A, B, C, arg0, arg1, arg2}});
  }
}