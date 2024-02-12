#include "tensor_compiler.h"
#include <string>
#include <cstring>
#include <map>
#include <algorithm>

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
    unsigned constant_id = constants.size();
    for (int i=0;i<constants.size();i++)
    {
      //to prevent any rounding errors we merge only exactly equal constants
      if (val == constants[i])
      {
        constant_id = i;
        break;
      }
    }
    if (constant_id == constants.size())
    {
      constants.push_back(val);
      vars[TensorProgram::CONSTS_VAR_ID].total_size = constants.size();
      vars[TensorProgram::CONSTS_VAR_ID].sizes[0] = constants.size();
    }
    commands.push_back({TensorProgram::COPY, {TensorProgram::CONSTS_VAR_ID, 0, id, constant_id, 0, 1u}});
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

  void TensorCompiler::inout(const TensorToken &t, std::string name)
  {
    input(t, name);
    output(t, name);
  }

  void TensorCompiler::start_program()
  {
    TensorToken::tp = this;

    vars.clear();
    commands.clear();
    constants.clear();
    input_vars.clear();
    output_vars.clear();

    vars.emplace_back();//leave 0 variable empty to let 0 indicate no argument in command
    vars.emplace_back();//reserve first variable for array of all constants
    vars[TensorProgram::CONSTS_VAR_ID].Dim = 1;
    vars[TensorProgram::CONSTS_VAR_ID].is_input = true;
    vars[TensorProgram::CONSTS_VAR_ID].is_output = true;
    add_command(TensorProgram::NOOP);
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
                                                 commands[i].args[3], commands[i].args[4], commands[i].args[5], commands[i].args[6], commands[i].args[7]}});
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

  void TensorCompiler::remove_noop()
  {
    std::vector<TensorProgram::Command> new_commands;
    for (unsigned i=1;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      new_commands.push_back(commands[i]);
    }
    commands = new_commands;
  }

  void TensorCompiler::calculate_variable_usage_intervals()
  {
    fm = std::vector<unsigned>(vars.size(), commands.size());
    fu = std::vector<unsigned>(vars.size(), commands.size());
    lm = std::vector<unsigned>(vars.size(), 0);
    lu = std::vector<unsigned>(vars.size(), 0);
    
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
  }

  bool TensorCompiler::optimize_unused_cycle()
  {
    calculate_variable_usage_intervals();

    //remove commands that use uninitialized variables or 
    //create variables that won't be used
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

  bool TensorCompiler::have_same_scheme(const Variable &a, const Variable &b)
  {
    bool same_scheme = a.Dim == b.Dim;
    for (int i=0;i<a.Dim;i++)
      same_scheme = same_scheme && a.sizes[i] == b.sizes[i];
    return same_scheme;
  }

  bool TensorCompiler::is_self_applicable_command(TensorProgram::CommandType type)
  {
    return TensorProgram::cmd_properties[type].is_self_applicable == TensorProgram::SELF_APPLICABLE_YES;
  }

  void TensorCompiler::replace_output_var(unsigned old_id, unsigned new_id)
  {
    if (!vars[old_id].is_output)
      return;
    vars[old_id].is_output = false;
    vars[new_id].is_output = true;
    for (auto &p : output_vars)
      if (p.second == old_id)
        p.second = new_id;
  }

  void TensorCompiler::reset_alias_rec(unsigned alias_id, unsigned master_id, unsigned base_offset)
  {
    for (auto &child_alias : vars[alias_id].aliases)
    {
      reset_alias_rec(child_alias, master_id, base_offset + vars[child_alias].alias_range_from);
      vars[child_alias].alias_master_id = master_id;
      vars[child_alias].alias_range_from += base_offset;
      vars[child_alias].alias_range_to += base_offset;
      vars[master_id].aliases.push_back(child_alias);
    }
    vars[alias_id].aliases = {};
  }

  void TensorCompiler::set_alias(unsigned alias_id, unsigned master_id, unsigned from, unsigned to)
  {
    assert(alias_id > 0 && master_id > 0);
    assert(vars[alias_id].is_alias == false);

    vars[alias_id].is_alias = true;
    vars[alias_id].alias_master_id = master_id;
    vars[alias_id].alias_range_from = from;
    vars[alias_id].alias_range_to = to;

    vars[master_id].aliases.push_back(alias_id);
    reset_alias_rec(alias_id, master_id, vars[alias_id].alias_range_from);
  }

  bool TensorCompiler::optimize_renaming_moves()
  {
    bool made_change = false;
    for (unsigned i=1;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      unsigned A = commands[i].args[0];
      unsigned C = commands[i].args[2];

      bool is_mov = (commands[i].type == TensorProgram::MOV);
      bool is_full_copy = (commands[i].type == TensorProgram::COPY && commands[i].args[3] == 0 && 
                           commands[i].args[4] == 0 && commands[i].args[5] == vars[A].total_size);
      bool last_A = lu[A] <= i && lm[A] <= i;
      bool first_C = fu[C] >= i && fm[C] >= i;
      bool same_scheme = have_same_scheme(vars[A], vars[C]);
      if ((is_mov || is_full_copy) && last_A && same_scheme)
      {
        //this command only renames variable, no operation is needed
        commands[i].type = TensorProgram::NOOP;
        //we use A instead of C everywhere
        lu[A] = std::max(lu[A], lu[C]);
        lm[A] = std::max(lm[A], lm[C]);
        //if C is output, me make A output variable instead
        replace_output_var(C, A);
        for (unsigned j=i+1;j<commands.size();j++)
        {
          if (commands[j].args[0] == C)
            commands[j].args[0] = A;
          if (commands[j].args[1] == C)
            commands[j].args[1] = A;
          if (commands[j].args[2] == C)
            commands[j].args[2] = A;
        }
        made_change = true;
      }
      else if (is_mov)
      {
        //replace move with copy for uniformity
        commands[i].type = TensorProgram::COPY;
        commands[i].args[3] = 0;
        commands[i].args[4] = 0;
        commands[i].args[5] = vars[A].total_size;
      }
    }
    return made_change;
  }

  bool TensorCompiler::optimize_self_applicable_commands()
  {
    bool made_change = false;
    for (unsigned i=1;i<commands.size();i++)
    {
      if (commands[i].type == TensorProgram::NOOP)
        continue;
      unsigned A = commands[i].args[0];
      unsigned B = commands[i].args[1];
      unsigned C = commands[i].args[2];

      bool last_A = lu[A] <= i && lm[A] <= i;
      bool last_B = lu[B] <= i && lm[B] <= i;
      bool first_C = fu[C] >= i && fm[C] >= i;

      unsigned rp = 0;
      if (is_self_applicable_command(commands[i].type))
      {
        if (A > 0 && last_A && first_C && have_same_scheme(vars[A], vars[C]))
          rp = A; //replace C = A x B with A = A x B
        else if (B > 0 && last_B && first_C && have_same_scheme(vars[B], vars[C]))
          rp = B; //replace C = A x B with B = A x B
      }
      if (rp > 0)
      {
        commands[i].args[2] = rp;
        //we use rp instead of C everywhere
        lu[rp] = std::max(lu[rp], lu[C]);
        lm[rp] = std::max(lm[rp], lm[C]);
        //if C is output, me make rp output variable instead
        replace_output_var(C, rp);
        for (unsigned j=i+1;j<commands.size();j++)
        {
          if (commands[j].args[0] == C)
            commands[j].args[0] = rp;
          if (commands[j].args[1] == C)
            commands[j].args[1] = rp;
          if (commands[j].args[2] == C)
            commands[j].args[2] = rp;
        }
        made_change = true;
      }
    }
    return made_change;
  }

  void TensorCompiler::optimize_copy_to_aliases()
  {
    struct Usage
    {
      unsigned cmd_id;
      unsigned begin;
      unsigned end;
    };

    std::vector<std::vector<Usage>> usages, modifications;
    usages.resize(vars.size());
    modifications.resize(vars.size());

    for (unsigned i=1;i<vars.size();i++)
      if (vars[i].is_input)
        modifications[i].push_back({0, 0, vars[i].total_size});

    for (unsigned i=1;i<commands.size();i++)
    {
      unsigned A = commands[i].args[0];
      unsigned B = commands[i].args[1];
      unsigned C = commands[i].args[2];

      if (commands[i].type == TensorProgram::COPY)
      {
        unsigned A_begin = commands[i].args[3];
        unsigned A_end = A_begin + commands[i].args[5];
        unsigned C_begin = commands[i].args[4];
        unsigned C_end = C_begin + commands[i].args[5];

        usages[A].push_back({i, A_begin, A_end});
        modifications[C].push_back({i, C_begin, C_end});
      }
      else
      {
        if (A > 0)
          usages[A].push_back({i, 0, vars[A].total_size});
        if (B > 0 && B != A)
          usages[B].push_back({i, 0, vars[B].total_size});
        if (C > 0)
          modifications[C].push_back({i, 0, vars[C].total_size});
      }
    }

    for (unsigned i=1;i<vars.size();i++)
      if (vars[i].is_output)
        usages[i].push_back({(unsigned)commands.size()+1, 0, vars[i].total_size});

    auto nextM = [&](unsigned var_id, unsigned start_index, unsigned begin, unsigned end) -> unsigned
    {
      for (auto &m : modifications[var_id])
        if (m.cmd_id > start_index && std::max(begin, m.begin) < std::min(end, m.end))
          return m.cmd_id;
      return commands.size()+2;
    };
    auto nextMA = [&](unsigned var_id, unsigned start_index, unsigned begin, unsigned end) -> unsigned
    {
      unsigned minM = nextM(var_id, start_index, begin, end);
      for (auto &aid : vars[var_id].aliases)
      {
        if (std::max(begin, vars[aid].alias_range_from) < std::min(end, vars[aid].alias_range_to))
          minM = std::min(minM, nextM(aid, start_index,
                          std::max(0,(int)begin-(int)vars[aid].alias_range_from),
                          std::max(0,(int)end-(int)vars[aid].alias_range_from)));
      }
      return minM;
    };

    auto nextU = [&](unsigned var_id, unsigned start_index, unsigned begin, unsigned end) -> unsigned
    {
      for (auto &m : usages[var_id])
        if (m.cmd_id > start_index && std::max(begin, m.begin) < std::min(end, m.end))
          return m.cmd_id;
      return commands.size()+1;
    };
    auto nextUA = [&](unsigned var_id, unsigned start_index, unsigned begin, unsigned end) -> unsigned
    {
      unsigned minU = nextU(var_id, start_index, begin, end);
      for (auto &aid : vars[var_id].aliases)
      {
        if (std::max(begin, vars[aid].alias_range_from) < std::min(end, vars[aid].alias_range_to))
          minU = std::min(minU, nextU(aid, start_index,
                          std::max(0,(int)begin-(int)vars[aid].alias_range_from),
                          std::max(0,(int)end-(int)vars[aid].alias_range_from)));
      }
      return minU;
    };


    auto prevU = [&](unsigned var_id, unsigned start_index, unsigned begin, unsigned end) -> unsigned
    {
      for (int i=usages[var_id].size()-1;i>=0;i--)
      {
        auto &m = usages[var_id][i];
        if (m.cmd_id < start_index && std::max(begin, m.begin) < std::min(end, m.end))
          return m.cmd_id;
      }
      return 0;
    };
    auto prevUA = [&](unsigned var_id, unsigned start_index, unsigned begin, unsigned end) -> unsigned
    {
      unsigned maxU = prevU(var_id, start_index, begin, end);
      for (auto &aid : vars[var_id].aliases)
      {
        if (std::max(begin, vars[aid].alias_range_from) < std::min(end, vars[aid].alias_range_to))
          maxU = std::max(maxU, prevU(aid, start_index,
                          std::max(0,(int)begin-(int)vars[aid].alias_range_from),
                          std::max(0,(int)end-(int)vars[aid].alias_range_from)));
      }
      return maxU;
    };

    auto set_aliases = [&](bool output_pass)
    {
    for (unsigned i=1;i<commands.size();i++)
    {
      if (commands[i].type != TensorProgram::COPY)
        continue;
      
      unsigned A = commands[i].args[0];
      unsigned A_begin = commands[i].args[3];
      unsigned A_end = A_begin + commands[i].args[5];
      unsigned C = commands[i].args[2];
      unsigned C_begin = commands[i].args[4];
      unsigned C_end = C_begin + commands[i].args[5];

      //printf("COPY%u %u(%u) %u(%u) - %u %u %u\n", i, A,vars[A].total_size, C,vars[C].total_size, commands[i].args[3], commands[i].args[4], commands[i].args[5]);
      //C is created with this copy
      if (commands[i].args[5] == vars[C].total_size && modifications[C][0].cmd_id == i && vars[C].is_alias == false)
      {
        unsigned master_id = A;
        unsigned from = A_begin;
        unsigned to = A_end;
        while (vars[master_id].is_alias)
        {
          from += vars[master_id].alias_range_from;
          to += vars[master_id].alias_range_from;
          master_id = vars[master_id].alias_master_id;
        }

        //when A and C coexist, C and (part of A that was copied to C)
        //both have the same value. So we can avoid this copy
        if (vars[C].is_output == output_pass &&
            nextUA(master_id, i, from, to) < nextMA(C, i, C_begin, C_end) && 
            nextUA(C, i, C_begin, C_end) < nextMA(master_id, i, from, to))
        {
          commands[i].type = TensorProgram::NOOP;
          set_alias(C, master_id, from, to);
          //printf("%u can be an alias of %u (%u %u) (%u %u)\n", C, master_id, nextUA(A, i, from, to), nextMA(C, i, C_begin, C_end),
          //       nextUA(C, i, C_begin, C_end), nextMA(A, i, from, to));
        }
      }
      else if (vars[A].is_output == output_pass && //A is no longer used after this copy
               commands[i].args[5] == vars[A].total_size && usages[A].back().cmd_id == i && vars[A].is_alias == false)
      {
        unsigned master_id = C;
        unsigned from = C_begin;
        unsigned to = C_end;
        while (vars[master_id].is_alias)
        {
          from += vars[master_id].alias_range_from;
          to += vars[master_id].alias_range_from;
          master_id = vars[master_id].alias_master_id;
        }

        if (prevUA(master_id, i, from, to) < modifications[A][0].cmd_id)
        {
          commands[i].type = TensorProgram::NOOP;
          set_alias(A, C, from, to);
          //printf("2 %u can be an alias of %u\n", A, master_id);
        }
        else if (prevUA(master_id, i, from, to) == modifications[A][0].cmd_id && 
                 is_self_applicable_command(commands[modifications[A][0].cmd_id].type))
        {
          printf("modific %s\n", TensorProgram::cmd_properties[commands[modifications[A][0].cmd_id].type].name.c_str());
          commands[i].type = TensorProgram::NOOP;
          set_alias(A, C, from, to);
        }
      }
    }
    };

    set_aliases(false);
    set_aliases(true);
  }

  void TensorCompiler::optimize_program()
  {
    //initial optimization
    //removes redundant MOV operations
    bool has_change = true;
    while (has_change)
    {
      has_change = false;
      has_change = has_change || optimize_unused_cycle();
      calculate_variable_usage_intervals();
      has_change = has_change || optimize_renaming_moves();
      has_change = has_change || optimize_self_applicable_commands();
      //printf("OPTIMIZE!!!\n");
    }
    compactify();
    optimize_copy_to_aliases();
    remove_noop();
  }

  unsigned TensorCompiler::calculate_memory_layout_naive()
  {
    unsigned total_memory = 0;
    for (auto &var : vars)
    {
      if (var.is_alias == false)
      {
        var.offset = total_memory;
        total_memory += var.total_size;
        for (auto &aid : var.aliases)
          vars[aid].offset = var.offset + vars[aid].alias_range_from;
      }
    }

    return total_memory;
  }

  void TensorCompiler::calculate_variable_usage_interval_with_aliases(unsigned v_id)
  {
    for (auto &a_id : vars[v_id].aliases)
    {
      calculate_variable_usage_interval_with_aliases(a_id);
      fu[v_id] = std::min(fu[a_id], fu[v_id]);
      fm[v_id] = std::min(fm[a_id], fm[v_id]);

      lu[v_id] = std::max(lu[a_id], lu[v_id]);
      lm[v_id] = std::max(lm[a_id], lm[v_id]);
    }
  }

  struct Chunk
  {
    unsigned t_start = 0;
    unsigned t_end = 0;
    unsigned offset = 0;
    unsigned size = 0;
    unsigned rnd = 0;
    unsigned var_id = 0;
  };

  unsigned TensorCompiler::calculate_memory_layout_interval_coloring()
  {
    /*
    Relatively simple and unefficient approach to find optimal memory layout
    It is based on Greedy Interval Coloring Algorithm
    see https://homepages.gac.edu/~sskulrat/Courses/2015F-375/lectures/g2.pdf 
    It is optimal if all chunks of memory are the same size but gurantees
    nothing othrwise (in our case)
    In most neural networks I tested this method can reduce memory requirements 
    by 15-35%
    */
    unsigned total_memory = 0;
    for (auto &var : vars)
    {
      if (var.is_alias == false)
      {
        var.offset = total_memory;
        total_memory += var.total_size;
        for (auto &aid : var.aliases)
          vars[aid].offset = var.offset + vars[aid].alias_range_from;
      }
    }
    calculate_variable_usage_intervals();
    std::vector<Chunk> chunks;
    for (int i=0;i<vars.size();i++)
    {
      if (vars[i].is_alias == false && vars[i].total_size > 0)
      {
        calculate_variable_usage_interval_with_aliases(i);
        chunks.push_back(Chunk{std::min(fu[i],fm[i]), std::max(lu[i], lm[i]) + 1u, 0u, vars[i].total_size, (unsigned)i , (unsigned)i});
        //printf("R%u [%u %u] size %u\n", (unsigned)(i), chunks.back().t_start, chunks.back().t_end, chunks.back().size);
      }
    }

    std::vector<std::map<unsigned, unsigned>> regions(commands.size()+1); //regions[i][offset] = free_size, size = 0 means last region
    for (int i=0;i<regions.size();i++)
      regions[i][0] = total_memory;

    std::sort(chunks.begin(), chunks.end(), [&](const Chunk & a, const Chunk & b) -> bool 
            {    
              if (a.t_start != b.t_start)
                return a.t_start < b.t_start;
              else if (a.t_end != b.t_end)
                return a.t_end > b.t_end;
              else
                return a.rnd < b.rnd;
            });

    for (int i=0;i<chunks.size();i++)
    {
      /*
      printf("Chunk %u(%u) [%u %u] size %u\n", (unsigned)(i), chunks[i].var_id, chunks[i].t_start, chunks[i].t_end, chunks[i].size);
      int r = 0;
      for (auto &R : regions)
      {
        printf("region %d size %d:", r, (int)R.size()); 
        for (auto &p : R)
          printf("(%u %u)", p.first, p.second);
        printf("\n");
        r++;
      }
      */

      bool found_fit_region = false; 
      std::vector<unsigned> fit_starts;
      for (auto &p : regions[chunks[i].t_start])
      {
        if (p.second < chunks[i].size)
          continue;

        //printf("region%u [%u %u] is candidate for placement\n", chunks[i].t_start, p.first, p.second);
        fit_starts = {p.first};
        bool fit_all = true;
        for (int r = chunks[i].t_start + 1; r < chunks[i].t_end; r++)
        {
          bool r_fit = false;
          // find region in this time that can fit this chunk
          for (auto &next_p : regions[r])
          {
            //printf("region%u [%u %u] is ??? for placement\n", r, next_p.first, next_p.second);
            if (next_p.first <= p.first && next_p.first + next_p.second >= p.first + chunks[i].size)
            {
              //printf("region%u [%u %u] is ok for placement\n", r, next_p.first, next_p.second);
              r_fit = true;
              fit_starts.push_back(next_p.first);
              break;
            }
          }
          if (!r_fit)
          {
            fit_all = false;
            break;
          }
        }
        
        if (!fit_all)
          continue;

        found_fit_region = true;
      }

      assert(found_fit_region);

      chunks[i].offset = fit_starts[0];

      //update regions
      for (int r = chunks[i].t_start; r < chunks[i].t_end; r++)
      {
        unsigned start = fit_starts[r-chunks[i].t_start];
        unsigned size = regions[r].at(start);
        if (chunks[i].offset > start)
          regions[r].at(start) = chunks[i].offset - start;
        else
          regions[r].erase(start);
        if (chunks[i].offset + chunks[i].size < start + size)
          regions[r][chunks[i].offset + chunks[i].size] = start + size - (chunks[i].offset + chunks[i].size);
      }
    }

    //for (int i=0;i<chunks.size();i++)
    //  printf("Chunk %u [%u %u] size %u offset %u\n", (unsigned)(i), chunks[i].t_start, chunks[i].t_end, chunks[i].size, chunks[i].offset);

    //printf("####\n");
    //for (int i=0;i<vars.size();i++)
    //  if (vars[i].is_alias == false)
    //    printf("var%d [%u %u]\n", i, vars[i].offset, vars[i].offset + vars[i].total_size);


    unsigned comp_memory = 0;
    for (auto &ch : chunks)
    {
      auto &var = vars[ch.var_id];
      //printf("%u offset %u %u\n", ch.var_id, var.total_size, var.offset);
      var.offset = ch.offset;
      comp_memory = std::max(comp_memory, var.offset + var.total_size);
      for (auto &aid : var.aliases)
        vars[aid].offset = var.offset + vars[aid].alias_range_from;
    }

    //printf("####\n");
    //for (int i=0;i<vars.size();i++)
    //  if (vars[i].is_alias == false)
    //    printf("var%d [%u %u]\n", i, vars[i].offset, vars[i].offset + vars[i].total_size);

    //printf("memory %u/%u\n",comp_memory, total_memory);

    return comp_memory;
  }

  TensorProgram TensorCompiler::finish_program(bool print_program)
  {
    optimize_program();
    unsigned total_memory_req = calculate_memory_layout_interval_coloring();

    TensorProgram pr;
    pr.commands = commands;
    pr.constants = constants;
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

    if (print_program)
    {
      printf("finished recording tensor program\n");
      printf("requires %d bytes of memory\n", (int)(sizeof(float)*total_memory_req));

      printf("%d variables\n", (int)vars.size());
      for (unsigned vid = 0; vid < vars.size(); vid++)
      {
        auto &var = vars[vid];
        printf("V%-2u:[%5u] %u %5u [%3u %3u %3u %3u %3u %3u] %s %s ", vid, var.offset, var.Dim, var.total_size, 
               var.sizes[0], var.sizes[1], var.sizes[2], var.sizes[3], var.sizes[4], var.sizes[5],
               var.is_input ? " input" : "      ", var.is_output ? "output" : "      ");
        if (var.is_alias)
        {
          printf("alias of %2u [%2u %2u]", var.alias_master_id, var.alias_range_from, var.alias_range_to);
        }
        else if (var.aliases.size() > 0)
        {
          printf("master of { ");
          for (auto &aid : var.aliases)
            printf("%2u ", aid);
          printf("}");
        }


        printf("\n");
      }

      printf("%d commands\n", (int)commands.size());
      unsigned cid = 0;
      for (auto &cmd : commands)
      {
        const char *cmd_name = TensorProgram::cmd_properties[cmd.type].name.c_str();
        printf("Cmd %2u: %-8s %2u %2u %2u - %3u %3u %3u %3u %3u\n", cid, cmd_name, cmd.args[0],cmd.args[1],cmd.args[2],cmd.args[3],
                                                                    cmd.args[4],cmd.args[5],cmd.args[6],cmd.args[7]);
        cid++;
      }
    }
    return pr;
  }

  void TensorCompiler::add_command(TensorProgram::CommandType type, unsigned A, unsigned B, unsigned C, 
                                   unsigned arg0, unsigned arg1, unsigned arg2, unsigned arg3, unsigned arg4)
  {
    commands.push_back({type, {A, B, C, arg0, arg1, arg2, arg3, arg4}});
    //printf("%s add command %u %u %u %u %u %u (%d %d) -> (%d %d)\n", TensorProgram::cmd_properties[type].name.c_str(),
    //A, B, C, arg0, arg1, arg2, vars[A].sizes[0], vars[A].sizes[1], vars[C].sizes[0], vars[C].sizes[1]);
  }
}