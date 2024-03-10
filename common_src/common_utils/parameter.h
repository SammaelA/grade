#pragma once

#include "common_utils/utility.h"
#include "common_utils/blk.h"
#include <vector>

struct CategorialParameter
{
    int val = 0;
    std::vector<int> possible_values;
    operator int () {return val;}
    CategorialParameter(int _val) 
    {
        val = _val;
        possible_values = {_val};
    }
    bool fixed() { return possible_values.size() == 1;}
};

struct OrdinalParameter
{
    int val = 0;
    int min_val = 0;
    int max_val = 0;
    operator int () {return val;}
    OrdinalParameter(int _val) 
    {
        val = _val;
        min_val = val;
        max_val = val;
    }
    bool fixed() { return max_val <= min_val;}
};

struct ContinuousParameter
{
    float val = 0;
    float min_val = 0;
    float max_val = 0;
    operator float () {return val;}
    ContinuousParameter(int _val) 
    {
        val = _val;
        min_val = val;
        max_val = val;
    }
    ContinuousParameter(float _val) 
    {
        val = _val;
        min_val = val;
        max_val = val;
    }
    bool fixed() { return max_val <= min_val;}
};

struct ParameterList
{
    std::map<std::string, CategorialParameter> categorialParameters;
    std::map<std::string, OrdinalParameter> ordinalParameters;
    std::map<std::string, ContinuousParameter> continuousParameters;

    void print();
    void load_borders_from_blk(Block &b);
    void to_simple_list(std::vector<float> &list, bool normalized = false, bool remove_fixed_params = false);
    void from_simple_list(std::vector<float> &list, bool normalized = false, bool remove_fixed_params = false);
    float diff(ParameterList &list, bool normalized = false, bool remove_fixed_params = false);
};

struct ParameterSet
{
    enum SaveLoadMode
    {
      BLK_SAVE,
      BLK_LOAD,
      PAR_LIST_SAVE,
      PAR_LIST_LOAD
    };
    virtual ~ParameterSet() = default;
    virtual float3 get_tree_max_size() = 0;
    virtual ParameterSet *copy() { return nullptr;};
    virtual float get_scale_factor() {return 1;}
    virtual void save_to_blk(Block &b)
    {
      ParameterList list;
      save_load_define(SaveLoadMode::BLK_SAVE, b, list);
    }
    virtual void load_from_blk(Block &b)
    {
      ParameterList list;
      save_load_define(SaveLoadMode::BLK_LOAD, b, list);
    }
    virtual void write_parameter_list(ParameterList &list)
    { 
      Block b;
      save_load_define(SaveLoadMode::PAR_LIST_LOAD, b, list);
    }
    virtual void read_parameter_list(ParameterList &list)
    {
      Block b;
      save_load_define(SaveLoadMode::PAR_LIST_SAVE, b, list);
    }
private:
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) {};
};