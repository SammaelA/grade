#include "weber_penn_parameters.h"

void WeberPennParameters::get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                                             ParameterVariablesSet v_set) 
{
    list = {
        {{ParameterMaskValues::CONSTANT,"shape"}, shape},
        {{ParameterMaskValues::CONSTANT,"g_scale"}, g_scale},
        {{ParameterMaskValues::CONSTANT,"g_scale_v"}, g_scale_v},
        {{ParameterMaskValues::CONSTANT,"levels"}, levels},
        {{ParameterMaskValues::CONSTANT,"ratio"}, ratio},
        {{ParameterMaskValues::CONSTANT,"ratio_power"}, ratio_power},
        {{ParameterMaskValues::CONSTANT,"flare"}, flare},
        {{ParameterMaskValues::CONSTANT,"base_splits"}, base_splits},
        {{ParameterMaskValues::CONSTANT,"base_size"}, base_size},
        {{ParameterMaskValues::CONSTANT,"down_angle"}, down_angle},
        {{ParameterMaskValues::CONSTANT,"down_angle_v"}, down_angle_v},
        {{ParameterMaskValues::CONSTANT,"rotate"}, rotate},
        {{ParameterMaskValues::CONSTANT,"rotate_v"}, rotate_v},
        {{ParameterMaskValues::CONSTANT,"branches"}, branches},
        {{ParameterMaskValues::CONSTANT,"length"}, length},
        {{ParameterMaskValues::CONSTANT,"length_v"}, length_v},
        {{ParameterMaskValues::CONSTANT,"taper"}, taper},
        {{ParameterMaskValues::CONSTANT,"seg_splits"}, seg_splits},
        {{ParameterMaskValues::CONSTANT,"split_angle"}, split_angle},
        {{ParameterMaskValues::CONSTANT,"split_angle_v"}, split_angle_v},
        {{ParameterMaskValues::CONSTANT,"curve_res"}, curve_res},
        {{ParameterMaskValues::CONSTANT,"curve"}, curve},
        {{ParameterMaskValues::CONSTANT,"curve_back"}, curve_back},
        {{ParameterMaskValues::CONSTANT,"curve_v"}, curve_v},
        {{ParameterMaskValues::CONSTANT,"bend_v"}, bend_v},
        {{ParameterMaskValues::CONSTANT,"branch_dist"}, branch_dist},
        {{ParameterMaskValues::CONSTANT,"leaf_blos_num"}, leaf_blos_num},
        {{ParameterMaskValues::CONSTANT,"leaf_shape"}, leaf_shape},
        {{ParameterMaskValues::CONSTANT,"leaf_scale"}, leaf_scale},
        {{ParameterMaskValues::CONSTANT,"leaf_bend"}, leaf_bend}
    };
}

std::string WeberPennParameters::convert_to_python_list()
{
    #define IS_INT(a) (abs((a) - round(a))<1e-6)

    std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> list;
    get_parameter_list(list);

    std::string pylist;
    pylist = pylist + "\"\"\" "+name+" \"\"\"\n";
    pylist = pylist + "params = {\n";
    for (auto &p : list)
    {
        pylist = pylist + "\'"+p.first.name+"\': ";
        if (p.second.state_qs.empty())
        {
            if IS_INT(p.second.baseValue)
                pylist = pylist + std::to_string((int)p.second.baseValue);
            else 
                pylist = pylist + std::to_string(p.second.baseValue);
        }
        else 
        {
            pylist = pylist + "[";
            for (int i = 0;i<p.second.state_qs.size();i++)
            {
                float val = p.second.baseValue*p.second.state_qs[i];
                if IS_INT(val)
                    pylist = pylist + std::to_string((int)val);
                else 
                    pylist = pylist + std::to_string(val);
                if (i + 1 != p.second.state_qs.size())
                    pylist = pylist + ", ";
            }
            pylist = pylist + "]";
        }
        pylist = pylist + ",\n";
    }
    pylist = pylist + "}\n";

    return pylist;
}