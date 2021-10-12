#pragma once
#include "parameter.h"


struct WeberPennParameters : public ParametersSet
{
    Parameter<float> shape, g_scale, g_scale_v, levels, ratio,
                     ratio_power, flare, base_splits, base_size, down_angle,
                     down_angle_v, rotate, rotate_v, branches, length, 
                     length_v, taper, seg_splits, split_angle, split_angle_v,
                     curve_res, curve, curve_back, curve_v, bend_v,
                     branch_dist, leaf_blos_num, leaf_shape, leaf_scale, leaf_bend;
    std::string name;
    bool settings_already_in_file = false;
    virtual void set_state(int state) override
    {
        shape.set_state(state);
        g_scale.set_state(state);
        g_scale_v.set_state(state);
        levels.set_state(state);
        ratio.set_state(state);
        ratio_power.set_state(state);
        flare.set_state(state);
        base_splits.set_state(state);
        base_size.set_state(state);
        down_angle.set_state(state);
        down_angle_v.set_state(state);
        rotate.set_state(state);
        rotate_v.set_state(state);
        branches.set_state(state);
        length.set_state(state);
        length_v.set_state(state);
        taper.set_state(state);
        seg_splits.set_state(state);
        split_angle.set_state(state);
        split_angle_v.set_state(state);
        curve_res.set_state(state);
        curve.set_state(state);
        curve_back.set_state(state);
        curve_v.set_state(state);
        bend_v.set_state(state);
        branch_dist.set_state(state);
        leaf_blos_num.set_state(state);
        leaf_shape.set_state(state);
        leaf_scale.set_state(state);
        leaf_bend.set_state(state);
    }
    WeberPennParameters():
    shape(4),
    g_scale(10),
    g_scale_v(1),
    levels(3),
    ratio(0.025),
    ratio_power(1.5),
    flare(0.6),
    base_splits(-2),
    base_size(std::vector<float>{0.1, 0.4, 0.02, 0.02}),
    down_angle(std::vector<float>{-0, 50, 50, 45} ),
    down_angle_v(std::vector<float>{-0, 5, 5, 10} ),
    rotate(std::vector<float>{-0, 140, 140, 77} ),
    rotate_v(std::vector<float>{-0, 0, 0, 0} ),
    branches(std::vector<float>{1, 6, 20, 5} ),
    length(std::vector<float>{1, 0.7, 0.3, 0} ),
    length_v(std::vector<float>{0, 0.05, 0.05, 0} ),
    taper(std::vector<float>{1, 1, 1, 1} ),
    seg_splits(std::vector<float>{1.5, 1.5, 0, 0} ),
    split_angle(std::vector<float>{50, 50, 0, 0} ),
    split_angle_v(std::vector<float>{5, 5, 0, 0} ),
    curve_res(std::vector<float>{6, 5, 3, 0} ),
    curve(std::vector<float>{0, 0, 0, 0} ),
    curve_back(std::vector<float>{0, 0, 0, 0} ),
    curve_v(std::vector<float>{200, 100, 100, 0} ),
    bend_v(std::vector<float>{-0, 50, 0, 0} ),
    branch_dist(std::vector<float>{-0, 0, 0, 0} ),
    leaf_blos_num(30),
    leaf_shape(5),
    leaf_scale(0.2),
    leaf_bend(0.3)

    {

    }
    virtual void get_parameter_list(std::vector<std::pair<ParameterTinyDesc,Parameter<float> &>> &list,
                                    ParameterVariablesSet v_set = ParameterVariablesSet::ALL_VALUES) override;
    std::string convert_to_python_list();
};