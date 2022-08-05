#include "weber_penn_generator.h"
#include "common_utils/parameter_save_load_defines.h"

void WeberPennParametersNative::save_load_define(SaveLoadMode mode, Block &b, ParameterList &list)
{
    P_INT(shape, mode);
    P_FLOAT(g_scale, mode);
    P_FLOAT(g_scale_v, mode);
    P_INT(levels, mode);
    P_FLOAT(ratio, mode);
    P_FLOAT(ratio_power, mode);
    P_FLOAT(flare, mode);
    P_INT(base_splits, mode);

    P_VEC4(base_size, mode);
    P_VEC4(down_angle, mode);
    P_VEC4(down_angle_v, mode);
    P_VEC4(rotate, mode);
    P_VEC4(rotate_v, mode);
    P_VEC4(branches, mode);
    P_VEC4(length, mode);
    P_VEC4(length_v, mode);
    P_VEC4(taper, mode);
    P_VEC4(seg_splits, mode);
    P_VEC4(split_angle, mode);
    P_VEC4(split_angle_v, mode);
    P_VEC4(curve_res, mode);
    P_VEC4(curve, mode);
    P_VEC4(curve_back, mode);
    P_VEC4(curve_v, mode);
    P_VEC4(bend_v, mode); 
    P_VEC4(branch_dist, mode);
    P_VEC4(radius_mod, mode);

    P_INT(leaf_blos_num, mode);
    P_INT(leaf_shape, mode);
    P_FLOAT(leaf_scale, mode);
    P_FLOAT(leaf_scale_x, mode);
    P_FLOAT(leaf_bend, mode);
    P_INT(blossom_shape, mode);
    P_FLOAT(blossom_scale, mode);
    P_FLOAT(blossom_rate, mode);
    P_FLOAT(leaf_rate, mode);
    P_VEC3(tropism, mode);
    P_FLOAT(prune_ratio, mode);
    P_FLOAT(prune_width, mode);
    P_FLOAT(prune_width_peak, mode);
    P_FLOAT(prune_power_high, mode); 
    P_FLOAT(prune_power_low, mode);
}