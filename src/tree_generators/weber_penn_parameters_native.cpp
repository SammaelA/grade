#include "weber_penn_generator.h"

void WeberPennParametersNative::save_to_blk(Block &b)
{
    b.set_int("shape", shape);
    b.set_double("g_scale", g_scale);
    b.set_double("g_scale_v", g_scale_v);
    b.set_int("levels", levels);
    b.set_double("ratio", ratio);
    b.set_double("ratio_power", ratio_power);
    b.set_double("flare", flare);
    b.set_int("base_splits", base_splits);
    b.set_vec4("base_size", base_size);
    b.set_vec4("down_angle", down_angle);
    b.set_vec4("down_angle_v", down_angle_v);
    b.set_vec4("rotate", rotate);
    b.set_vec4("rotate_v", rotate_v);
    b.set_vec4("branches", branches);
    b.set_vec4("length", length);
    b.set_vec4("length_v", length_v);
    b.set_vec4("taper", taper);
    b.set_vec4("seg_splits", seg_splits);
    b.set_vec4("split_angle", split_angle);
    b.set_vec4("split_angle_v", split_angle_v);
    b.set_vec4("curve_res", curve_res);
    b.set_vec4("curve", curve);
    b.set_vec4("curve_back", curve_back);
    b.set_vec4("curve_v", curve_v);
    b.set_vec4("bend_v", bend_v);
    b.set_vec4("branch_dist", branch_dist);
    b.set_vec4("radius_mod", radius_mod);
    b.set_int("leaf_blos_num", leaf_blos_num);
    b.set_int("leaf_shape", leaf_shape);
    b.set_double("leaf_scale", leaf_scale);
    b.set_double("leaf_scale_x", leaf_scale_x);
    b.set_double("leaf_bend", leaf_bend);
    b.set_int("blossom_shape", blossom_shape);
    b.set_double("blossom_scale", blossom_scale);
    b.set_double("blossom_rate", blossom_rate);
    b.set_vec3("tropism", tropism);
    b.set_double("prune_ratio", prune_ratio);
    b.set_double("prune_width", prune_width);
    b.set_double("prune_width_peak", prune_width_peak);
    b.set_double("prune_power_high", prune_power_high);
    b.set_double("prune_power_low", prune_power_low);
}

void WeberPennParametersNative::load_from_blk(Block &b)
{
    shape = b.get_int("shape", shape);
    g_scale = b.get_double("g_scale", g_scale);
    g_scale_v = b.get_double("g_scale_v", g_scale_v);
    levels = b.get_int("levels", levels);
    ratio = b.get_double("ratio", ratio);
    ratio_power = b.get_double("ratio_power", ratio_power);
    flare = b.get_double("flare", flare);
    base_splits = b.get_int("base_splits", base_splits);
    base_size = b.get_vec4("base_size", base_size);
    down_angle = b.get_vec4("down_angle", down_angle);
    down_angle_v = b.get_vec4("down_angle_v", down_angle_v);
    rotate = b.get_vec4("rotate", rotate);
    rotate_v = b.get_vec4("rotate_v", rotate_v);
    branches = b.get_vec4("branches", branches);
    length = b.get_vec4("length", length);
    length_v = b.get_vec4("length_v", length_v);
    taper = b.get_vec4("taper", taper);
    seg_splits = b.get_vec4("seg_splits", seg_splits);
    split_angle = b.get_vec4("split_angle", split_angle);
    split_angle_v = b.get_vec4("split_angle_v", split_angle_v);
    curve_res = b.get_vec4("curve_res", curve_res);
    curve = b.get_vec4("curve", curve);
    curve_back = b.get_vec4("curve_back", curve_back);
    curve_v = b.get_vec4("curve_v", curve_v);
    bend_v = b.get_vec4("bend_v", bend_v);
    branch_dist = b.get_vec4("branch_dist", branch_dist);
    radius_mod = b.get_vec4("radius_mod", radius_mod);
    leaf_blos_num = b.get_int("leaf_blos_num", leaf_blos_num);
    leaf_shape = b.get_int("leaf_shape", leaf_shape);
    leaf_scale = b.get_double("leaf_scale", leaf_scale);
    leaf_scale_x = b.get_double("leaf_scale_x", leaf_scale_x);
    leaf_bend = b.get_double("leaf_bend", leaf_bend);
    blossom_shape = b.get_int("blossom_shape", blossom_shape);
    blossom_scale = b.get_double("blossom_scale", blossom_scale);
    blossom_rate = b.get_double("blossom_rate", blossom_rate);
    tropism = b.get_vec3("tropism", tropism);
    prune_ratio = b.get_double("prune_ratio", prune_ratio);
    prune_width = b.get_double("prune_width", prune_width);
    prune_width_peak = b.get_double("prune_width_peak", prune_width_peak);
    prune_power_high = b.get_double("prune_power_high", prune_power_high);
    prune_power_low = b.get_double("prune_power_low", prune_power_low);
}

void WeberPennParametersNative::RW_parameter_list(bool write, ParameterList &list)
{
    #define ORD1(par) if (write) {list.ordinalParameters.emplace(#par, par);} else {par = list.ordinalParameters.at(#par);}
    #define ORD2(par_name, par) if (write) {list.ordinalParameters.emplace(#par_name, par);} else {par = list.ordinalParameters.at(#par_name);}

    #define CAT1(par) if (write) {list.categorialParameters.emplace(#par, par);} else {par = list.categorialParameters.at(#par);}
    #define CAT2(par_name, par) if (write) {list.categorialParameters.emplace(#par_name, par);} else {par = list.categorialParameters.at(#par_name);}

    #define CON1(par) if (write) {list.continuousParameters.emplace(#par, par);} else {par = list.continuousParameters.at(#par);}
    #define CON2(par_name, par) if (write) {list.continuousParameters.emplace(#par_name, par);} else {par = list.continuousParameters.at(#par_name);}
    #define CON_V4(par) if (write) {\
    list.continuousParameters.emplace(std::string(#par)+"_x", par.x);\
    list.continuousParameters.emplace(std::string(#par)+"_y", par.y);\
    list.continuousParameters.emplace(std::string(#par)+"_z", par.z);\
    list.continuousParameters.emplace(std::string(#par)+"_w", par.w);\
    } else {\
    par.x = list.continuousParameters.at(std::string(#par)+"_x");\
    par.y = list.continuousParameters.at(std::string(#par)+"_y");\
    par.z = list.continuousParameters.at(std::string(#par)+"_z");\
    par.w = list.continuousParameters.at(std::string(#par)+"_w");}

    ORD1(shape);
    CON1(g_scale);
    CON1(g_scale_v);
    ORD1(levels);
    CON1(ratio);
    CON1(ratio_power);
    CON1(flare);
    ORD1(base_splits);

    CON_V4(base_size);
    CON_V4(down_angle);
    CON_V4(down_angle_v);
    CON_V4(rotate);
    CON_V4(rotate_v);
    CON_V4(branches);
    CON_V4(length);
    CON_V4(length_v);
    CON_V4(taper);
    CON_V4(seg_splits);
    CON_V4(split_angle);
    CON_V4(split_angle_v);
    CON_V4(curve_res);
    CON_V4(curve);
    CON_V4(curve_back);
    CON_V4(curve_v);
    CON_V4(bend_v); 
    CON_V4(branch_dist);
    CON_V4(radius_mod);

    ORD1(leaf_blos_num);
    ORD1(leaf_shape);
    CON1(leaf_scale);
    CON1(leaf_scale_x);
    CON1(leaf_bend);
    ORD1(blossom_shape);
    CON1(blossom_scale);
    CON1(blossom_rate);
    CON2(tropism_x, tropism.x);
    CON2(tropism_y, tropism.y);
    CON2(tropism_z, tropism.z);
    CON1(prune_ratio);
    CON1(prune_width);
    CON1(prune_width_peak);
    CON1(prune_power_high); 
    CON1(prune_power_low);
}