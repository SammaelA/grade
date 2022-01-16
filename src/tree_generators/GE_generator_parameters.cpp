#include "GE_generator.h"

void GETreeParameters::save_to_blk(Block &b)
{
    b.set_double("lambda", 0.52);
    b.set_double("k", 0.75);
    b.set_int("tau", 6);
    b.set_double("ro", 1.0);
    b.set_double("X0", 2);
    b.set_double("Xm", 100);
    b.set_double("r", 0.37);
    b.set_int("alpha", 4);
    b.set_double("sigma", 0.5);
    b.set_double("mu", 1.5);
    b.set_double("b_min", 1.8);
    b.set_double("b_max", 2.2);
    b.set_double("r_s", 0.05);
    b.set_double("rs_size_factor", 0.02);
    b.set_int("remove_min_level", 2);

    b.set_double("base_r", 0.025);
    b.set_int("max_branches", 1);
    b.set_int("occlusion_pyramid_d", 10);
    b.set_double("r_pow", 2.2);
    b.set_int("sp_points_base", 16);
    b.set_double("branching_angle_min", 0);
    b.set_double("branching_angle_max", PI / 3);
    b.set_int("max_iterations", 100);
    b.set_double("leaf_size_mult", 3.5);
    b.set_double("leaves_cnt", 1.0);
    b.set_int("max_joints_in_branch", 16);
    b.set_double("resource_mult", 7.5);
    b.set_double("res_q", 1.0);
    b.set_double("leaves_max_r", 5);
    b.set_double("leaves_angle_a", 0.3);
    b.set_double("leaves_angle_b", 0.4);
    b.set_int("root_type", 1);
    b.set_vec2("tropism_min_max", tropism_min_max);
    b.set_vec4("tropism_params", tropism_params);
}

void GETreeParameters::load_from_blk(Block &b)
{
    lambda = b.get_double("lambda", 0.52);
    k = b.get_double("k", 0.75);
    tau = b.get_int("tau", 6);
    ro = b.get_double("ro", 1.0);
    X0 = b.get_double("X0", 2);
    Xm = b.get_double("Xm", 100);
    r = b.get_double("r", 0.37);
    alpha = b.get_int("alpha", 4);
    sigma = b.get_double("sigma", sigma);
    mu = b.get_double("mu", mu);
    b_min = b.get_double("b_min", b_min);
    b_max = b.get_double("b_max", b_max);
    r_s = b.get_double("r_s", r_s);
    rs_size_factor = b.get_double("rs_size_factor", rs_size_factor);
    remove_min_level = b.get_int("remove_min_level", remove_min_level);

    base_r = b.get_double("base_r", base_r);
    max_branches = b.get_int("max_branches", max_branches);
    occlusion_pyramid_d = b.get_int("occlusion_pyramid_d", occlusion_pyramid_d);
    r_pow = b.get_double("r_pow", r_pow);
    sp_points_base = b.get_int("sp_points_base", sp_points_base);
    branching_angle_min = b.get_double("branching_angle_min", branching_angle_min);
    branching_angle_max = b.get_double("branching_angle_max", branching_angle_max);
    max_iterations = b.get_int("max_iterations", max_iterations);
    leaf_size_mult = b.get_double("leaf_size_mult", leaf_size_mult);
    leaves_cnt = b.get_double("leaves_cnt", leaves_cnt);
    max_joints_in_branch = b.get_int("max_joints_in_branch", max_joints_in_branch);
    resource_mult = b.get_double("resource_mult", resource_mult);
    res_q = b.get_double("res_q", res_q);
    leaves_max_r = b.get_double("leaves_max_r", leaves_max_r);
    leaves_angle_a = b.get_double("leaves_angle_a", leaves_angle_a);
    leaves_angle_b = b.get_double("leaves_angle_b", leaves_angle_b);
    root_type = b.get_int("root_type", root_type);
    tropism_min_max = b.get_vec2("tropism_min_max", tropism_min_max);
    tropism_params = b.get_vec4("tropism_params", tropism_params);
    branching_tropims_mult = b.get_double("branching_tropims_mult", branching_tropims_mult);
    tropism_level_base = b.get_int("tropism_level_base", tropism_level_base);
    res_decrease_step = b.get_double("res_decrease_step", res_decrease_step);
    res_decrease_min = b.get_double("res_decrease_min", res_decrease_min);
}
void GETreeParameters::RW_parameter_list(bool write, ParameterList &list)
{
    #define ORD1(par) if (write) {list.ordinalParameters.emplace(#par, par);} else {par = list.ordinalParameters.at(#par);}
    #define ORD2(par_name, par) if (write) {list.ordinalParameters.emplace(#par_name, par);} else {par = list.ordinalParameters.at(#par_name);}

    #define CAT1(par) if (write) {list.categorialParameters.emplace(#par, par);} else {par = list.categorialParameters.at(#par);}
    #define CAT2(par_name, par) if (write) {list.categorialParameters.emplace(#par_name, par);} else {par = list.categorialParameters.at(#par_name);}

    #define CON1(par) if (write) {list.continuousParameters.emplace(#par, par);} else {par = list.continuousParameters.at(#par);}
    #define CON2(par_name, par) if (write) {list.continuousParameters.emplace(#par_name, par);} else {par = list.continuousParameters.at(#par_name);}

    CAT1(root_type);

    ORD1(tau);
    ORD1(alpha);
    ORD1(remove_min_level);
    ORD1(max_branches);
    ORD1(occlusion_pyramid_d);
    ORD1(sp_points_base);
    ORD1(max_iterations);
    ORD1(max_joints_in_branch);
    ORD1(tropism_level_base);

    CON1(lambda);
    CON1(k);
    CON1(ro);
    CON1(X0);
    CON1(Xm);
    CON1(r);
    CON1(sigma);
    CON1(mu);
    CON1(b_min);
    CON1(b_max);
    CON1(r_s);
    CON1(rs_size_factor);

    CON1(base_r);
    CON1(r_pow);
    CON1(branching_angle_min);
    CON1(branching_angle_max);
    CON1(leaf_size_mult);
    CON1(leaves_cnt);
    CON1(resource_mult);
    CON1(res_q);
    CON1(leaves_max_r);
    CON1(leaves_angle_a);
    CON1(leaves_angle_b);

    CON2(tropism_params_x, tropism_params.x);
    CON2(tropism_params_y, tropism_params.y);
    CON2(tropism_params_z, tropism_params.z);
    CON2(tropism_params_w, tropism_params.w);

    CON2(tropism_min_max_x, tropism_min_max.x);
    CON2(tropism_min_max_y, tropism_min_max.y);

    CON1(branching_tropims_mult);
    CON1(res_decrease_step);
    CON1(res_decrease_min);
}
