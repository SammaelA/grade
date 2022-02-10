#include "GE_generator.h"

void GETreeParameters::save_to_blk(Block &b)
{
    b.set_double("lambda", lambda);
    b.set_double("k", k);
    b.set_int("tau", tau);
    b.set_double("ro", ro);
    b.set_double("X0", X0);
    b.set_double("Xm", Xm);
    b.set_double("r", r);
    b.set_int("alpha", alpha);
    b.set_double("sigma", sigma);
    b.set_double("mu", mu);
    b.set_double("b_min", b_min);
    b.set_double("b_max", b_max);
    b.set_double("r_s", r_s);
    b.set_double("rs_size_factor", rs_size_factor);
    b.set_int("remove_min_level", remove_min_level);
    b.set_double("base_r", base_r);
    b.set_int("max_branches", max_branches);
    b.set_int("occlusion_pyramid_d", occlusion_pyramid_d);
    b.set_double("r_pow", r_pow);
    b.set_int("sp_points_base", sp_points_base);
    b.set_double("branching_angle_min", branching_angle_min);
    b.set_double("branching_angle_max", branching_angle_max);
    b.set_int("max_iterations", max_iterations);
    b.set_double("leaf_size_mult", leaf_size_mult);
    b.set_double("leaves_cnt", leaves_cnt);
    b.set_int("max_joints_in_branch", max_joints_in_branch);
    b.set_double("resource_mult", resource_mult);
    b.set_double("res_q", res_q);
    b.set_double("leaves_max_r", leaves_max_r);
    b.set_double("leaves_angle_a", leaves_angle_a);
    b.set_double("leaves_angle_b", leaves_angle_b);
    b.set_int("root_type", root_type);
    b.set_vec2("tropism_min_max", tropism_min_max);
    b.set_vec4("tropism_params", tropism_params);
    b.set_double("branching_tropims_mult", branching_tropims_mult);
    b.set_int("tropism_level_base", tropism_level_base);
    b.set_double("res_decrease_step", res_decrease_step);
    b.set_double("res_decrease_min", res_decrease_min);
    b.set_vec2("initial_trunk_scale", initial_trunk_scale);
    b.set_double("trunk_bonus_radius", trunk_bonus_radius);
    b.set_double("trunk_bonus_radius_mod", trunk_bonus_radius_mod);
    b.set_vec2("distance_fine_start", distance_fine_start);
    b.set_vec2("distance_fine_slope", distance_fine_slope);

    b.set_double("L0",L0);
    b.set_double("R0",R0);
    b.set_int("B_cnt",B_cnt);
    b.set_double("dX",dX);
    b.set_vec4("iR",iR);
    b.set_vec4("dR",dR);
    b.set_vec4("Lt",Lt);
    b.set_vec4("Lb",Lb);
    b.set_vec4("phi",phi);
    b.set_vec4("psi",psi);
}

void GETreeParameters::load_from_blk(Block &b)
{
    lambda = b.get_double("lambda", lambda);
    k = b.get_double("k", k);
    tau = b.get_int("tau", tau);
    ro = b.get_double("ro", ro);
    X0 = b.get_double("X0", X0);
    Xm = b.get_double("Xm", Xm);
    r = b.get_double("r", r);
    alpha = b.get_int("alpha", alpha);
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
    initial_trunk_scale = b.get_vec2("initial_trunk_scale", initial_trunk_scale);
    trunk_bonus_radius = b.get_double("trunk_bonus_radius", trunk_bonus_radius);
    trunk_bonus_radius_mod = b.get_double("trunk_bonus_radius_mod", trunk_bonus_radius_mod);
    distance_fine_start = b.get_vec2("distance_fine_start", distance_fine_start);
    distance_fine_slope = b.get_vec2("distance_fine_slope", distance_fine_slope);

    L0 = b.get_double("L0",L0);
    R0 = b.get_double("R0",R0);
    B_cnt = b.get_int("B_cnt",B_cnt);
    dX = b.get_double("dX",dX);
    iR = b.get_vec4("iR",iR);
    dR = b.get_vec4("dR",dR);
    Lt = b.get_vec4("Lt",Lt);
    Lb = b.get_vec4("Lb",Lb);
    phi = b.get_vec4("phi",phi);
    psi = b.get_vec4("psi",psi);
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
    ORD1(B_cnt);

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

    CON2(initial_trunk_scale_x, initial_trunk_scale.x);
    CON2(initial_trunk_scale_y, initial_trunk_scale.y);
    CON1(trunk_bonus_radius);
    CON1(trunk_bonus_radius_mod);

    CON2(distance_fine_start_x, distance_fine_start.x);
    CON2(distance_fine_start_y, distance_fine_start.y);
    CON2(distance_fine_slope_x, distance_fine_slope.x);
    CON2(distance_fine_slope_y, distance_fine_slope.y);

    CON1(L0);
    CON1(R0);
    CON1(dX);


    CON2(iR_x, iR.x);
    CON2(iR_y, iR.y);
    CON2(iR_z, iR.z);
    CON2(iR_w, iR.w);

    CON2(dR_x, dR.x);
    CON2(dR_y, dR.y);
    CON2(dR_z, dR.z);
    CON2(dR_w, dR.w);

    CON2(Lt_x, Lt.x);
    CON2(Lt_y, Lt.y);
    CON2(Lt_z, Lt.z);
    CON2(Lt_w, Lt.w);

    CON2(Lb_x, Lb.x);
    CON2(Lb_y, Lb.y);
    CON2(Lb_z, Lb.z);
    CON2(Lb_w, Lb.w);

    CON2(phi_x, phi.x);
    CON2(phi_y, phi.y);
    CON2(phi_z, phi.z);
    CON2(phi_w, phi.w);

    CON2(psi_x, psi.x);
    CON2(psi_y, psi.y);
    CON2(psi_z, psi.z);
    CON2(psi_w, psi.w);
}
