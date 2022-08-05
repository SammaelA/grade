#include "GE_generator.h"
#include "common_utils/parameter_save_load_defines.h"

void GETreeParameters::save_load_define(SaveLoadMode mode, Block &b, ParameterList &list)
{
    P_CAT(root_type, mode);

    P_INT(tau, mode);
    P_INT(alpha, mode);
    P_INT(remove_min_level, mode);
    P_INT(max_branches, mode);
    P_INT(occlusion_pyramid_d, mode);
    P_INT(sp_points_base, mode);
    P_INT(max_iterations, mode);
    P_INT(max_joints_in_branch, mode);
    P_INT(tropism_level_base, mode);
    P_INT(B_cnt, mode);

    P_FLOAT(lambda, mode);
    P_FLOAT(k, mode);
    P_FLOAT(ro, mode);
    P_FLOAT(X0, mode);
    P_FLOAT(Xm, mode);
    P_FLOAT(r, mode);
    P_FLOAT(sigma, mode);
    P_FLOAT(mu, mode);
    P_FLOAT(b_min, mode);
    P_FLOAT(b_max, mode);
    P_FLOAT(r_s, mode);
    P_FLOAT(rs_size_factor, mode);

    P_FLOAT(base_r, mode);
    P_FLOAT(r_pow, mode);
    P_FLOAT(branching_angle_min, mode);
    P_FLOAT(branching_angle_max, mode);
    P_FLOAT(leaf_size_mult, mode);
    P_FLOAT(leaves_cnt, mode);
    P_FLOAT(resource_mult, mode);
    P_FLOAT(res_q, mode);
    P_FLOAT(leaves_max_r, mode);
    P_FLOAT(leaves_angle_a, mode);
    P_FLOAT(leaves_angle_b, mode);
    P_VEC4(tropism_params, mode);
    P_VEC2(tropism_min_max, mode);
    P_FLOAT(branching_tropims_mult, mode);
    P_FLOAT(res_decrease_step, mode);
    P_FLOAT(res_decrease_min, mode);
    P_VEC2(initial_trunk_scale, mode);
    P_FLOAT(trunk_bonus_radius, mode);
    P_FLOAT(trunk_bonus_radius_mod, mode);
    P_VEC2(distance_fine_start, mode);
    P_VEC2(distance_fine_slope, mode);
    P_FLOAT(L0, mode);
    P_FLOAT(R0, mode);
    P_FLOAT(dX, mode);

    P_VEC4(iR, mode);
    P_VEC4(dR, mode);
    P_VEC4(Lt, mode);
    P_VEC4(Lb, mode);
    P_VEC4(phi, mode);
    P_VEC4(psi, mode);
}
