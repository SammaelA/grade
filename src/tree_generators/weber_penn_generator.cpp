#include "weber_penn_generator.h"
#include "common_utils/turtle.h"
#include <atomic>
#include <cmath>

using glm::vec3;
using glm::vec4;
using glm::mat4;
using glm::degrees;

float declination(vec3 vec)
{
    return degrees(atan2(sqrt(SQR(vec.x) + SQR(vec.y)), vec.z));
}

void vec_to_quat(glm::quat &q, const vec3 vec, short axis, const short upflag)
{
  //it is taken from blender/python/mathutils/mathutils_Vector.c with mathutils vecs
  //and quats replaced with glm ones
  const float eps = 1e-4f;
  vec3 nor, tvec;
  float angle, si, co, len;

  /* first set the quat to unit */
  q = glm::quat(1,0,0,0);

  len = length(vec);

  if (len == 0.0f) {
    return;
  }

  /* rotate to axis */
  if (axis > 2) {
    tvec = vec;
    axis = (short)(axis - 3);
  }
  else {
    tvec = -vec;
  }

  /* nasty! I need a good routine for this...
   * problem is a rotation of an Y axis to the negative Y-axis for example.
   */

  if (axis == 0) { /* x-axis */
    nor[0] = 0.0;
    nor[1] = -tvec[2];
    nor[2] = tvec[1];

    if (fabsf(tvec[1]) + fabsf(tvec[2]) < eps) {
      nor[1] = 1.0f;
    }

    co = tvec[0];
  }
  else if (axis == 1) { /* y-axis */
    nor[0] = tvec[2];
    nor[1] = 0.0;
    nor[2] = -tvec[0];

    if (fabsf(tvec[0]) + fabsf(tvec[2]) < eps) {
      nor[2] = 1.0f;
    }

    co = tvec[1];
  }
  else { /* z-axis */
    nor[0] = -tvec[1];
    nor[1] = tvec[0];
    nor[2] = 0.0;

    if (fabsf(tvec[0]) + fabsf(tvec[1]) < eps) {
      nor[0] = 1.0f;
    }

    co = tvec[2];
  }
  co /= len;

  nor = normalize(nor);

  q = glm::angleAxis(acosf(co), nor);
  
  if (axis != upflag) {
    glm::mat3x3 mat;
    glm::quat q2;
    const vec3 &fp = mat[2];
    mat = glm::mat3x3(q);

    if (axis == 0) {
      if (upflag == 1) {
        angle = 0.5f * atan2f(fp[2], fp[1]);
      }
      else {
        angle = -0.5f * atan2f(fp[1], fp[2]);
      }
    }
    else if (axis == 1) {
      if (upflag == 0) {
        angle = -0.5f * atan2f(fp[2], fp[0]);
      }
      else {
        angle = 0.5f * atan2f(fp[0], fp[2]);
      }
    }
    else {
      if (upflag == 0) {
        angle = 0.5f * atan2f(-fp[1], -fp[0]);
      }
      else {
        angle = -0.5f * atan2f(-fp[0], -fp[1]);
      }
    }

    co = cosf(angle);
    si = sinf(angle) / len;
    q2 = glm::quat(co, tvec[0] * si, tvec[1] * si, tvec[2] * si);
    q = q2*q;
  }
}

glm::quat to_track_quat_ZY(glm::vec3 vec)
{
    glm::quat q;
    vec_to_quat(q, -vec, 2, 1);
    return q;
}

    struct SetUpBranchRetStruct
    {
        CHTurtle pos_tur;
        CHTurtle dir_tur;
        float radius_limit;
        float stem_offset;
    };

void WeberPennGenerator::Tree::init(WeberPennParametersNative _param, bool _generate_leaves)
{
    param = _param;
    generate_leaves = _generate_leaves;
    if (!generate_leaves)
        param.leaf_blos_num = 0;
}

void WeberPennGenerator::Tree::make()
{
    create_branches();
    //if (generate_leaves)
    //    create_leaf_mesh();
}

void WeberPennGenerator::Tree::create_branches()
{
    for (int level_depth = 0; level_depth < param.levels; level_depth++)
    {
        Curve level_curve;
        level_curve.resolution_u = param.curve_res[level_depth];
        branch_curves.push_back(level_curve);
    }
    std::vector<std::pair<vec3, float>> points;
    if (param.branches[0] > 0)
        points_for_floor_split(points);
    
    for (int ind = 0; ind<param.branches[0]; ind++)
    {
        tree_scale = param.g_scale + random_uniform(-1, 1) * param.g_scale_v;
        tree_scale = MAX(tree_scale, 0.1*param.g_scale);
        CHTurtle turtle = CHTurtle();
        turtle.pos = vec3(0,0,0);
        turtle.dir = vec3(0,0,1);
        turtle.right = vec3(1,0,0);
        if (param.branches[0] > 1)
        {
            //position randomly at base and rotate to face out
            auto &point = points[ind];
            turtle.roll_right(degrees(point.second) - 90);
            turtle.pos = point.first;
        }
        else
        {
            //start at random rotation
            turtle.roll_right(random_uniform(0, 360));
        }
        branch_curves[0].splines.data.emplace_back();
        branch_curves[0].splines.data.back().resolution_u = 1;
        branch_curves[0].splines.data.back().bezier_points[0].co = vec3(0,0,-0.1);
        branch_curves[0].splines.data.back().bezier_points[0].radius = 0.01;
        branch_curves[0].splines.data.back().bezier_points.emplace_back();
        branch_curves[0].splines.data.back().bezier_points[1].co = vec3(0,0,0);
        branch_curves[0].splines.data.back().bezier_points[1].radius = 0.01;
        stems.push_back(new Stem(0, branch_curves[0].splines.data.size()-1));
        root = stems.back();

        branch_curves[0].splines.data.emplace_back();
        auto &trunk = branch_curves[0].splines.data.back();
        trunk.resolution_u = param.curve_res[0];

        stems.push_back(new Stem(0, branch_curves[0].splines.data.size()-1));
        Stem *stem = stems.back();
        try
        {
            make_stem(turtle, *stem);
        }
        catch(const std::exception& e)
        {
            //branch_curves = {};
            //leaves_array = {};
            clear();   
            root = nullptr;
            //std::cerr << e.what() << '\n';
        }
        if (root)
            stem->parent = root;
        for (auto *s : stems)
        {
            s->children = {};
        }
        for (auto *s : stems)
        {
            //logerr("%uud s",s);
            //logerr("%uud",s->parent);
            if (s->copied_from)
            {
                s->copied_from->children.push_back(s);
            }
            else if (s->parent)
            {
                s->parent->children.push_back(s);
            }
        }

        /*
        int rec_points = 0;
        std::vector<Stem *> stk = {root};
        std::vector<Stem *> stk_n = {};
        if (root)
        {
            while (stk.size() > 0)
            {
                stk_n = {};
                for (Stem *s : stk)
                {
                    if (s->parent)
                        debug("stem (%d %d) parent (%d %d)\n", s->depth, s->spline_pos, s->parent->depth, s->parent->spline_pos);
                    else
                        debug("stem (%d %d) root\n", s->depth, s->spline_pos);
                    auto &c = branch_curves[s->depth].splines.data[s->spline_pos];
                    if (c.bezier_points.empty())
                        debug("curve is empty\n");
                    else
                    {
                        debug("curve %d points %f %f %f --> %f %f %f\n",c.bezier_points.size(),
                            c.bezier_points.front().co.x, c.bezier_points.front().co.y, c.bezier_points.front().co.z,
                            c.bezier_points.back().co.x, c.bezier_points.back().co.y, c.bezier_points.back().co.z);
                        for (auto &p : c.bezier_points)
                        {
                            debug("(%.4f %.4f %.4f)",p.co.x, p.co.y, p.co.z);
                        }
                        debugnl();
                    }
                    rec_points += c.bezier_points.size();
                    debug(" children {");
                    for (Stem *ch : s->children)
                    {
                        debug("(%d %d)", ch->depth, ch->spline_pos);
                        stk_n.push_back(ch);
                    }
                    debug("}\n");
                    s->radius_limit = -1000;
                }
                stk = stk_n;
            }
        }
        for (auto *s : stems)
        {
            //logerr("stem %u parent %u checked %d",s,s->parent,(int)s->radius_limit);
        }
        
        int curve_points = 0;
        for (auto &bl : branch_curves)
        {
            for (auto &c : bl.splines.data)
            {
                curve_points += c.bezier_points.size();
                for (auto &p : c.bezier_points)
                {
                    //logerr("%f %f %f",p.co.x, p.co.y, p.co.z);
                }
            }
        }
        logerr("created tree with %d/%d points", rec_points, curve_points);
        */
        int curve_points = 0;
        for (auto &bl : branch_curves)
        {
            for (auto &c : bl.splines.data)
            {
                curve_points += c.bezier_points.size();
                for (auto &p : c.bezier_points)
                {
                    //logerr("%f %f %f",p.co.x, p.co.y, p.co.z);
                }
            }
        }
        //logerr("created tree with %d points", curve_points);
    }
}

void WeberPennGenerator::Tree::create_leaf_mesh()
{
    if (leaves_array.empty())
        return;
    
    BaseLeafMesh base_leaf_shape;
    Leaf::get_shape(param.leaf_shape, tree_scale / param.g_scale,
                                           param.leaf_scale, param.leaf_scale_x, base_leaf_shape);
    BaseLeafMesh base_blossom_shape;
    Leaf::get_shape(-param.blossom_shape, tree_scale / param.g_scale, param.blossom_scale, 1, base_blossom_shape);
    logerr("leaf sz %f %f %f %f ",tree_scale, param.g_scale,
                                           param.leaf_scale, param.leaf_scale_x);
    std::vector<glm::vec3> leaf_verts;
    std::vector<std::vector<int>> leaf_faces;
    int leaf_index = 0;
    std::vector<glm::vec3> blossom_verts;
    std::vector<std::vector<int>> blossom_faces;
    int blossom_index = 0;
    int counter = 0;

    for (auto &leaf : leaves_array)
    {
        if (param.blossom_rate && random_uniform(0,1) < param.blossom_rate)
        {
            leaf.get_mesh(param.leaf_bend, base_blossom_shape, blossom_index,blossom_verts,blossom_faces);
            blossom_index++;
        }
        else
        {
            leaf.get_mesh(param.leaf_bend, base_leaf_shape, leaf_index,leaf_verts,leaf_faces);
            leaf_index++;
        }
        counter++;
    }
    //logerr("%d leaves added to mesh ", counter);
}

void WeberPennGenerator::Tree::points_for_floor_split(std::vector<std::pair<glm::vec3, float>> &points)
{
    points = {};
    tree_scale = param.g_scale + param.g_scale_v;
    Stem stem(0, 0);
    stem.length = calc_stem_length(stem);
    float rad = 2.5 * calc_stem_radius(stem);
    for (int i=0;i<param.branches[0];i++)
    {
        bool point_ok = false;
        while (!point_ok)
        {
            //distance from center proportional for number of splits, tree scale and stem radius
            float dis = sqrt(random_uniform(0.01,1) * param.branches[0] / 2.5 * param.g_scale * param.ratio); 
            //angle random in circle
            float theta = random_uniform(0, 2 * PI);
            vec3 pos = vec3(dis * cos(theta), dis * sin(theta), 0);
            //test point against those already in array to ensure it will not intersect
            bool point_m_ok = true;
            for (auto &point : points)
            {
                if (length(point.first - pos) < rad)
                {
                    point_m_ok = false;
                    break;
                }
            }
            if (point_m_ok)
            {
                point_ok = true;
                points.push_back(std::pair<glm::vec3, float>(pos, theta));
            }

        }
    }
}

float WeberPennGenerator::Tree::random_uniform(float from, float to)
{
    seed = seed * 1103515245 + 12345;
    float r =  ((unsigned int)(seed / 65536) % 32768) / 32768.0f;
    float res = from >= to ? from : from + r*(to - from);
    return res;
}

void WeberPennGenerator::Tree::make_stem(CHTurtle &turtle, Stem &stem, int start, float split_corr_angle, 
                                         float num_branches_factor, float clone_prob, CHTurtle *pos_corr_turtle, 
                                         CHTurtle *cloned_turtle)
{
    //logerr("make stem %d %d with turtle %f %f %f - %f %f %f st %d", stem.depth, stem.spline_pos, turtle.pos.x, turtle.pos.y,
    //       turtle.pos.z, turtle.dir.x, turtle.dir.y, turtle.dir.z, start);
    //logerr("r lim %f", stem.radius_limit);
    //if the stem is so thin as to be invisible then don't bother to make it
    if (stem.radius_limit >= 0 && stem.radius_limit < 0.0001)
        return;
    //use level 3 parameters for any depth greater than this
    int depth = stem.depth;
    int d_plus_1 = depth + 1;
    if (d_plus_1 > 3)
        d_plus_1 = 3;
    //calc length and radius for this stem (only applies for non clones)
    if (start == 0)     
    {       
        stem.length_child_max = param.length[d_plus_1] + random_uniform(-1, 1) * param.length_v[d_plus_1];
        stem.length_child_max = MAX(stem.length_child_max, 0.1*param.length[d_plus_1]);
        stem.length = calc_stem_length(stem);
        stem.radius = calc_stem_radius(stem);
        if (depth == 0)
            base_length = stem.length * param.base_size[0];
    }
    //if the branch origin needs to be repositioned so bevel doesnt sit outside parent
    if (pos_corr_turtle)
    {
        //pos_corr_turtle currently positioned on circumference so subtract this branch radius
        //to ensure open ends are never visible
        pos_corr_turtle->move(-MIN(stem.radius, stem.radius_limit));
        turtle.pos = pos_corr_turtle->pos;
    }
    //apply pruning, not required if is a clone, as this will have been tested already
    if (cloned_turtle == nullptr && param.prune_ratio > 0)
    {
        //save start length and random state
        float start_length = stem.length;
        unsigned long r_state = random_getstate();
        auto split_err_state = split_num_error;

        //iteratively scale length by 0.9 until it fits, or remove entirely if we get to 80% reduction
        CHTurtle copy_turtle = CHTurtle(turtle);
        bool in_pruning_envelope = test_stem(copy_turtle, stem, start, split_corr_angle, clone_prob);
        int cnt = 0;
        while (!in_pruning_envelope)
        {
            stem.length *= 0.9;
            if (stem.length <= 0.15 * start_length || cnt > 100)
            {
                //too short to look good so remove allow for semi prune with 0 length
                if (param.prune_ratio < 1)
                {
                    stem.length = 0;
                    break;
                }
                else
                    return;
            }
            random_setstate(r_state);
            split_num_error = split_err_state;
            CHTurtle ct = CHTurtle(turtle);
            in_pruning_envelope = test_stem(ct, stem, start, split_corr_angle, clone_prob);
            cnt++;
        }
        float fitting_length = stem.length;
        //apply reduction scaled by prune ratio
        stem.length = start_length * (1 - param.prune_ratio) + fitting_length * param.prune_ratio;
        //recalculate stem radius for new length
        stem.radius = calc_stem_radius(stem);
        //restore random state
        random_setstate(r_state);
        split_num_error = split_err_state;
    }
    //get parameters
    int curve_res = MAX(1,int(param.curve_res[depth]));
    float seg_splits = param.seg_splits[depth];
    float seg_length = stem.length / curve_res;

    //calc base segment
    int base_seg_ind = ceil(param.base_size[0] * int(param.curve_res[0]));

    int leaf_count = 0;
    int branch_count = 0;
    float f_branches_on_seg = 0;
    float f_leaves_on_seg = 0;
    if (depth == param.levels - 1 && depth > 0 && param.leaf_blos_num != 0)
    {
        //calc base leaf count
        leaf_count = calc_leaf_count(stem);
        //correct leaf count for start position along stem
        leaf_count *= 1 - start / (float)curve_res;
        //divide by curve_res to get no per seg
        f_leaves_on_seg = leaf_count / (float)curve_res;
    }
    else
    {
        //calc base branch count
        float f_branch_count = calc_branch_count(stem);
        //correct branch Count for start position along stem
        f_branch_count *= 1 - start /(float)curve_res;
        //correct for reduced number on clone branches
        f_branch_count *= num_branches_factor;
        //divide by curve_res to get no per seg
        f_branches_on_seg = f_branch_count /(float)curve_res;
        branch_count = f_branch_count;
    }
    //higher point resolution for flared based
    int max_points_per_seg = ceil(MAX(1.0, 100 /(float)curve_res));

    //set up FS error values
    float branch_num_error = 0;
    float leaf_num_error = 0;

    //decide on start rotation for branches/leaves
    //use array to allow other methods to modify the value (otherwise passed by value)
    std::vector<float> prev_rotation_angle = {0};//is it necessary to have vector here?
    if (param.rotate[d_plus_1] >= 0)
    {
        //start at random rotation
        prev_rotation_angle[0] = random_uniform(0, 360);
    }
    else
    {
        //on this case prev_rotation_angle used as multiplier to alternate side of branch
        prev_rotation_angle[0] = 1;
    }

    //calc helix parameters if needed
    glm::vec3 hel_p_0 = glm::vec3(0,0,0);
    glm::vec3 hel_p_1 = glm::vec3(0,0,0);
    glm::vec3 hel_p_2 = glm::vec3(0,0,0);
    glm::vec3 hel_axis = glm::vec3(1,0,0);

    if (param.curve_v[depth] < 0)
    {
        float tan_ang = tan(glm::radians(90 - abs(param.curve_v[depth])));
        float hel_pitch = 2 * stem.length / curve_res * random_uniform(0.8, 1.2);
        float hel_radius = 3 * hel_pitch / (16 * tan_ang) * random_uniform(0.8, 1.2);

        //apply full tropism if not trunk/main branch and horizontal tropism if is
        if (depth > 1)
            apply_tropism(turtle, param.tropism);
        else if (depth > 0)
            apply_tropism(turtle, vec3(param.tropism[0], param.tropism[1], 0));

        calc_helix_points(turtle, hel_radius, hel_pitch, hel_p_0, hel_p_1, hel_p_2, hel_axis);
    }
    
    //point resolution for this seg, max_points_per_seg if base, 1 otherwise
    int points_per_seg;
    if (depth == 0 || param.taper[depth] > 1)
        points_per_seg = max_points_per_seg;
    else
        points_per_seg = 2;
    
    for (int seg_ind=start;seg_ind<curve_res + 1;seg_ind++)
    {
        //logerr("creating stem %d, segments [%d %d]", depth, start,curve_res);
        int remaining_segs = curve_res + 1 - seg_ind;
    
        //set up next bezier point
        Point *new_point = nullptr;
        if (param.curve_v[depth] < 0)
        {
            //negative curve_v so helix branch
            glm::vec3 pos = turtle.pos;
            if (seg_ind == 0)
            {
                new_point = &(spline(stem).bezier_points[0]);
                new_point->co = pos;
                new_point->handle_right = hel_p_0 + pos;
                new_point->handle_left = pos;
                turtle.pos = new_point->co;
                turtle.dir = glm::normalize(new_point->handle_right);
            }
            else
            {
                spline(stem).add();
                total_points_cnt++;
                if (total_points_cnt > AbstractTreeGenerator::joints_limit)
                {
                    //logerr("too many joints");
                    throw std::exception();
                }
                new_point = &(spline(stem).bezier_points.back());

                if (seg_ind == 1)
                {
                    new_point->co = hel_p_2 + pos;
                    new_point->handle_left = hel_p_1 + pos;
                    new_point->handle_right = 2.0f * new_point->co - new_point->handle_left;
                }
                else
                {
                    Point &prev_point = spline(stem).bezier_points[spline(stem).bezier_points.size() - 2];
                    new_point->co = glm::angleAxis((seg_ind - 1) * PI, hel_axis) * hel_p_2;
                    new_point->co += prev_point.co;
                    glm::vec3 dif_p = glm::angleAxis((seg_ind - 1) * PI, hel_axis) * (hel_p_2 - hel_p_1);
                    new_point->handle_left = new_point->co - dif_p;
                    new_point->handle_right = 2.0f * new_point->co - new_point->handle_left;
                }
                turtle.pos = new_point->co;
                turtle.dir = glm::normalize(new_point->handle_right);
                if (stem.depth == 0)//TODO: remove me
                   turtle.dir = glm::vec3(0,0,1); 
            }
        }
        else
        {
            if (stem.depth == 0)//TODO: remove me
                   turtle.dir = glm::vec3(0,0,1); 
            //normal curved branch
            //get/make new point to be modified
            if (seg_ind != start)
            {
                turtle.move(seg_length);
                spline(stem).add();
                total_points_cnt++;
                if (total_points_cnt > AbstractTreeGenerator::joints_limit)
                {
                    //logerr("too many joints");
                    throw std::exception();
                }
            }
            new_point = seg_ind == start ? &(spline(stem).bezier_points[0]) : &(spline(stem).bezier_points.back());
            
            //set position and handles of new point
            //if this is a clone then correct initial direction to match original to make
            //split smoother
            new_point->co = turtle.pos;
            //if (stem.depth == 0)
            //    logerr("add pos %d %f %f %f seg len %f", stem.depth, turtle.dir.x, turtle.dir.y, turtle.dir.z, seg_length);
            if (cloned_turtle and seg_ind == start)
            {
                new_point->handle_left = turtle.pos - cloned_turtle->dir * (stem.length / (curve_res * 3));
                new_point->handle_right = turtle.pos + cloned_turtle->dir * (stem.length / (curve_res * 3));
            }
            else
            {
                new_point->handle_left = turtle.pos - turtle.dir * (stem.length / (curve_res * 3));
                new_point->handle_right = turtle.pos + turtle.dir * (stem.length / (curve_res * 3));
            }
        }

        //set radius of new point
        float actual_radius = radius_at_offset(stem, seg_ind /(float)curve_res);
        new_point->radius = actual_radius;
        int branches_on_seg = 0;
        int leaves_on_seg  = 0;
        if (seg_ind > start)
        {

            //calc number of splits at this seg (N/A for helix)
            int num_of_splits = 0;
            if (param.curve_v[depth] >= 0)
            {
                num_of_splits = 0;

                if (param.base_splits > 0 && depth == 0 && seg_ind == base_seg_ind)
                {
                    //if base_seg_ind and has base splits then override with base split number
                    //take random number of splits up to max of base_splits if negative
                    if (param.base_splits < 0)
                        num_of_splits = int(random_uniform(0, 1) * (abs(param.base_splits) + 0.5));
                    else
                        num_of_splits = int(param.base_splits);
                    num_of_splits = MIN(2,num_of_splits);
                }
                else if (seg_splits > 0 && seg_ind < curve_res && (depth > 0 || seg_ind > base_seg_ind))
                {
                    //otherwise get number of splits from seg_splits and use floyd-steinberg to
                    //fix non-integer values only clone with probability clone_prob
                    if (random_uniform(0, 1) <= clone_prob)
                    {
                        num_of_splits = int(seg_splits + split_num_error[depth]);
                        num_of_splits = MIN(2,num_of_splits);
                        split_num_error[depth] -= num_of_splits - seg_splits;

                        //reduce clone/branch propensity
                        clone_prob /= num_of_splits + 1;
                        num_branches_factor /= num_of_splits + 1;
                        num_branches_factor = MAX(0.8, num_branches_factor);

                        //TODO do this better?
                        //if depth != param.levels - 1:
                        branch_count *= num_branches_factor;
                        f_branches_on_seg = branch_count / (float)curve_res;
                    }
                }
            }
            //add branches/leaves for this seg
            //if below max level of recursion then draw branches, otherwise draw leaves
            auto r_state = random_getstate();
            if (abs(branch_count) > 0 and depth < param.levels - 1)
            {
                if (branch_count < 0)
                {
                    //fan branches
                    if (seg_ind == curve_res)
                        branches_on_seg = int(branch_count);
                    else
                        branches_on_seg = 0;
                }
                else
                {
                    //get FS corrected branch number
                    branches_on_seg = int(f_branches_on_seg + branch_num_error);
                    branch_num_error -= branches_on_seg - f_branches_on_seg;
                }
                //add branches
                //logerr("branch_count %d %d %d %d %d %d %f %f",branch_count,depth, param.levels,seg_ind,curve_res,branches_on_seg,
                //f_branches_on_seg,branch_num_error);
                if (abs(branches_on_seg) > 0)
                    make_branches(turtle, stem, seg_ind, branches_on_seg, prev_rotation_angle);
            }
            else if (abs(leaf_count) > 0 && depth > 0)
            {
                    if (leaf_count < 0)
                    {
                        //fan leaves
                        if (seg_ind == curve_res)
                            leaves_on_seg = leaf_count;
                        else
                            leaves_on_seg = 0;
                    }
                    else
                    {
                        //get FS corrected number of leaves
                        leaves_on_seg = int(f_leaves_on_seg + leaf_num_error);
                        leaf_num_error -= leaves_on_seg - f_leaves_on_seg;
                    }
                    //add leaves
                    //logerr("los %d",leaves_on_seg);
                    if (abs(leaves_on_seg) > 0)
                        make_leaves(turtle, stem, seg_ind, leaves_on_seg, prev_rotation_angle);
            }
            random_setstate(r_state);

            //perform cloning if needed, not allowed for helix (also don't curve/apply tropism as irrelevant)
                if (param.curve_v[depth] >= 0)
                {
                    if (num_of_splits > 0)
                    {
                        //calc angles for split
                        bool is_base_split = (param.base_splits > 0 and depth == 0 and seg_ind == base_seg_ind);
                        bool using_direct_split = param.split_angle[depth] < 0;
                        float spr_angle = 0;
                        float spl_angle = 0;
                        float split_corr_angle = 0;

                        if (using_direct_split)
                        {
                            spr_angle = abs(param.split_angle[depth]) + random_uniform(-1, 1) * 
                                        param.split_angle_v[depth];
                            spl_angle = 0;
                            split_corr_angle = 0;
                        }
                        else
                        {
                            float decl = declination(turtle.dir);//warning - degrees
                            spl_angle = param.split_angle[depth] + random_uniform(-1, 1) * param.split_angle_v[depth] - 
                                        decl;
                            spl_angle = MAX(0, spl_angle);
                            split_corr_angle = spl_angle / remaining_segs;
                            spr_angle = - (20 + 0.75 * (30 + abs(decl - 90) * pow(random_uniform(0,1),2)));
                        }
                        //make clone branches
                        r_state = random_getstate();
                        make_clones(turtle, seg_ind, split_corr_angle, num_branches_factor, clone_prob, stem,
                                    num_of_splits, spl_angle, spr_angle, is_base_split);
                        random_setstate(r_state);

                        //apply split to base stem
                        turtle.pitch_down(spl_angle / 2);

                        //apply spread if splitting to 2 and not base split
                        if (!is_base_split && num_of_splits == 1)
                        {
                            if (using_direct_split)
                            {
                                turtle.turn_right(spr_angle / 2);
                            }
                            else
                            {
                                turtle.dir = glm::angleAxis(glm::radians(-spr_angle / 2), vec3(0, 0, 1))*turtle.dir;
                                normalize(turtle.dir);
                                turtle.right = glm::angleAxis(glm::radians(-spr_angle / 2), vec3(0, 0, 1))*turtle.right;
                                normalize(turtle.right);
                            }
                        }
                    }
                    else
                    {
                        //just apply curve and split correction
                        turtle.turn_left(random_uniform(-1, 1) * param.bend_v[depth] / curve_res);
                        float curve_angle = calc_curve_angle(depth, seg_ind);
                        turtle.pitch_down(curve_angle - split_corr_angle);
                    }
                    //apply full tropism if not trunk/main branch and horizontal tropism if is
                    if (depth > 1)
                        apply_tropism(turtle, param.tropism);
                    else
                        apply_tropism(turtle, vec3(param.tropism[0], param.tropism[1], 0));
                }
                //increase point resolution at base of trunk and apply flaring effect
                if (points_per_seg > 2)
                    increase_bezier_point_res(stem, seg_ind, points_per_seg);
        }
    }
    //scale down bezier point handles for flared base of trunk
    if (points_per_seg > 2)
        scale_bezier_handles_for_flare(stem, max_points_per_seg);
    stem_index++;
}

float WeberPennGenerator::Tree::calc_stem_length(Stem &stem)
{
    float result = 0;
    //Calculate length of this stem as defined in paper
    if (stem.depth == 0) // trunk
    {
        float len = (param.length[0] + random_uniform(-1, 1) * param.length_v[0]);
        len = MAX(len, 0.1*param.length[0]);
        result = tree_scale * len;
        trunk_length = result;
    }
    else if (stem.depth == 1) //first level
    {
        result = stem.parent->length * stem.parent->length_child_max * 
                 shape_ratio(param.shape, (stem.parent->length - stem.offset) / (stem.parent->length - base_length));
    }
    else
        result = stem.parent->length_child_max * (stem.parent->length - 0.7 * stem.offset);
    return MAX(1e-4, result);
}

float WeberPennGenerator::Tree::calc_stem_radius(Stem &stem)
{
    float result = 0;
    //Calculate radius of this stem as defined in paper
    if (stem.depth == 0) // trunk
        result = stem.length * param.ratio * param.radius_mod[0];
    else
    {
        result = param.radius_mod[stem.depth] * stem.parent->radius * 
                 pow((stem.length / stem.parent->length), param.ratio_power);
        result = MAX(0.005, result);
        result = MIN(stem.radius_limit, result);
    }
    return result;
}

bool WeberPennGenerator::Tree::test_stem(CHTurtle &turtle, Stem &stem, int start, float split_corr_angle, 
                                         float clone_prob)
{
    //Test if stem is inside pruning envelope
    //use level 3 parameters for any depth greater than this
    int depth = stem.depth;
    int d_plus_1 = depth + 1;
    if (d_plus_1 > 3)
        d_plus_1 = 3;
    //get parameters
    int curve_res = int(param.curve_res[depth]);
    float seg_splits = param.seg_splits[depth];
    float seg_length = stem.length / curve_res;

    //calc base segment
    int base_seg_ind = ceil(param.base_size[0] * int(param.curve_res[0]));

    //decide on start rotation for branches/leaves
    //use array to allow other methods to modify the value (otherwise passed by value)
    std::vector<float> prev_rotation_angle = {0};//is it necessary to have vector here?
    if (param.rotate[d_plus_1] >= 0)
    {
        //start at random rotation
        prev_rotation_angle[0] = random_uniform(0, 360);
    }
    else
    {
        //on this case prev_rotation_angle used as multiplier to alternate side of branch
        prev_rotation_angle[0] = 1;
    }

    //calc helix parameters if needed
    glm::vec3 hel_p_0 = glm::vec3(0,0,0);
    glm::vec3 hel_p_1 = glm::vec3(0,0,0);
    glm::vec3 hel_p_2 = glm::vec3(0,0,0);
    glm::vec3 hel_axis = glm::vec3(1,0,0);

    if (param.curve_v[depth] < 0)
    {
        float tan_ang = tan(glm::radians(90 - abs(param.curve_v[depth])));
        float hel_pitch = 2 * stem.length / curve_res * random_uniform(0.8, 1.2);
        float hel_radius = 3 * hel_pitch / (16 * tan_ang) * random_uniform(0.8, 1.2);

        //apply full tropism if not trunk/main branch and horizontal tropism if is
        if (depth > 1)
            apply_tropism(turtle, param.tropism);
        else
            apply_tropism(turtle, vec3(param.tropism[0], param.tropism[1], 0));

        calc_helix_points(turtle, hel_radius, hel_pitch, hel_p_0, hel_p_1, hel_p_2, hel_axis);
    }

    glm::vec3 previous_helix_point = glm::vec3(0,0,0);
    for (int seg_ind = start; seg_ind< curve_res + 1;seg_ind++)
    {
        int remaining_segs = curve_res + 1 - seg_ind;

        //set up next bezier point
        if (param.curve_v[depth] < 0)
        {
            //negative curve_v so helix branch
            vec3 pos = turtle.pos;
            if (seg_ind == 0)
                turtle.pos = pos;
            else
            {
                if (seg_ind == 1)
                    turtle.pos = hel_p_2 + pos;
                else
                {
                    hel_p_2 = glm::angleAxis((seg_ind - 1)*PI, hel_axis)*hel_p_2;
                    turtle.pos = hel_p_2 + previous_helix_point;
                }
            }
            previous_helix_point = turtle.pos;
        }
        else
        {
            //normal curved branch
            //move turtle
            if (seg_ind != start)
            {
                turtle.move(seg_length);
                if  (!(stem.depth == 0 && start < base_seg_ind) && !point_inside(turtle.pos))
                    return false;
            }
        }
        if (seg_ind > start)
        {
            //calc number of splits at this seg (N/A for helix)
            if (param.curve_v[depth] >= 0)
            {
                int num_of_splits = 0;
                if (param.base_splits > 0 && depth == 0 && seg_ind == base_seg_ind)
                {
                    //if base_seg_ind and has base splits then override with base split number
                    //take random number of splits up to max of base_splits
                    num_of_splits = int(random_uniform(0,1) * (param.base_splits + 0.5));
                }
                else if (seg_splits > 0 && seg_ind < curve_res && (depth > 0 || seg_ind > base_seg_ind))
                {
                   //otherwise get number of splits from seg_splits and use Floyd-Steinberg to
                   //fix non-integer values only clone with probability clone_prob 
                   if (random_uniform(0,1) < clone_prob)
                   {
                       num_of_splits = int(seg_splits + split_num_error[depth]);
                       split_num_error[depth] -= num_of_splits - seg_splits;
                       //reduce clone/branch propensity
                       clone_prob /= num_of_splits + 1;
                   }
                }

                //perform cloning if needed, not allowed for helix (also don't curve/apply tropism as irrelevant)
                if (num_of_splits > 0)
                {
                    //calc angles for split  
                    bool is_base_split = (param.base_splits > 0 && depth == 0 && seg_ind == base_seg_ind); 
                    bool using_direct_split = param.split_angle[depth] < 0;
                    float spr_angle = 0;
                    float spl_angle = 0;
                    float split_corr_angle = 0;

                    if (using_direct_split)
                    {
                        spr_angle = abs(param.split_angle[depth]) + random_uniform(-1, 1) * param.split_angle_v[depth];
                        spl_angle = 0;
                        split_corr_angle = 0;
                    }
                    else
                    {
                        float decl = declination(turtle.dir);
                        spl_angle = param.split_angle[depth] + random_uniform(-1, 1) * param.split_angle_v[depth] - decl;
                        spl_angle = MAX(0, spl_angle);
                        split_corr_angle = spl_angle / remaining_segs;
                        spr_angle = - (20 + 0.75 * (30 + abs(decl - 90) * pow(random_uniform(0,1),2)));
                    }
                    //apply split to base stem
                    turtle.pitch_down(spl_angle / 2);
                    // apply spread if splitting to 2 and not base split
                    if (!is_base_split && num_of_splits == 1)
                    {
                        if (using_direct_split)
                        {
                            turtle.turn_right(spr_angle / 2);
                        }
                        else
                        {
                            turtle.dir = glm::angleAxis(glm::radians(-spr_angle / 2), vec3(0, 0, 1)) * turtle.dir;
                            normalize(turtle.dir);
                            turtle.right = glm::angleAxis(glm::radians(-spr_angle / 2), vec3(0, 0, 1)) * turtle.right;
                            normalize(turtle.right);
                        }
                    }
                }
                else
                {
                    //just apply curve and split correction
                    turtle.turn_left(random_uniform(-1, 1) * param.bend_v[depth] / curve_res);
                    float curve_angle = calc_curve_angle(depth, seg_ind);
                    turtle.pitch_down(curve_angle - split_corr_angle);                
                }
                //apply full tropism if not trunk/main branch and horizontal tropism if is
                if (depth > 1)
                    apply_tropism(turtle, param.tropism);
                else
                    apply_tropism(turtle, glm::vec3(param.tropism[0], param.tropism[1], 0));
            }
        }
    }

    return point_inside(turtle.pos);
}

int WeberPennGenerator::Tree::calc_leaf_count(Stem &stem)
{
    // Calculate leaf count of this stem as defined in paper
    int result = 0;
    if (param.leaf_blos_num >= 0)
    {
        // scale number of leaves to match global scale and taper
        float leaves = param.leaf_blos_num * tree_scale / param.g_scale;
        float div = (stem.parent->length_child_max * stem.parent->length);
        
        result = div > 1e-4 ? leaves * (stem.length / div) : leaves;
        //logerr("l cnt = %d %f %f %f %f %d",param.leaf_blos_num , tree_scale , param.g_scale,
        //stem.parent->length_child_max, stem.parent->length, result);
    }
    else // fan leaves
        return param.leaf_blos_num;
    return result;
}

float WeberPennGenerator::Tree::calc_branch_count(Stem &stem)
{
    // Calculate branch count of this stem as defined in paper
    int d_p_1 = MIN(stem.depth + 1, 3);
    float result = 0;
    if (stem.depth == 0)
        result = param.branches[d_p_1] * (random_uniform(0, 1) * 0.2 + 0.9);
    else
    {
        if (param.branches[d_p_1] < 0)
            result = param.branches[d_p_1];
        else if (stem.depth == 1)
            result = param.branches[d_p_1] * (0.2 + 0.8 * (stem.length / stem.parent->length) / stem.parent->length_child_max);
        else
            result = param.branches[d_p_1] * (1.0 - 0.5 * stem.offset / stem.parent->length);
    }
    return result / (1 - param.base_size[stem.depth]);
}

void WeberPennGenerator::Tree::apply_tropism(CHTurtle &turtle, glm::vec3 tropism_vector)
{
    //Apply tropism_vector to turtle direction
    auto h_cross_t = glm::cross(turtle.dir, tropism_vector);
    //calc angle to rotate by (from ABoP) multiply to achieve accurate results from WP attractionUp param
    float alpha = 10 * length(h_cross_t);
    normalize(h_cross_t);
    //rotate by angle about axis perpendicular to turtle direction and tropism vector
    turtle.dir = glm::normalize(glm::angleAxis(glm::radians(alpha), h_cross_t) * turtle.dir);
    turtle.right = glm::normalize(glm::angleAxis(glm::radians(alpha), h_cross_t) * turtle.right);
}

void WeberPennGenerator::Tree::calc_helix_points(CHTurtle &turtle, float rad, float pitch, 
                                                 glm::vec3 &hel_p_0, glm::vec3 &hel_p_1, glm::vec3 &hel_p_2, 
                                                 glm::vec3 &hel_axis)
{
    //calculates required points to produce helix bezier curve with given radius and pitch in direction of turtle
    // alpha = radians(90)
    // pit = pitch/(2*PI)
    // a_x = rad*cos(alpha)
    // a_y = rad*sin(alpha)
    // a = pit*alpha*(rad - a_x)*(3*rad - a_x)/(a_y*(4*rad - a_x)*tan(alpha))
    // b_0 = Vector([a_x, -a_y, -alpha*pit])
    // b_1 = Vector([(4*rad - a_x)/3, -(rad - a_x)*(3*rad - a_x)/(3*a_y), -a])
    // b_2 = Vector([(4*rad - a_x)/3, (rad - a_x)*(3*rad - a_x)/(3*a_y), a])
    // b_3 = Vector([a_x, a_y, alpha*pit])
    // axis = Vector([0, 0, 1])

    // simplifies greatly for case inc_angle = 90
    vec3 points[4] = {vec3(0, -rad, -pitch / 4),
                           vec3((4 * rad) / 3, -rad, 0),
                           vec3((4 * rad) / 3, rad, 0),
                           vec3(0, rad, pitch / 4)};

    // align helix points to turtle direction and randomize rotation around axis
    auto trf = to_track_quat_ZY(turtle.dir);
    float spin_ang = random_uniform(0, 2 * PI);
    auto rot_quat = glm::angleAxis(glm::radians(spin_ang), vec3(0, 0, 1));

    for (auto &p :points)
    {
        p = rot_quat*p;
        p = trf*p;
    }
    hel_p_0 =  points[1] - points[0];
    hel_p_1 =  points[2] - points[0];
    hel_p_2 =  points[3] - points[0];
    hel_axis =  turtle.dir;
}

float WeberPennGenerator::Tree::radius_at_offset(Stem &stem, float z_1)
{
    // calculate radius of stem at offset z_1 along it
    float n_taper = param.taper[stem.depth];
    float unit_taper = 0;
    float radius = 0;
    if (n_taper < 1)
        unit_taper = n_taper;
    else if (n_taper < 2)
        unit_taper = 2 - n_taper;
    else
        unit_taper = 0;
    float taper = stem.radius * (1 - unit_taper * z_1);

    if (n_taper < 1)
        radius = taper;
    else
    {
        float z_2 = (1 - z_1) * stem.length;
        float z_3 = 0;
        float depth = 0;
        if (n_taper < 2 || z_2 < taper)
            depth = 1;
        else
            depth = n_taper - 2;
        if (n_taper < 2)
            z_3 = z_2;
        else
            z_3 = abs(z_2 - 2 * taper * int(z_2 / (2 * taper) + 0.5));
        if (n_taper < 2 && z_3 >= taper)
            radius = taper;
        else
            radius = (1 - depth) * taper + depth * sqrt(pow(taper, 2) - pow((z_3 - taper), 2));
    }
    if (stem.depth == 0)
    {
        float y_val = MAX(0, 1 - 8 * z_1);
        float flare = param.flare * ((pow(100, y_val) - 1) / 100) + 1;
        radius *= flare;
    }
    return CLAMP(radius,0,1000);
}

void WeberPennGenerator::Tree::make_branches(CHTurtle &turtle, Stem &stem, int seg_ind, int branches_on_seg,
                                             std::vector<float> &prev_rotation_angle, bool is_leaves)
{
    // Make the required branches for a segment of the stem
    //logerr("make branches %d", branches_on_seg);
    Point &start_point = spline(stem).bezier_points[spline(stem).bezier_points.size() - 2];
    Point &end_point = spline(stem).bezier_points[spline(stem).bezier_points.size() - 1];
    std::vector<SetUpBranchRetStruct> branches_array = {};
    int d_plus_1 = MIN(3, stem.depth + 1);

    if (branches_on_seg < 0) // fan branches
    {
        //logerr("bos %d",branches_on_seg);
        for (int branch_ind = 0; branch_ind < abs(int(branches_on_seg)); branch_ind++)
        {
            int stem_offset = 1;
            //logerr("%d %d",branch_ind,branches_on_seg);
            branches_array.push_back(set_up_branch(turtle, stem, BranchMode::fan, 1, start_point, end_point,
                                                   stem_offset, branch_ind, prev_rotation_angle, abs(branches_on_seg)));
        }
    }
    else
    {
        base_length = stem.length * param.base_size[stem.depth];
        float branch_dist = param.branch_dist[d_plus_1];
        int curve_res = int(param.curve_res[stem.depth]);

        if (branch_dist > 1) // whorled branches
        {
            // calc number of whorls, will result in a rounded number of branches rather than the
            // exact amount specified by branches_on_seg
            int num_of_whorls = int(branches_on_seg / (branch_dist + 1));
            float branches_per_whorl = branch_dist + 1;
            float branch_whorl_error = 0;

            for (int whorl_num = 0; whorl_num < num_of_whorls; whorl_num++)
            {
                // calc whorl offset in segment and on stem
                float offset = MIN(MAX(0.0, whorl_num / (float)num_of_whorls), 1.0);
                float stem_offset = (((seg_ind - 1) + offset) / curve_res) * stem.length;

                // if not in base area then make the branches
                if (stem_offset > base_length)
                {
                    // calc FS corrected num of branches this whorl
                    int branches_this_whorl = int(branches_per_whorl + branch_whorl_error);
                    branch_whorl_error -= branches_this_whorl - branches_per_whorl;

                    // set up these branches
                    for (int branch_ind = 0; branch_ind < branches_this_whorl; branch_ind++)
                    {
                        branches_array.push_back(set_up_branch(turtle, stem, BranchMode::whorled, offset, start_point,
                                                               end_point, stem_offset, branch_ind, prev_rotation_angle,
                                                               branches_this_whorl));
                    }
                }
                // rotate start angle for next whorl
                prev_rotation_angle[0] += param.rotate[d_plus_1];
            }
        }
        else // alternating or opposite branches
        {
            // ensure even number of branches on segment if near opposite
            for (int branch_ind = 0; branch_ind < branches_on_seg; branch_ind++)
            {
                // calc offset in segment and on stem
                float offset = 0;
                if (branch_ind % 2 == 0)
                    offset = MIN(MAX(0, branch_ind / (float)branches_on_seg), 1);
                else
                    offset = MIN(MAX(0, (branch_ind - branch_dist) / (float)branches_on_seg), 1);

                float stem_offset = (((seg_ind - 1) + offset) / curve_res) * stem.length;
                //logerr("%d %d %f off %f",branch_ind,branches_on_seg,branch_dist,offset);
                // if not in base area then set up the branch
                if (stem_offset > base_length)
                {
                    branches_array.push_back(set_up_branch(turtle, stem, BranchMode::alt_opp, offset,
                                                           start_point, end_point, stem_offset, branch_ind,
                                                           prev_rotation_angle));
                }
            }
        }
    }
    // make all new branches from branches_array, passing pos_corr_turtle which will be used to
    // set the position of branch_turtle in this call
    if (is_leaves)
    {
        for (auto &b : branches_array)
        {
            if (random_uniform(0, 1) < param.leaf_rate)
            {
                stem.leaves.push_back(leaves_array.size());
                leaves_array.push_back(Leaf(b.pos_tur.pos, b.dir_tur.dir, b.dir_tur.right));
                //logerr("too many leaves");
                if (leaves_array.size() > 10*AbstractTreeGenerator::joints_limit)
                    throw std::exception();
            }
        }
    }
    else
    {
        for (auto &b : branches_array)
        {
            branch_curves[d_plus_1].splines.data.emplace_back();
            Spline &new_spline = branch_curves[d_plus_1].splines.data.back();
            new_spline.resolution_u = param.curve_res[d_plus_1];
            stems.push_back(new Stem(d_plus_1, branch_curves[d_plus_1].splines.data.size()-1, &stem, b.stem_offset, b.radius_limit));
            Stem *new_stem = stems.back();
            //logerr("created stem (%d %d)", new_stem->depth, new_stem->spline_pos);
            make_stem(b.dir_tur, *new_stem, 0, 0, 1, 1, &b.pos_tur, nullptr);
        }
    }
}

void WeberPennGenerator::Tree::make_leaves(CHTurtle &turtle, Stem &stem, int seg_ind, int leaves_on_seg, 
                                           std::vector<float> &prev_rotation_angle)
{
    make_branches(turtle, stem, seg_ind, leaves_on_seg,prev_rotation_angle, true);
}

float WeberPennGenerator::Tree::calc_curve_angle(int depth, int seg_ind)
{
    // Calculate curve angle for segment number seg_ind on a stem
    float curve = param.curve[depth];
    float curve_v = param.curve_v[depth];
    float curve_back = param.curve_back[depth];
    int curve_res = int(param.curve_res[depth]);
    float curve_angle = 0;
    if (curve_back == 0)
        curve_angle = curve / curve_res;
    else
    {
        if (seg_ind < curve_res / 2.0)
            curve_angle = curve / (curve_res / 2.0);
        else
            curve_angle = curve_back / (curve_res / 2.0);
    }
    curve_angle += random_uniform(-1, 1) * (curve_v / curve_res);
    return curve_angle;
}

void WeberPennGenerator::Tree::make_clones(CHTurtle &turtle, int seg_ind, float split_corr_angle, float num_branches_factor,
                                           float clone_prob, Stem &stem, int num_of_splits, float spl_angle, float spr_angle,
                                           bool is_base_split)
{
    // make clones of branch used if seg_splits or base_splits > 0

    bool using_direct_split = param.split_angle[stem.depth] < 0;
    float stem_depth = param.split_angle_v[stem.depth];

    if (!is_base_split && num_of_splits > 2 && using_direct_split)
    {
        logerr("Only splitting up to 3 branches is supported");
        return;
    }
    for (int split_index = 0; split_index < num_of_splits; split_index++)
    {
        // copy turtle for new branch
        CHTurtle n_turtle = CHTurtle(turtle);
        // tip branch down away from axis of stem
        n_turtle.pitch_down(spl_angle / 2);

        // spread out clones
        float eff_spr_angle = 0;
        if (is_base_split && !using_direct_split)
        {
            eff_spr_angle = (split_index + 1) * (360.0 / (num_of_splits + 1)) + random_uniform(-1, 1) * stem_depth;
        }
        else
        {
            if (split_index == 0)
                eff_spr_angle = spr_angle / 2;
            else
                eff_spr_angle = -spr_angle / 2;
        }
        if (using_direct_split)
            n_turtle.turn_left(eff_spr_angle);

        else
        {
            auto quat = glm::angleAxis(glm::radians(eff_spr_angle), vec3(0, 0, 1));

            n_turtle.dir = glm::normalize(quat * n_turtle.dir);
            n_turtle.right = glm::normalize(quat * n_turtle.right);
        }
        // create new clone branch and set up then recurse
        branch_curves[stem.depth].splines.data.emplace_back();
        auto &split_stem = branch_curves[stem.depth].splines.data.back();
        split_stem.resolution_u = spline(stem).resolution_u;
        //split_stem.bezier_points = spline(stem).bezier_points;
        //split_stem.
        Stem *new_stem = new Stem(stem);
        stems.push_back(new_stem);
        new_stem->spline_pos = branch_curves[stem.depth].splines.data.size()-1;
        //logerr("copied stem (%d %d) from (%d %d)", new_stem->depth, new_stem->spline_pos,
        //stem.depth, stem.spline_pos);
        CHTurtle *cloned = nullptr;
        if (param.split_angle_v[stem.depth] >= 0)
            cloned = &turtle;

        make_stem(n_turtle, *new_stem, seg_ind, split_corr_angle, num_branches_factor, clone_prob,
                  nullptr, cloned);
    }
}

void WeberPennGenerator::Tree::increase_bezier_point_res(Stem &stem, int seg_ind, int points_per_seg)
{
    // need a copy of the end point as it is moved during the process, but also used for calculations throughout
    int curve_res = int(param.curve_res[stem.depth]);

    Point &seg_end_point = spline(stem).bezier_points[spline(stem).bezier_points.size() - 1];
    Point end_point = seg_end_point;
    Point &seg_start_point = spline(stem).bezier_points[spline(stem).bezier_points.size() - 2];
    Point start_point = seg_start_point;
    for (int k = 0; k < points_per_seg; k++)
    {
        // add new point and position
        // at this point the normals are left over-sized in order to allow for evaluation of the
        // original curve in later steps
        // once the stem is entirely built we then go back and scale the handles
        float offset = k / (float)(points_per_seg - 1);
        Point *curr_point = nullptr;
        if (k == 0)
            curr_point = &seg_start_point;
        else
        {
            if (k == 1)
                curr_point = &seg_end_point;
            else
            {
                spline(stem).add();
                total_points_cnt++;
                if (total_points_cnt > AbstractTreeGenerator::joints_limit)
                {
                    //logerr("too many joints");
                    throw std::exception();
                }
                curr_point = &(spline(stem).bezier_points[spline(stem).bezier_points.size() - 1]);
            }   
            if (k == points_per_seg - 1)
            {
                curr_point->co = end_point.co;
                curr_point->handle_left = end_point.handle_left;
                curr_point->handle_right = end_point.handle_right;
                //if (stem.depth == 0)
                //    logerr("add pos 1 %f %f %f", curr_point->co.x, curr_point->co.y, curr_point->co.z);
            }
            else
            {
                curr_point->co = calc_point_on_bezier(offset, start_point, end_point);
                // set handle to match direction of curve
                vec3 tangent = normalize(calc_tangent_to_bezier(offset, start_point, end_point));
                // and set the magnitude to match other control points
                float dir_vec_mag = length(end_point.handle_left - end_point.co);
                curr_point->handle_left = curr_point->co - tangent * dir_vec_mag;
                curr_point->handle_right = curr_point->co + tangent * dir_vec_mag;
                //if (stem.depth == 0)
                //{
                //    logerr("add pos 2 %f %f %f OFF %f", curr_point->co.x, curr_point->co.y, curr_point->co.z, offset);
                //    logerr("stp %f %f %f", start_point.co.x,start_point.co.y,start_point.co.z);
                //    logerr("enp %f %f %f", end_point.co.x,end_point.co.y,end_point.co.z);
                //}
            }
        }
        curr_point->radius = radius_at_offset(stem, (offset + seg_ind - 1) / curve_res);
    }
}

void WeberPennGenerator::Tree::scale_bezier_handles_for_flare(Stem &stem, int max_points_per_seg)
{
    //Reduce length of bezier handles to account for increased density of points on curve for
    //flared base of trunk
    for (auto &point : spline(stem).bezier_points)
    {
        point.handle_left = point.co + (point.handle_left - point.co) / (float)max_points_per_seg;
        point.handle_right = point.co + (point.handle_right - point.co) / (float)max_points_per_seg;
    }
}

float WeberPennGenerator::Tree::shape_ratio(int shape, float ratio)
{
    // Calculate shape ratio as defined in paper
    float result = 0;
    if (shape == 1) // spherical
        result = 0.2 + 0.8 * sin(PI * ratio);
    else if (shape == 2) // hemispherical
        result = 0.2 + 0.8 * sin(0.5 * PI * ratio);
    else if (shape == 3) // cylindrical
        result = 1.0;
    else if (shape == 4) // tapered cylindrical
        result = 0.5 + 0.5 * ratio;
    else if (shape == 5) // flame
    {
        if (ratio <= 0.7)
            result = ratio / 0.7;
        else
            result = (1.0 - ratio) / 0.3;
    }
    else if (shape == 6) // inverse conical
        result = 1.0 - 0.8 * ratio;
    else if (shape == 7) // tend flame
    {
        if (ratio <= 0.7)
            result = 0.5 + 0.5 * ratio / 0.7;
        else
            result = 0.5 + 0.5 * (1.0 - ratio) / 0.3;
    }
    else if (shape == 8) // envelope
    {
        if (ratio < 0 || ratio > 1)
            result = 0.0;
        else if (ratio < 1 - param.prune_width_peak)
            result = pow(ratio / (1 - param.prune_width_peak), param.prune_power_high);
        else
            result = pow((1 - ratio) / (1 - param.prune_width_peak), param.prune_power_low);
    }
    else // conical (0)
        result = 0.2 + 0.8 * ratio;
    return result;
}

bool WeberPennGenerator::Tree::point_inside(glm::vec3 point)
{
    // Check if point == inside pruning envelope, from WP 4.6
    // return point_in_cube(Vector([point.x, point.y, point.z - base_length]))
    float dist = sqrt(SQR(point.x) + SQR(point.y));
    float ratio = (tree_scale - point.z) / (tree_scale * (1 - param.base_size[0]));
    bool inside = (dist / tree_scale) < (param.prune_width * shape_ratio(8, ratio));
    // inside = inside and (point.x > -0.7 or point.z > 5.3)
    return inside;
}

SetUpBranchRetStruct
WeberPennGenerator::Tree::set_up_branch(CHTurtle &turtle, Stem &stem, BranchMode branch_mode, float offset,
                                        Point &start_point, Point &end_point, float stem_offset, int branch_ind,
                                        std::vector<float> &prev_rot_ang, int branches_in_group)
{
    if (!std::isfinite(start_point.co.x))
    {
        throw std::exception();
    }
    //Set up a new branch, creating the new direction and position turtle and orienting them
    //correctly and adding the required info to the list of branches to be made

        int d_plus_1 = MIN(3, stem.depth + 1);

        //make branch direction turtle
        auto branch_dir_turtle = make_branch_dir_turtle(turtle, param.curve_v[stem.depth] < 0, offset, start_point,
                                                        end_point);
        float radius_limit = 0;
        //calc rotation angle
        if (branch_mode == BranchMode::fan)
        {
            float t_angle = 0;
            if (branches_in_group == 1)
                t_angle = 0;
            else
            {
                t_angle = (param.rotate[d_plus_1] * (
                    (branch_ind / (branches_in_group - 1.0)) - 1.0 / 2.0)) + random_uniform(-1, 1) * param.rotate_v[
                    d_plus_1];
            }
            branch_dir_turtle.turn_right(t_angle);
            radius_limit = 0;
        }
        else
        {
            float r_angle = 0;
            if (branch_mode == BranchMode::whorled)
            {
                r_angle = prev_rot_ang[0] + (360 * branch_ind / (float)branches_in_group) + random_uniform(-1, 1) *
                    param.rotate_v[d_plus_1];
            }
            else
            {
                r_angle = calc_rotate_angle(d_plus_1, prev_rot_ang[0]);
                if (param.rotate[d_plus_1] >= 0)
                    prev_rot_ang[0] = r_angle;
                else
                    prev_rot_ang[0] = -prev_rot_ang[0];
            }
            //orient direction turtle to correct rotation
            branch_dir_turtle.roll_right(r_angle);
            radius_limit = radius_at_offset(stem, stem_offset / MAX(stem.length,1e-6));
        }
        //make branch position turtle in appropriate position on circumference
        auto branch_pos_turtle = make_branch_pos_turtle(branch_dir_turtle, offset, start_point,
                                                        end_point, radius_limit);

        //calc down angle
        float d_angle = calc_down_angle(stem, stem_offset);

        //orient direction turtle to correct declination
        branch_dir_turtle.pitch_down(d_angle);

        //return branch info
        SetUpBranchRetStruct res;
        res.pos_tur = branch_pos_turtle;
        res.dir_tur = branch_dir_turtle;
        res.radius_limit = radius_limit;
        res.stem_offset = stem_offset;
        return res;
}

glm::vec3 WeberPennGenerator::Tree::calc_point_on_bezier(float offset, Point &start_point, Point &end_point)
{
    //Evaluate Bezier curve at offset between bezier_spline_points start_point and end_point
    if (offset < 0 || offset > 1)
    {
        logerr("Offset out of range: %f not between 0 and 1", offset);
        return vec3(0,0,0);
    }

    float one_minus_offset = 1 - offset;

    vec3 res = one_minus_offset*one_minus_offset*one_minus_offset * start_point.co + 
              3 * one_minus_offset*one_minus_offset * offset * start_point.handle_right + 
              3 * one_minus_offset * offset*offset * end_point.handle_left + 
              offset*offset*offset * end_point.co;

    return res;
}

glm::vec3 WeberPennGenerator::Tree::calc_tangent_to_bezier(float offset, Point &start_point, Point &end_point)
{
    //Calculate tangent to Bezier curve at offset between bezier_spline_points start_point and end_point
    if (offset < 0 || offset > 1)
    {
        logerr("Offset out of range: %f not between 0 and 1", offset);
        return vec3(1,0,0);
    }
    if (length(start_point.co - end_point.co)<1e-6)
    {
        return vec3(1,0,0);
    }
    float one_minus_offset = 1 - offset;
    vec3 start_handle_right = start_point.handle_right;
    vec3 end_handle_left = end_point.handle_left;
    vec3 res = 3 * one_minus_offset*one_minus_offset * (start_handle_right - start_point.co) + 
               6 * one_minus_offset * offset *(end_handle_left - start_handle_right) + 
               3 * offset*offset * (end_point.co - end_handle_left);
    
    return res;
}

CHTurtle WeberPennGenerator::Tree::make_branch_dir_turtle(CHTurtle &turtle, bool helix, float offset, Point &start_point,
                                                          Point &end_point)
{
    //Create and setup the turtle for the direction of a new branch

    CHTurtle branch_dir_turtle = CHTurtle();
    vec3 tangent = calc_tangent_to_bezier(offset, start_point, end_point);
    tangent = normalize(tangent);
    branch_dir_turtle.dir = tangent;

    if (helix)
    {
        //approximation to actual normal to preserve for helix
        vec3 tan_d = normalize(calc_tangent_to_bezier(CLAMP(offset + 0.0001,0,1), start_point, end_point));
        branch_dir_turtle.right = glm::cross(branch_dir_turtle.dir,tan_d);
    }
    else
    {
        //generally curve lines in plane define by turtle.right, so is fair approximation to take new right as being
        //parallel to this, ie find the turtle up vector (in the plane) and cross with tangent (assumed in the plane)
        //to get the new direction - this doesn't hold for the helix
        branch_dir_turtle.right = glm::cross(glm::cross(turtle.dir, turtle.right), branch_dir_turtle.dir);
    }
    if (!std::isfinite(branch_dir_turtle.dir.x))
    {
        throw std::exception();
    }
    return branch_dir_turtle;
}

CHTurtle WeberPennGenerator::Tree::make_branch_pos_turtle(CHTurtle &dir_turtle, float offset, Point &start_point,
                                                          Point &end_point, float radius_limit)
{
    //Create and setup the turtle for the position of a new branch, also returning the radius
    //of the parent to use as a limit for the child

    //logerr("offset %f", offset);
    dir_turtle.pos = calc_point_on_bezier(offset, start_point, end_point);
    CHTurtle branch_pos_turtle = CHTurtle(dir_turtle);
    //logerr("dir %f %f %f right %f %f %f",branch_pos_turtle.dir.x, branch_pos_turtle.dir.y, branch_pos_turtle.dir.z,
    //branch_pos_turtle.right.x, branch_pos_turtle.right.y, branch_pos_turtle.right.z);
    branch_pos_turtle.pitch_down(90);
    branch_pos_turtle.move(radius_limit);
    //logerr("%f %f",dir_turtle.pos.x, branch_pos_turtle.pos.x);
    if (!std::isfinite(start_point.co.x))
    {
        throw std::exception();
    }
    if (!std::isfinite(end_point.co.x))
    {
        throw std::exception();
    }
    if (!std::isfinite(branch_pos_turtle.pos.x))
    {
        //logerr("dir %f %f %f right %f %f %f",dir_turtle.dir.x, dir_turtle.dir.y, dir_turtle.dir.z,
        //dir_turtle.right.x, dir_turtle.right.y, dir_turtle.right.z);
        //logerr("%d",(int)radius_limit/0);
        throw std::exception();
    }
    return branch_pos_turtle;
}

float WeberPennGenerator::Tree::calc_rotate_angle(int depth, float prev_angle)
{
    //calc rotate angle as defined in paper, limit to 0-360
    float r_angle = 0;
    if (param.rotate[depth] >= 0)
        r_angle = int(prev_angle + param.rotate[depth] + random_uniform(-1, 1) * param.rotate_v[depth]) % 360;
    else
        r_angle = prev_angle * (180 + param.rotate[depth] + random_uniform(-1, 1) * param.rotate_v[depth]);
    return r_angle;
}

float WeberPennGenerator::Tree::calc_down_angle(Stem &stem, float stem_offset)
{
    //calc down angle as defined in paper
    int d_plus_1 = MIN(stem.depth + 1, 3);
    float d_angle = 0;
        if (param.down_angle_v[d_plus_1] >= 0)
            d_angle = param.down_angle[d_plus_1] + random_uniform(-1, 1) * param.down_angle_v[d_plus_1];
        else
        {
            d_angle = param.down_angle[d_plus_1] + (param.down_angle_v[d_plus_1] * (
                1 - 2 * shape_ratio(0, (stem.length - stem_offset) / (stem.length * (
                    1 - param.base_size[stem.depth])))));
            //introduce some variance to improve visual result
            d_angle += random_uniform(-1, 1) * abs(d_angle * 0.1);
        }
        return d_angle;    
}

void WeberPennGenerator::Tree::clear()
{
    for (auto *s : stems)
    {
        if (s)
            delete s;
    }
    stems = {};
}

WeberPennGenerator::Stem::Stem(int _depth, int _spline_pos, Stem *_parent, float _offset, float _radius_limit)
{
    spline_pos = _spline_pos;
    depth = _depth;
    parent = _parent;
    offset = _offset;
    radius_limit = _radius_limit;
    children = {};
    length = 0;
    radius = 0;
    length_child_max = 0;

    if (parent)
    {
        parent->children.push_back(this);
    }
}

WeberPennGenerator::Stem::Stem(Stem &other)
{
    spline_pos = other.spline_pos;
    depth = other.depth;
    parent = other.parent;
    offset = other.offset;
    radius_limit = other.radius_limit;
    children = {};
    length = other.length;
    radius = other.radius;
    length_child_max = other.length_child_max;
    leaves = other.leaves;
    other.children.push_back(this);
    copied_from = &other;

}

void WeberPennGenerator::create_grove(GroveGenerationData ggd, ::Tree *trees_external, Heightmap &h)
{
    logerr("create_grove is not yet implemented for WeberPennGenerator");
}

void WeberPennGenerator::plant_tree(glm::vec3 pos, TreeTypeData *type)
{
    types.push_back(type);
    positions.push_back(pos);
}

void WeberPennGenerator::finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels)
{
    if (types.empty() || positions.size() != types.size())
    {
        return;
    }
    WeberPennParametersNative dummy_params;
    for (int i=0;i<types.size();i++)
    {
        vec3 pos = positions[i];
        WeberPennParametersNative *params = dynamic_cast<WeberPennParametersNative *>(types[i]->get_params());
        if (!params)
        {
            logerr("WeberPenn tree generator got wrong tree type");
            params = &dummy_params;
        }
        for (int j=0;j<params->levels;j++)
        {
            BranchHeap *br = new BranchHeap();
            trees_external[i].branchHeaps.push_back(br);
        }

        trees_external[i].leaves = new LeafHeap();
        trees_external[i].id = tree_next_id.fetch_add(1);
        trees_external[i].pos = pos;
        trees_external[i].type = types[i];
        trees_external[i].valid = true;

        Tree tree;
        tree.init(*params,true);
        tree.make();
        convert(tree, trees_external[i]);
        tree.clear();
    }
}

void WeberPennGenerator::convert(Tree &src, ::Tree &dst)
{
    if (!src.root)
        return;
    dst.root = dst.branchHeaps[0]->new_branch();
    dst.root->type_id = dst.type->type_id;
    dst.root->self_id = branch_next_id.fetch_add(1);
    dst.root->level = 0;
    dst.root->dead = false;
    dst.root->center_self = dst.pos;
    dst.root->center_par = vec3(0, 0, 0);
    dst.root->plane_coef = vec4(1, 0, 0, -dst.pos.x);
    dst.root->id = dst.id;
    
    convert(src, dst, src.root, dst.root);
    //logerr("convert %u", src.root, dst.branchHeaps[0]->branches.size());
    //logerr("convert %d leaves", src.leaves_array.size());
}

void WeberPennGenerator::convert(Tree &src, ::Tree &dst, Stem *src_br, ::Branch *dst_br)
{
    Spline &spline = src.branch_curves[src_br->depth].splines.data[src_br->spline_pos];
    float prev_r = 0;
    vec3 prev_pos = vec3(0,0,0);
    int same_cnt = 0;
    for (Point &p : spline.bezier_points)
    {
        vec3 pos = dst.pos + 10.0f*vec3(p.co.x, p.co.z, -p.co.y);
        if (length(pos - prev_pos) < 1e-4)
        {
            same_cnt++;
            pos.y += 1e-4*same_cnt;
        }
        else
        {
            same_cnt = 0;
        }
        //if (src_br->depth == 0)
        //    logerr("adding %f %f %f r= %f", pos.x, pos.y, pos.z, 10.0f*p.radius);
        float radius = 10.0f*p.radius;
        if (!dst_br->joints.empty())
        {
            dst_br->segments.emplace_back();
            auto &prev_j = dst_br->joints.back();
            auto &s = dst_br->segments.back();
            if (abs(radius - prev_r)>length(pos - prev_j.pos))
                radius = prev_r;
            s.rel_r_begin = prev_r;
            s.rel_r_end = radius;
            s.begin = prev_j.pos;
            s.end = pos;
        }
        dst_br->joints.emplace_back();
        auto &j = dst_br->joints.back();
        j.pos = pos;
        prev_r = radius;
        prev_pos = dst.pos + 10.0f*vec3(p.co.x, p.co.z, -p.co.y);
    }

    for (int &l_ind : src_br->leaves)
    {
        if (!dst_br->joints.back().leaf)
        {
            ::Leaf *l = dst.leaves->new_leaf();
            dst_br->joints.back().leaf = l;
        }
        std::vector<glm::vec3> out_verts;
        std::vector<std::vector<int>> out_indicies;

        BaseLeafMesh base_leaf_shape;
        Leaf::get_shape(9, src.tree_scale / src.param.g_scale,
                        src.param.leaf_scale, src.param.leaf_scale_x, base_leaf_shape);
        src.leaves_array[l_ind].get_mesh(src.param.leaf_bend, base_leaf_shape, 0, out_verts, out_indicies);

        for (auto &v : out_indicies)
        {
            ::Leaf *l = dst_br->joints.back().leaf;
            for (auto &ind : v)
            {
                vec3 pos = dst.pos + 10.0f*vec3(out_verts[ind].x, out_verts[ind].z, -out_verts[ind].y);
                l->edges.push_back(pos);
            }
            l->pos = l->edges[0];
            l->type = 0;
        }
    }

            for (Stem *s : src_br->children)
        {

          float best_dist = 1000;
          Joint *best_joint = nullptr;
          vec3 root_pos = vec3(0,0,0);
          int i=0;
          Spline &child_spline = src.branch_curves[s->depth].splines.data[s->spline_pos];
                      if (child_spline.bezier_points.empty())
            continue;
          if (!child_spline.bezier_points.empty())
          {
              //logerr("dist %f",length(child_spline.bezier_points[0].co - p.co));
          } 
          for (auto &j : dst_br->joints) 
          {
            glm::vec3 &init_p = spline.bezier_points[i].co;
            float len = length(init_p - child_spline.bezier_points[0].co);
            if (len < best_dist)
            {
                best_dist = len;
                best_joint = &j;
                root_pos = init_p;
            }
            i++;
          }
          //logerr("best dist %f curve (%d %d) child (%d %d)", best_dist, src_br->depth, src_br->spline_pos,
          //s->depth, s->spline_pos);
          if (best_dist < 1000)
          {
            float level = dst_br->level + 1;
            //level = src_br->depth + 1;
            if (dst.branchHeaps.size() == level)
            {
                BranchHeap *br = new BranchHeap();
                dst.branchHeaps.push_back(br);
            }
            child_spline.bezier_points.front().co = root_pos;
            s->already_used = true;
            Branch *br = dst.branchHeaps[level]->new_branch();
            br->type_id = dst.type->type_id;
            br->self_id = branch_next_id.fetch_add(1);
            br->level = level;
            br->dead = false;
            br->center_self = best_joint->pos;
            br->center_par = dst_br->center_self;
            br->plane_coef = vec4(1, 0, 0, -best_joint->pos.x);
            br->id = dst.id; 

            best_joint->childBranches.push_back(br);
            convert(src, dst, s, br);
          }
        }
}