#include "GPU_clusterization.h"
#include "tinyEngine/engine.h"

int min_count = 0;
double error_full = 0;
int missed = 0;
GPUClusterizationHelper::GPUClusterizationHelper():
distCompute({"dist_compute.comp"},{})
{

}
void GPUClusterizationHelper::prepare_ddt(std::vector<BranchWithData *> &_branches,
                                          DistDataTable &ddt, ClassicStructureSimilarityParams &cp)
{

    branches_size = _branches.size();
    debugl(13,"started GPU clusterization %d elements\n", branches_size);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    params.bwd_rotations = cp.bwd_rotations;
    params.delta = cp.delta;
    params.different_types_tolerance = cp.different_types_tolerance;
    params.ignore_structure_level = cp.ignore_structure_level;
    params.light_importance = cp.light_importance;
    params.max_individual_dist = cp.max_individual_dist;
    params.min_clusters = cp.min_clusters;
    params.voxels_size_mult = cp.voxels_size_mult;
    params.voxelized_structure = cp.voxelized_structure;
    params.structure_voxels_size_mult = cp.structure_voxels_size_mult;
    for (int i = 0; i<MAX_BRANCH_LEVELS; i++)
    {
        params.light_weights[i] = cp.light_weights[MIN(i,cp.light_weights.size()-1)];
        params.r_weights[i] = cp.r_weights[MIN(i,cp.r_weights.size()-1)];
        params.weights[i] = cp.weights[MIN(i,cp.weights.size()-1)];
    }
    float rt = (2*PI)/params.bwd_rotations;
    for (int i=0;i<params.bwd_rotations;i++)
    {
        float3 axis = float3(1,0,0);
        rotates_transforms.push_back(LiteMath::rotate(float4x4(),i*rt,axis));
    }

    positions.emplace_back();
    joint_rs.emplace_back();
    sticks.emplace_back();
    joints.emplace_back();
    cur_voxels_pointer = 1;
    cur_structure_voxels_pointer = 1;
    int all_voxels_count = 0;
    int all_structure_voxels_count = 0;
    bool voxels_needed = params.light_importance > IMPORTANCE_EPS;
    bool voxelized_structure = params.voxelized_structure && (1 - params.light_importance > IMPORTANCE_EPS);

    dist_data = safe_new<float>(2*branches_size*branches_size, "dist_data");
    if (voxels_needed)
    {
        int3 vsz = _branches.front()->leavesDensity.front()->get_vox_sizes();
        all_voxels_count = branches_size*_branches.front()->leavesDensity.size()*(2*vsz.x + 1)*(2*vsz.y + 1)*(2*vsz.z + 1);
        all_voxels = safe_new<float>(all_voxels_count + 1, "all_voxels");
    }
    if (voxelized_structure)
    {
        int3 vsz = _branches.front()->voxelizedStructures.front()->get_vox_sizes();
        all_structure_voxels_count = branches_size*_branches.front()->voxelizedStructures.size()*(2*vsz.x + 1)*(2*vsz.y + 1)*(2*vsz.z + 1);
        all_structure_voxels = safe_new<float>(all_structure_voxels_count + 1, "all_structure_voxels");
    }
    for (auto &bwd : _branches)
    {
        fill_branch_data(*bwd, voxels_needed, voxelized_structure);
    }
    int m_pos = positions.size(),m_sti = sticks.size(),m_j = joint_rs.size();
    int m_vox = all_voxels_count + 1,m_bd = branches.size(),m_dd = 2*branches_size*branches_size;
    int m_struct_vox = all_structure_voxels_count + 1;
    int bytes = m_pos*sizeof(float3) + m_sti*sizeof(gBranch) + m_j*sizeof(float) + m_vox*sizeof(float) +
    m_bd*sizeof(gBranchWithData) + m_dd*sizeof(float) + m_struct_vox*sizeof(float);

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    debugl(13,"Intermediate data allocated %4.2f Mbytes\n",bytes/powf(2,20));
    debugl(13,"preparation finished. Took %u [ms]\n",std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
    int hardness = 0;
    if ((_branches.front()->b->level == 0 && params.ignore_structure_level > 0) ||
       (params.ignore_structure_level > 2))
    {
        hardness = 1;
    }
    calculate_distances(hardness, false, false);

        for (int i=0;i<branches_size;i++)
        {
            for (int j=0;j<branches_size;j++)
            {    
                Answer a;
                DistData d;
                a.exact = true;
                a.from = dd_dist(i,j);
                if (a.from > 0.999*params.max_individual_dist)
                    a.from = 1e9;
                a.to = a.from;
                d.dist = a.from;
                d.rotation = dd_r(i,j);
                ddt.set(i,j,a,d);
            }
        }
}
void GPUClusterizationHelper::fill_branch_data(BranchWithData &branch, bool voxels_needed, bool voxelized_structure)
{
    branches.emplace_back();
    auto &br = branches.back();
    br.branch_id = branch.id;
    br.dead = false;
    br.type_id = branch.original->type_id;
    br.level = branch.original->level;

    for (int i = 0; i<MAX_BRANCH_LEVELS; i++)
    {
        br.joints_counts[i] = branch.joint_counts[MIN(i,branch.joint_counts.size()-1)];
    }

    if (voxelized_structure)
    {
        int3 sizes;
        float *ptr = nullptr;
        br.structure_voxels_offset = cur_structure_voxels_pointer;
        branch.voxelizedStructures.front()->get_data(&ptr,sizes);
        br.voxels_xyz_rots.w = branch.voxelizedStructures.size();
        br.structure_voxels_xyz_and_size = uvec4(sizes.x,sizes.y,sizes.z,
                                           (2*sizes.x + 1)*(2*sizes.y + 1)*(2*sizes.z + 1)*branch.voxelizedStructures.size());

        for (int k=0;k<branch.voxelizedStructures.size();k++)
        {
            branch.voxelizedStructures[k]->get_data(&ptr,sizes);
            int sz = br.structure_voxels_xyz_and_size.w/br.voxels_xyz_rots.w;
            memcpy(&(all_structure_voxels[br.structure_voxels_offset + k*sz]),ptr,sizeof(float)*sz);
        }
        
        cur_structure_voxels_pointer += br.structure_voxels_xyz_and_size.w;
    }
    else if (1 - params.light_importance > IMPORTANCE_EPS)
    {
        br.branch_pos = fill_branch_data(branch.b);
        calc_cumulative_weight(br.branch_pos,branch.b->level);
    }

    if (voxels_needed)
    {
        int3 sizes;
        float *ptr = nullptr;
        br.voxels_offset = cur_voxels_pointer;
        branch.leavesDensity.front()->get_data(&ptr,sizes);
        br.voxels_xyz_rots = uvec4(sizes.x,sizes.y,sizes.z,branch.leavesDensity.size());
        br.voxels_size = (2*sizes.x + 1)*(2*sizes.y + 1)*(2*sizes.z + 1)*branch.leavesDensity.size();

        for (int k=0;k<branch.leavesDensity.size();k++)
        {
            branch.leavesDensity[k]->get_data(&ptr,sizes);
            int sz = br.voxels_size/br.voxels_xyz_rots.w;
            memcpy(&(all_voxels[br.voxels_offset + k*sz]),ptr,sizeof(float)*sz);
        }
        
        cur_voxels_pointer += br.voxels_size;
    }
}
float GPUClusterizationHelper::calc_cumulative_weight(uint stick_id, int level)
{
    float sd = 0;
    float kd = 0;
    for (int i=0;i<sticks[stick_id].joint_count;i++)
    {
        #define jj(i) (joints[sticks[stick_id].joint_offset + i])
            sd += params.weights[level];
            for (int k=0;k<MAX_CHILD_BRANCHES;k++)
            {
                if (jj(i).child_branches_ids[k])
                {
                    kd += calc_cumulative_weight(jj(i).child_branches_ids[k],level+1);
                }
            }
    }
    kd += sd;
    sticks[stick_id].self_weight = sd;
    sticks[stick_id].cumulative_weight = kd;

    return kd;
}
uint GPUClusterizationHelper::fill_branch_data(Branch *branch)
{
    uint pos = sticks.size();
    sticks.emplace_back();

    int i = 0;
    auto s_mult = branch->segments.begin();
    sticks[pos].joint_offset = joints.size();
    sticks[pos].joint_count = branch->joints.size();
    for (Joint &j : branch->joints)
    {
        joints.emplace_back();
    }
    for (Joint &j : branch->joints)
    {
        if (i >= MAX_JOINTS)
        {
            logerr("GPU clusterization error: branch has too many joints %d",branch->joints.size());
            break;
        }  
        joints[sticks[pos].joint_offset + i].pos_id = positions.size();

        joints[sticks[pos].joint_offset + i].r_id = joint_rs.size();
        for (int k = 0;k<params.bwd_rotations;k++)
        {
            positions.push_back(float4(rotates_transforms[k]*to_float4(j.pos,1.0f)));
        }
        for (int k = 0;k<params.r_samples;k++)
        {
            float phi = ((float)k)*2*PI/params.r_samples;
            joint_rs.push_back(Branch::get_r_mult(phi,s_mult->mults));
        }
        if (branch->level < params.ignore_structure_level)
        {
            if (j.childBranches.size() > MAX_CHILD_BRANCHES)
            {
                logerr("GPU clusterization error: joint has too many child branches");
            }
            auto chbit = j.childBranches.begin();
            for (int k = 0;k<MIN(j.childBranches.size(),MAX_CHILD_BRANCHES);k++)
            {
                joints[sticks[pos].joint_offset + i].child_branches_ids[k] = fill_branch_data(*chbit);
                chbit++;
            }
        }
        i++;
        if (std::next(s_mult) != branch->segments.end())
        {
            s_mult++;
        }


    }
    return pos;
}

GPUClusterizationHelper::~GPUClusterizationHelper()
{
    safe_delete<float>(all_voxels, "all_voxels");
    safe_delete<float>(all_structure_voxels, "all_structure_voxels");
    safe_delete<float>(dist_data, "dist_data");

    #define DELBUF(a) if (a) { delete_buffer((a)); a = 0;}

    DELBUF(pos_buf);
    DELBUF(voxels_buf);
    DELBUF(structure_voxels_buf);
    DELBUF(rs_buf);
    DELBUF(sticks_buf);
    DELBUF(branches_buf);
    DELBUF(dist_data_buf);
    DELBUF(params_buf);
    DELBUF(joints_buf);
}
    void GPUClusterizationHelper::calculate_distances(int hardness, bool cpu_only, bool cpu_check)
    {   
        std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
        if (!cpu_only)
        {
            setup_buffers();
            distCompute.use();
            distCompute.uniform("branches_size",branches_size);
            static int max_dispatches = hardness > 0 ? 8 : 16;
            static int threads = 8;
            static int branches_per_thread = 1;
            int step = (max_dispatches*threads*branches_per_thread);
            int iters = ceil((float)branches_size/step);
            int a = 0;
            for (int i=0;i<iters;i++)
            {
                for (int j=0;j<iters;j++)
                {
                    a++;
                    int x_sz = (i == iters - 1) ? branches_size - i*step : step;
                    int y_sz = (j == iters - 1) ? branches_size - j*step : step;
                    int start_x = i*step;
                    int start_y = j*step;
                    
                    if (start_y + y_sz < x_sz)
                    continue;
                    distCompute.uniform("start_x",start_x);
                    distCompute.uniform("start_y",start_y);
                    debugl(13,"GPU clusterization in process %d/%d\n",a,iters*iters);
                    glDispatchCompute(ceil((float)x_sz/threads),ceil((float)y_sz/threads), 1);

                    engine::view->next_frame();
                }
            }
            
            glMemoryBarrier(GL_ALL_BARRIER_BITS);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, dist_data_buf);
            GLvoid* p = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
            memcpy(dist_data,p,sizeof(float)*2*branches_size*branches_size);
        }
        std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        if (cpu_only)
        {
            for (int i=0;i<branches_size;i++)
            {
                for (int j=0;j<branches_size;j++)
                {
                    if (j > i)
                        calculate_dist(i,j);
                }
            }
        }
        else if (cpu_check)
        {
            for (int i=0;i<branches_size;i++)
            {
                for (int j=0;j<branches_size;j++)
                {
                    if (j > i)
                    {
                        float d = dd_dist(i,j);
                        float r = dd_r(i,j);
                        calculate_dist(i,j);
                        if (abs(d - dd_dist(i,j)) > 1e-4)
                        logerr("%d %d from %d compute shader error %f %f --> %f %f",i,j,branches_size,d,r,dd_dist(i,j),dd_r(i,j));
                    }
                }
            }
        }
        for (int i=0;i<branches_size;i++)
        {
            for (int j=0;j<branches_size;j++)
            {
                if (i == j)
                {
                    dd_r(i,i) = 0;
                    dd_dist(i,i) = 0;
                }
                else if (j < i)
                {
                    dd_r(i,j) = dd_r(j,i);
                    dd_dist(i,j) = dd_dist(j,i);
                }
            }
        }

        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        debugl(17,"GPU distance calculation finished. Took %u [ms]\n",std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count());
        debugl(1,"CPU distance calculation finished. Took %u [ms]\n",std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
        
        positions.clear();
        safe_delete<float>(all_voxels, "all_voxels");
        all_voxels = nullptr;
        joint_rs.clear();
        sticks.clear();
        branches.clear();
        joints.clear();
    }
    #define uvec3 uint3
    #define vec3 float3
    #define vec2 float2
    
    float GPUClusterizationHelper::calculate_dist_structure(int i, int j, int rot, float max_dist)
    {
        float e = -11;
        #define STACK_SIZE (MAX_CHILD_BRANCHES*MAX_JOINTS*MAX_JOINTS)
        uvec3 ring_stack[STACK_SIZE];//id1,id2,level
        uint joints_passed[MAX_BRANCH_LEVELS];
        float matches[MAX_BRANCH_LEVELS];
        bool pair_found[MAX_JOINTS];
        int bottom = 0;
        int top = 1;
        float kd = 0;
        float fd = 0;
        bool more_than_max = false;
        int rot1 = 0;
        int rot2 = rot;
        ring_stack[bottom].x = branches[i].branch_pos;
        ring_stack[bottom].y = branches[j].branch_pos;
        ring_stack[bottom].z = branches[i].level;
        for (int k=0;k<MAX_BRANCH_LEVELS;k++)
        {
            joints_passed[k] = 0;
            matches[k] = 0;
        }
        for (int k=0;k<MAX_BRANCH_LEVELS;k++)
        {
            if (k > params.ignore_structure_level)
                break;
            fd += params.weights[k]*(branches[i].joints_counts[k] + branches[j].joints_counts[k]);
        }
        while (top != bottom && !more_than_max)
        {
            for (int j1 = 0;j1<MAX_JOINTS;j1++)
            {
                pair_found[j1] = false;
            }
            int bid_1 = ring_stack[bottom].x;
            int bid_2 = ring_stack[bottom].y;

            #define POS1(j1) positions[joints[sticks[bid_1].joint_offset + j1].pos_id + rot1]
            #define POS2(j2) positions[joints[sticks[bid_2].joint_offset + j2].pos_id + rot2]
            #define JJ1(j1)  joints[sticks[bid_1].joint_offset + j1]
            #define JJ2(j2)  joints[sticks[bid_2].joint_offset + j2]
            int b_level = ring_stack[bottom].z;
            int j1_sz = sticks[bid_1].joint_count;
            int j2_sz = sticks[bid_2].joint_count;

            float not_matched = j1_sz + j2_sz;

            if (j1_sz < j2_sz)
            {
                int t = j1_sz;
                j1_sz = j2_sz;
                j2_sz = t;
                bid_1 = ring_stack[bottom].y;
                bid_2 = ring_stack[bottom].x;
                rot1 = rot;
                rot2 = 0;
            }
            if (j1_sz <= 1)
            {
                if (j2_sz == j1_sz)
                {
                    matches[b_level] += j1_sz + j2_sz;
                    not_matched = 0;
                }
                
            }
            else
            {
                float l1 = length(POS1(j1_sz - 1) - POS1(0));
                float l2 = length(POS2(j2_sz - 1) - POS2(0));
                float delta = params.delta*(l1+l2)/2;

                for (int j1 = 0; j1 < j1_sz; j1++)
                {
                    float min_dist = delta;
                    int min_pos = -1;
                    vec3 pos1 = to_float3(POS1(j1));
                    for (int j2 = 0; j2 < j2_sz; j2++)
                    {
                        if (pair_found[j2])
                            continue;
                        float l = length(pos1 - to_float3(POS2(j2)));
                        if (l < min_dist)
                        {
                            min_dist = l;
                            min_pos = j2;
                        }
                    }
                    if (min_pos >= 0)
                    {
                        pair_found[min_pos] = true;
                        float r_diff = 0;
                        if (JJ1(j1).r_id && JJ2(min_pos).r_id && params.r_samples)
                        {
                            for (int t = 0;t<params.r_samples;t++)
                            {
                                float r1 = joint_rs[JJ1(j1).r_id + t];
                                float r2 = joint_rs[JJ2(min_pos).r_id + t];
                                r_diff += sqrt(abs(r1-r2)/(r1+r2+1e-6));
                            }
                            r_diff /= params.r_samples;
                        }
                        float match = 2*(1-r_diff);
                        matches[b_level] += match;
                        not_matched -= match;
                        if (b_level < params.ignore_structure_level)
                        {
                            int chbid1 = JJ1(j1).child_branches_ids[0];
                            int chbid2 = JJ2(min_pos).child_branches_ids[0];
                            if (chbid1 && chbid2)
                            {
                                ring_stack[top].x = chbid1;
                                ring_stack[top].y = chbid2;
                                ring_stack[top].z = b_level + 1;
                                top = (top + 1) % STACK_SIZE;
                            }
                            else if (chbid1)
                            {
                                kd += sticks[chbid1].cumulative_weight;
                            }
                            else if (chbid2)
                            {
                                kd += sticks[chbid2].cumulative_weight;
                            }
                        }
                    }
                    else
                    {
                        int chbid1 = JJ1(j1).child_branches_ids[0];
                        if (chbid1)
                        {
                            kd += sticks[chbid1].cumulative_weight;
                        }
                    }
                }
                for (int j2 = 0; j2 < j2_sz; j2++)
                {
                    if (!pair_found[j2])
                    {
                        int chbid2 = JJ2(j2).child_branches_ids[0];
                        if (chbid2)
                        {
                            kd += sticks[chbid2].cumulative_weight;
                        }
                    }
                }
            }
            kd += params.weights[b_level]*not_matched;
            if (kd/fd > max_dist)
            {
                if (e < 0)
                    e = kd/fd;
                more_than_max = true;
            }
            bottom = (bottom + 1) % STACK_SIZE;
        }
        if (more_than_max)
        {
            return 1;
        }
        else
        {
            
            float r = 0;
            float d = 0;
            for (int k=0;k<MAX_BRANCH_LEVELS;k++)
            {
                if (k > params.ignore_structure_level)
                    break;
                int jk = (branches[i].joints_counts[k] + branches[j].joints_counts[k]);
                if (jk > 0)
                {
                    r += params.weights[k]*matches[k];
                    d += params.weights[k]*jk;
                }
            }

            return kd/fd;
        }
    }
    float GPUClusterizationHelper::calculate_dist_light(int i, int j, int rot, float max_dist)
    {
        if (branches[i].voxels_size != branches[j].voxels_size || 
            branches[i].voxels_xyz_rots.w != branches[j].voxels_xyz_rots.w)
            return 1;
        int vs = branches[i].voxels_size/branches[i].voxels_xyz_rots.w;
        float a = 0;
        float b = 0.01;
        int off1 = branches[i].voxels_offset;
        int off2 = branches[j].voxels_offset + rot*vs;

        for (int k = 0; k<vs;k++)
        {
            a += abs(all_voxels[off1 + k] - all_voxels[off2 + k]);
            b += (all_voxels[off1 + k] + all_voxels[off2 + k]);
        }
        return a/b;
    }
    float GPUClusterizationHelper::calculate_dist_voxelized_structure(int i, int j, int rot, float max_dist)
    {
        if (branches[i].structure_voxels_xyz_and_size.w != branches[j].structure_voxels_xyz_and_size.w || 
            branches[i].voxels_xyz_rots.w != branches[j].voxels_xyz_rots.w)
            return 1;
        int vs = branches[i].structure_voxels_xyz_and_size.w/branches[i].voxels_xyz_rots.w;
        float a = 0;
        float b = 0.01;
        int off1 = branches[i].structure_voxels_offset;
        int off2 = branches[j].structure_voxels_offset + rot*vs;

        for (int k = 0; k<vs;k++)
        {
            a += abs(all_structure_voxels[off1 + k] - all_structure_voxels[off2 + k]);
            b += (all_structure_voxels[off1 + k] + all_structure_voxels[off2 + k]);
        }
        //logerr("vox dist %f",a/b);
        return a/b;
    }
    vec2 GPUClusterizationHelper::calculate_dist_simple(int i, int j, int rot, float min_dist_struct)
    {
        //empirically got that if structure distance is > 1.1*min_dist_struct then
        //the full distance is probably more that min_dist
        #define MAX_STR_DIST_MULT 1.1
        float MS = 0, ML = 0, a = params.light_importance;
        if (a > IMPORTANCE_EPS)
        {
            ML = calculate_dist_light(i,j,rot,1);
        }
        if (1 - a > IMPORTANCE_EPS)
        {
            if (params.voxelized_structure)
            {
                MS = calculate_dist_voxelized_structure(i,j,rot,1);
            }
            else
            {
                MS = calculate_dist_structure(i,j,rot,MAX_STR_DIST_MULT*min_dist_struct);
            }
        }
        if (a*ML + (1-a)*MS < 0)
        {
            logerr("ERROR %f %f %f",a,ML,MS);
        }
        //if (branches[i].level == 0)
        //    logerr("dist_simple[%d %d] = %f %f %f",i,j,a, ML, MS);
        return vec2(a*ML + (1-a)*MS,MS);
    }
    void GPUClusterizationHelper::calculate_dist(int _i, int _j)
    {
        int i = _i;
        int j = _j;

        if ( i<0 || i>=branches_size || j<0 || j>=branches_size)
            return;
        
        if  ((branches[i].type_id != branches[j].type_id && !params.different_types_tolerance)|| 
             branches[i].dead != branches[i].dead || 
             branches[i].level != branches[i].level)
        {
            dd_r(i,j) = 0;
            dd_dist(i,j) = 1000;
            return;
        }

        float min_dist = params.max_individual_dist;
        float min_structure_dist = 1;
        int rot = 0;
        for (int k =0;k<params.bwd_rotations;k++)
        {
            vec2 dist = calculate_dist_simple(i,j,k,min_structure_dist);
            if (dist.x<min_dist)
            {
                min_dist = dist.x;
                min_structure_dist = dist.y;
                rot = k;
            }
        }

        dd_r(i,j) = (2*PI*rot)/params.bwd_rotations;
        dd_dist(i,j) = min_dist;
        //if (branches[i].level == 0)
        //    logerr("dist[%d %d] = %f",i,j,min_dist);
    }
    void GPUClusterizationHelper::setup_buffers()
    {
        pos_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, pos_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float4)*positions.size(), positions.data(), GL_STATIC_DRAW);
    
        voxels_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, voxels_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*cur_voxels_pointer, all_voxels, GL_STATIC_DRAW);
            
        structure_voxels_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, structure_voxels_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*cur_structure_voxels_pointer, all_structure_voxels, GL_STATIC_DRAW);

        rs_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, rs_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*joint_rs.size(), joint_rs.data(), GL_STATIC_DRAW);

        sticks_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, sticks_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(gBranch)*sticks.size(), sticks.data(), GL_STATIC_DRAW);
        
        branches_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, branches_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(gBranchWithData)*branches.size(), branches.data(), GL_STATIC_DRAW);

        dist_data_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, dist_data_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*2*branches_size*branches_size, dist_data, GL_STATIC_DRAW);

        params_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, params_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(gClusterizationParams), &params, GL_STATIC_DRAW);

        joints_buf = create_buffer();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, joints_buf);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(gJoint)*joints.size(), joints.data(), GL_STATIC_DRAW);


    }