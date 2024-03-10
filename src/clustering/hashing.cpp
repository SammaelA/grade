#include "hasing.h"
#include "dist_data_table.h"
#include <boost/math/constants/constants.hpp>
struct HashingParams
{
    int hash_count = 4;
    int EV_hasing_voxels_per_cell = 17;
    int EV_hasing_cells = 5;
    float max_individual_dist = 1000;
    bool voxelized_structure = true;
    float light_importance = 0.5;
    int ignore_structure_level = 2;
    float wood_size_mult = 1;
    void load_from_block(Block *b);
    void load(Block *b);
};
void HashingParams::load(Block *b)
{
    if (!b)
        return;
    
    Block &def = get_default_block();
    load_from_block(&def);
    load_from_block(b);
}
void HashingParams::load_from_block(Block *b)
{
    hash_count = b->get_int("hash_count",hash_count);
    EV_hasing_voxels_per_cell = b->get_int("EV_hasing_voxels_per_cell",EV_hasing_voxels_per_cell);
    EV_hasing_cells = b->get_int("EV_hasing_cells",EV_hasing_cells);
    max_individual_dist = b->get_double("max_individual_dist",max_individual_dist);
    voxelized_structure = b->get_bool("voxelized_structure",voxelized_structure);
    light_importance = b->get_double("light_importance",light_importance);
    ignore_structure_level = b->get_int("ignore_structure_level",ignore_structure_level);
    wood_size_mult = b->get_double("wood_size_mult",wood_size_mult);
}

HashingParams isParams;
void set_eigen_values_hash(LightVoxelsCube *voxels, Hash &hash, int cells, int voxels_per_cell, int sz);
BranchClusteringData *HashBasedClusteringHelper::convert_branch_eigin_vectors(Block &settings, Branch *base, 
                                                                              ClusteringContext *ctx, 
                                                                              BaseBranchClusteringData &data)
{
    isParams.load(&settings);
    BranchHash *branchHash = new BranchHash();

    BranchHeap branchHeap;
    LeafHeap leafHeap;
    Branch *b = branchHeap.new_branch();
    b->deep_copy(base, branchHeap, &leafHeap);
    auto tr = glm::inverse(data.transform); 
    b->transform(tr, data.r_transform);

    glm::vec3 axis = b->joints.back().pos - b->joints.front().pos;
    glm::mat4 rot = LiteMath::rotate(glm::mat4(1.0f), 2 * PI / isParams.hash_count, axis);

    for (int i = 0; i < isParams.hash_count; i++)
    {
        b->transform(rot);

            int sz_per_cell = isParams.EV_hasing_voxels_per_cell; 
            int cells = isParams.EV_hasing_cells;
            int sz = 0.5*(sz_per_cell*cells - 1);
            auto *l = new LightVoxelsCube(
                            glm::vec3(0.5f*canonical_bbox().x,0,0),
                            glm::ivec3(sz,sz,sz),
                            0.5f*canonical_bbox().x/sz);
            
            branchHash->hashes.emplace_back();
            set_occlusion(b, l);
            set_eigen_values_hash(l, branchHash->hashes.back(), cells, sz_per_cell, sz);

            if (isParams.voxelized_structure)
            {
                auto *vb = new LightVoxelsCube(
                                glm::vec3(0.5f*canonical_bbox().x,0,0),
                                glm::ivec3(sz,sz,sz),
                                0.5f*canonical_bbox().x/sz);
                
                branchHash->hashes.back().start_points.push_back(branchHash->hashes.back().data.size());
                branchHash->hashes.back().weights.emplace_back();
                branchHash->hashes.back().weights[0] = isParams.light_importance;
                branchHash->hashes.back().weights[1] = 1 - isParams.light_importance;
                voxelize_original_branch(b, vb, isParams.ignore_structure_level, isParams.wood_size_mult);
                set_eigen_values_hash(vb, branchHash->hashes.back(), cells, sz_per_cell, sz);

                delete vb;
            }

            delete l;
    }

    return branchHash;
}

void HashBasedClusteringHelper::clear_branch_data(BranchClusteringData *base, ClusteringContext *ctx)
{
    delete base;
}

IntermediateClusteringData *HashBasedClusteringHelper::prepare_intermediate_data_ddt(Block &settings, 
                                                                                     std::vector<BranchClusteringData *> branches,
                                                                                     ClusteringContext *ctx)
{
    isParams.load(&settings);
    IntermediateClusteringDataDDT *data = new IntermediateClusteringDataDDT();
    data->branches = branches;
    data->ddt.create(branches.size());
    data->elements_count = branches.size();
    std::vector<BranchHash *> real_branches;
    for (int i = 0; i < branches.size(); i++)
    {
        real_branches.push_back(dynamic_cast<BranchHash *>(branches[i]));
        if (!real_branches.back())
        {
            logerr("CPUSSClusteringHelper error: wrong type of BranchClusteringData");
            return data;
        }
    }
    for (int i = 0; i < branches.size(); i++)
    {
        for (int j = 0; j < branches.size(); j++)
        {
            Answer a;
            DistData d;

            if (i == j)
            {
                a.exact = true;
                a.from = 0;
                a.to = 0;
                d.dist = 0;
                d.rotation = 0;
            }
            else if (j < i)
            {
                auto p = data->ddt.get(j,i);
                a = p.first;
                d = p.second;
            }
            else
            {
                int rots = real_branches[j]->hashes.size();
                float min_dist = 1000;
                int min_rot = 0;
                for (int r=0;r<rots;r++)
                {
                    float dist = Hash::L2_dist(real_branches[i]->hashes.front(),real_branches[j]->hashes[r]);
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        min_rot = r;
                    }
                }
                a.from = min_dist;
                a.to = min_dist;
                a.exact = true;
                d.dist = min_dist;
                d.rotation = (2*PI*min_rot)/rots;
            }
            data->ddt.set(i,j,a,d);
        }
    }

    return data;
}

IntermediateClusteringData *HashBasedClusteringHelper::prepare_intermediate_data_simple(Block &settings, 
                                                                                        std::vector<BranchClusteringData *> branches,
                                                                                        ClusteringContext *ctx)
{
    IntermediateClusteringDataVectorsList *data = new IntermediateClusteringDataVectorsList();
    data->branches = branches;
    data->elements_count = branches.size();
    std::vector<BranchHash *> real_branches;
    for (int i = 0; i < branches.size(); i++)
    {
        real_branches.push_back(dynamic_cast<BranchHash *>(branches[i]));
        if (!real_branches.back())
        {
            logerr("CPUSSClusteringHelper error: wrong type of BranchClusteringData");
            return data;
        }
    }
    for (int i = 0; i < real_branches.size(); i++)
    {
        data->feature_vectors.emplace_back();
        real_branches[i]->hashes[0].weight();
        std::vector<GLfloat> &hash = real_branches[i]->hashes[0].data;
        for (GLfloat &f : hash)
        {
            data->feature_vectors.back().push_back(f);
        }
    }
    return data;
}

template<typename Real>
std::vector<Real> eigen_values(const Real A[3][3])
{
    using boost::math::constants::third;
    using boost::math::constants::pi;
    using boost::math::constants::half;

    static_assert(std::numeric_limits<Real>::is_iec559,
                  "Template argument must be a floating point type.\n");

    std::vector<Real> eigs(3, std::numeric_limits<Real>::quiet_NaN());
    auto p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    auto diag_sq = A[0][0]*A[0][0] + A[1][1]*A[1][1] + A[2][2];
    if (p1 == 0 || 2*p1/(2*p1 + diag_sq) < std::numeric_limits<Real>::epsilon())
    {
        eigs[0] = A[0][0];
        eigs[1] = A[1][1];
        eigs[2] = A[2][2];
        return eigs;
    }

    auto q = third<Real>()*(A[0][0] + A[1][1] + A[2][2]);
    auto p2 = (A[0][0] - q)*(A[0][0] - q) + (A[1][1] - q)*(A[1][1] -q) + (A[2][2] -q)*(A[2][2] -q) + 2*p1;
    auto p = std::sqrt(p2/6);
    auto invp = 1/p;
    Real B[3][3];
    B[0][0] = A[0][0] - q;
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][1] = A[1][1] - q;
    B[1][2] = A[1][2];
    B[2][2] = A[2][2] - q;
    auto detB = B[0][0]*(B[1][1]*B[2][2] - B[1][2]*B[1][2])
              - B[0][1]*(B[0][1]*B[2][2] - B[1][2]*B[0][2])
              + B[0][2]*(B[0][1]*B[1][2] - B[1][1]*B[0][1]);
    auto r = invp*invp*invp*half<Real>()*detB;
    if (r >= 1)
    {
        eigs[0] = q + 2*p;
        eigs[1] = q - p;
        eigs[2] = 3*q - eigs[1] - eigs[0];
        return eigs;
    }

    if (r <= -1)
    {
        eigs[0] = q + p;
        eigs[1] = q - 2*p;
        eigs[2] = 3*q - eigs[1] - eigs[0];
        return eigs;
    }

    auto phi = third<Real>()*std::acos(r);
    eigs[0] = q + 2*p*std::cos(phi);
    eigs[1] = q + 2*p*std::cos(phi + 2*third<Real>()*pi<Real>());
    eigs[2] = 3*q - eigs[0] - eigs[1];

    return eigs;
}

void set_eigen_values_hash(LightVoxelsCube *voxels, Hash &hash, int cells, int voxels_per_cell, int sz)
{
    for (int x_0 = -sz; x_0<sz; x_0+=voxels_per_cell)
    {
        for (int y_0 = -sz; y_0<sz; y_0+=voxels_per_cell)
        {
            for (int z_0 = -sz; z_0<sz; z_0+=voxels_per_cell)
            {
                double e_1 = 0, e_2 = 0, e_3 = 0;
                #define EPS 0.01
                double x_c = 0, y_c = 0, z_c = 0, sum = 0;
                for (int x = x_0; x<x_0+voxels_per_cell;x++)
                {
                    for (int y = y_0; y<x_0+voxels_per_cell;y++)
                    {
                        for (int z = z_0; z<z_0+voxels_per_cell;z++)
                        {
                            float val = voxels->get_occlusion_voxel_unsafe(x,y,z);
                            x_c += x*val;
                            y_c += y*val;
                            z_c += z*val;
                            sum += val;
                        }                    
                    }
                }
                if (sum > EPS)
                {
                    x_c /= sum;
                    y_c /= sum;
                    z_c /= sum;
                    double cov_mat[3][3];
                    for (int i=0;i<3;i++)
                    {
                        for (int j=0;j<3;j++)
                        {
                            cov_mat[i][j] = 0;
                        }
                    }
                    for (int x = x_0; x<x_0+voxels_per_cell;x++)
                    {
                        for (int y = y_0; y<x_0+voxels_per_cell;y++)
                        {
                            for (int z = z_0; z<z_0+voxels_per_cell;z++)
                            {
                                float val = voxels->get_occlusion_voxel_unsafe(x,y,z);
                                cov_mat[0][0] += val*(x-x_c)*(x-x_c);
                                cov_mat[0][1] += val*(x-x_c)*(y-y_c);
                                cov_mat[0][2] += val*(x-x_c)*(z-z_c);

                                cov_mat[1][0] += val*(x-x_c)*(y-y_c);
                                cov_mat[1][1] += val*(y-y_c)*(y-y_c);
                                cov_mat[1][2] += val*(y-y_c)*(z-z_c);

                                cov_mat[2][0] += val*(z-z_c)*(x-x_c);
                                cov_mat[2][1] += val*(y-y_c)*(z-z_c);
                                cov_mat[2][2] += val*(z-z_c)*(z-z_c);
                            }                    
                        }
                    }
                    for (int i=0;i<3;i++)
                    {
                        for (int j=0;j<3;j++)
                        {
                            cov_mat[i][j] /= sum;
                        }
                    }
                    std::vector<double> e_vals = eigen_values(cov_mat);
                    //logerr("eigen values %f %f %f",e_vals[0],e_vals[1],e_vals[2]);
                    e_1 = e_vals[0];
                    e_2 = e_vals[1];
                    e_3 = e_vals[2];
                }

                hash.data.push_back(e_1);
                hash.data.push_back(e_2);
                hash.data.push_back(e_3);
            }
        }
    }
}