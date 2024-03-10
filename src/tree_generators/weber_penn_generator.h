#pragma once
#include "abstract_generator.h"
#include "common_utils/quaternion.h"

class CHTurtle;
struct SetUpBranchRetStruct;
struct WeberPennParametersNative : public ParameterSet
{
    int shape = 7;
    float g_scale = 13;
    float g_scale_v = 3;
    int levels = 3;
    float ratio = 0.015;
    float ratio_power = 1.2;
    float flare = 0.6;
    int base_splits = 0;
    float4 base_size= float4(0.15, 0.02, 0.02, 0.02);
    float4 down_angle= float4(-0, 60, 45, 45);
    float4 down_angle_v= float4(-0, -50, 10, 10);
    float4 rotate= float4(-0, 140, 140, 77);
    float4 rotate_v= float4(-0, 0, 0, 0);
    float4 branches= float4(1, 50, 60, 10);
    float4 length= float4(1, 0.3, 0.6, 0);
    float4 length_v= float4(0, 0, 0, 0);
    float4 taper= float4(1, 1, 1, 1);
    float4 seg_splits= float4(0, 0, 0, 0);
    float4 split_angle= float4(0, 0, 0, 0);
    float4 split_angle_v= float4(0, 0, 0, 0);
    float4 curve_res= float4(15, 5, 3, 1);
    float4 curve= float4(0, -40, -40, 0);
    float4 curve_back= float4(0, 0, 0, 0);
    float4 curve_v= float4(20, 50, 75, 0);
    float4 bend_v= float4(-0, 50, 0, 0);
    float4 branch_dist= float4(-0, 0, 0, 0);
    float4 radius_mod = float4(1,1,1,1); 
    int leaf_blos_num = 40;
    int leaf_shape = 10;
    float leaf_scale = 0.17;
    float leaf_scale_x = 1;
    float leaf_bend = 0.6;
    int blossom_shape = 1;
    float blossom_scale = 0;
    float blossom_rate = 0;
    float leaf_rate = 1;
    float3 tropism = float3(0,0,0.5);
    float prune_ratio = 0;
    float prune_width = 0.5;
    float prune_width_peak = 0.5;
    float prune_power_high = 0.5; 
    float prune_power_low = 0.5;

    virtual float3 get_tree_max_size() override
    {
        return float3(100,200,100);
    }
    virtual float get_scale_factor() override 
    {
        return 1;
    }
    virtual ParameterSet *copy() override
    { 
        auto Ps = new WeberPennParametersNative();
        *Ps = *this;
        return Ps;
    };
    
    virtual void save_load_define(SaveLoadMode mode, Block &b, ParameterList &list) override;
};

class WeberPennGenerator : public AbstractTreeGenerator
{
public:
    virtual bool iterate(LightVoxelsCube &voxels) override { return false;};
    virtual void plant_tree(float3 pos, const TreeTypeData *type) override;
    virtual void finalize_generation(::Tree *trees_external, LightVoxelsCube &voxels) override;
    virtual void set_seed(int _seed) override {seed = _seed;}
//private:
    struct Tree;
    struct Stem;
    void convert(Tree &src, ::Tree &dst);
    void convert(Tree &src, ::Tree &dst, Stem *src_br, ::Branch *dst_br);
    #define MAX_DEPTH 8
    struct BaseLeafMesh
    {
        std::vector<float3> verts;
        std::vector<std::vector<int>> indicies;
    };
    struct Leaf
    {
        float3 position;
        float3 direction;
        float3 right;
        Leaf(float3 pos,float3 dir, float3 r)
        {
            position = pos;
            direction = dir;
            right = r;
        }
        static void get_shape(int leaf_type, float g_scale, float scale, float scale_x, BaseLeafMesh &blm);
        void get_mesh(float bend, BaseLeafMesh &base_shape, int index, std::vector<float3> &out_verts, 
                      std::vector<std::vector<int>> &out_indicies);
        void calc_bend_trf(float bend, LiteMath::quat &bend_trf_1, LiteMath::quat &bend_trf_2);
    };
    struct Point
    {
        float3 co = float3(0,0,0);
        float3 handle_left = float3(0,0,0);
        float3 handle_right = float3(0,0,0);
        float radius = 1.0;
    };
    struct Spline
    {
        std::vector<Point> bezier_points = {Point()};
        float resolution_u = 1.0;
        void add() {bezier_points.push_back(Point());}
    };
    struct Splines
    {
        std::vector<Spline> data;
    };
    struct Curve
    {
        float resolution_u = 1;
        Splines splines = Splines();
    };
    struct Stem
    {
        //Spline &curve;
        bool already_used = false;
        int spline_pos;
        int depth;
        Stem *parent = nullptr;
        Stem *copied_from = nullptr;
        float offset;
        float radius_limit;
        std::vector<int> leaves;
        std::vector<Stem *> children;
        float length;
        float radius;
        float length_child_max;
        Stem(int depth, int _spline_pos, Stem *_parent = nullptr, float _offset = 0, float _radius_limit=-1);
        Stem& operator=(const Stem& other) = delete;
        Stem& operator=(Stem&& other)  = delete;
        Stem(Stem &other);//do not copy children!!
        //Stem(Stem &&other) = default;
        //Stem& operator=(const Stem& other)
        //{ return *this = Stem(other);}
        //Stem& operator=(Stem&& other) = default;
    };

    enum BranchMode
    {
        alt_opp = 1,
        whorled = 2,
        fan = 3
    };
    struct Tree
    {
        WeberPennParametersNative param = WeberPennParametersNative();

        void init(WeberPennParametersNative _param, bool _generate_leaves);
        void make();    
        void create_branches();
        void create_leaf_mesh();
        void points_for_floor_split(std::vector<std::pair<float3, float>> &points);
        float random_uniform(float from, float to);
        unsigned long random_getstate() {return seed;}
        void random_setstate(unsigned long _seed){ seed = _seed;}
        void make_stem(CHTurtle &turtle, Stem &stem, int start=0, float split_corr_angle=0, float num_branches_factor=1,
                       float clone_prob=1, CHTurtle *pos_corr_turtle=nullptr, CHTurtle *cloned_turtle=nullptr);
        float calc_stem_length(Stem &stem);
        float calc_stem_radius(Stem &stem);
        bool test_stem(CHTurtle &turtle, Stem &stem, int start=0, float split_corr_angle=0, float clone_prob=1);
        int calc_leaf_count(Stem &stem);
        float calc_branch_count(Stem &stem);
        void apply_tropism(CHTurtle &turtle, float3 tropism_vector);
        void calc_helix_points(CHTurtle &turtle, float hel_radius, float hel_pitch, 
                               float3 &hel_p_0, float3 &hel_p_1, float3 &hel_p_2, float3 &hel_axis);
        float radius_at_offset(Stem &stem, float z_1);
        void make_branches(CHTurtle &turtle, Stem &stem, int seg_ind, int branches_on_seg, 
                           std::vector<float> &prev_rotation_angle, bool is_leaves=false);
        void make_leaves(CHTurtle &turtle, Stem &stem, int seg_ind, int branches_on_seg, 
                         std::vector<float> &prev_rotation_angle);
        float calc_curve_angle(int depth, int seg_ind);
        void make_clones(CHTurtle &turtle, int seg_ind, float split_corr_angle,float num_branches_factor, 
                         float clone_prob, Stem &stem, int num_of_splits, float spl_angle, float spr_angle, 
                         bool is_base_split);
        void increase_bezier_point_res(Stem &stem, int seg_ind, int points_per_seg);
        void scale_bezier_handles_for_flare(Stem &stem, int points_per_seg);
        float shape_ratio(int shape, float ratio);
        bool point_inside(float3 point);
        SetUpBranchRetStruct set_up_branch(CHTurtle &turtle, Stem &stem, BranchMode branch_mode, float offset, 
                                           Point &start_point, Point &end_point, float stem_offset, int branch_ind, 
                                           std::vector<float> &prev_rot_ang, int branches_in_group=0);
        float3 calc_point_on_bezier(float offset, Point &start_point, Point &end_point);
        float3 calc_tangent_to_bezier(float offset, Point &start_point, Point &end_point);
        CHTurtle make_branch_dir_turtle(CHTurtle &turtle, bool helix, float offset, Point &start_point, Point &end_point);
        CHTurtle make_branch_pos_turtle(CHTurtle &turtle, float offset, Point &start_point, Point &end_point,
                                        float radius_limit);
        float calc_rotate_angle(int depth, float prev_angle);
        float calc_down_angle(Stem &stem, float stem_offset);
        Spline &spline(Stem &stem) {return branch_curves[stem.depth].splines.data[stem.spline_pos];}
        void clear();
        bool generate_leaves = true;
        std::vector<Leaf> leaves_array = {}; 
        int stem_index = 0;
        float tree_scale = 0;
        std::vector<Curve> branch_curves = {}; 
        float base_length = 0;
        std::vector<float> split_num_error = std::vector<float>(MAX_DEPTH,0);
        float trunk_length = 0;
        std::vector<Stem *> stems;
        Stem *root = nullptr;
        unsigned long seed = 11;
        int total_points_cnt = 0;
    };

    std::vector<const TreeTypeData *> types;
    std::vector<float3> positions;
    static WeberPennParametersNative defaultParameters;
    unsigned long seed = 0;

};
float declination(float3 vec);
LiteMath::quat to_track_quat_ZY(float3 vec);