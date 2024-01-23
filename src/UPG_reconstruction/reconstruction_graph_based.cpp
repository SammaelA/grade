#include "reconstruction_graph_based.h"
#include "sdf_rendering.h"
#include "common_utils/bbox.h"
#include "common_utils/distribution.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"

namespace upg
{
  struct GRPrimitive
  {
    UPGStructure structure;
    UPGParametersRaw parameters;
  };

  struct GREdge;
  struct GRNode
  {
    GRNode(const FieldSdfCompare &_opt_func):
      opt_func(_opt_func)
    {

    }
    FieldSdfCompare opt_func;
    std::vector<GRPrimitive> primitives;
    std::vector<GREdge> edges;
    GRNode *parent = nullptr;

    float quality = -1;
    int depth = 0;
  };

  struct GREdge
  {
    GRPrimitive start_prim;//starting parameters
    GRPrimitive opt_prim;//parameters after regular optimization
    GRPrimitive fine_prim;//parameters after fine-tuning
    float opt_quality = -1;
    float fine_quality = -1;
    GRNode *start = nullptr;
    std::shared_ptr<GRNode> end;
  };

  struct GROptimizationSettings
  {
    float distance_base_thr = 0.01f;
    float distance_fine_thr = 0.002f;
  };

  struct GROptimizationContext
  {
    GROptimizationSettings settings;
    Block opt_settings_blk;
    Block fine_opt_settings_blk;
    
    std::shared_ptr<GRNode> root;
    std::vector<glm::vec3> target_points;
    std::vector<float> target_distances;
    UPGStructure target_structure; //can be empty, it means we should find proper structure during optimization
    std::vector<UPGPart> target_structure_parts; //can be empty
    AABB point_cloud_bbox;
  };

  UPGStructure get_structure_known_structure(const GROptimizationContext &ctx, const GRNode &node)
  {
    UPGStructure part_structure;
    auto &part = ctx.target_structure_parts[node.depth];
    for (int j=part.s_range.first; j<part.s_range.second; j++)
      part_structure.s.push_back(ctx.target_structure.s[j]);
    logerr("depth %d created structure %d %d", node.depth, (int)part_structure.s[0], (int)part_structure.s[1]);
    return part_structure;
  }
  UPGStructure get_random_structure(const GROptimizationContext &ctx, const GRNode &node)
  {
    logerr("get_random_structure is not implemented!");
    return UPGStructure();
  }

  glm::vec3 get_position_min_distance(const GRNode &node, int tries = 25)
  {
    glm::vec3 best_pos;
    float best_dist = 1e6;
    for (int k = 0; k < tries; k++)
    {
      unsigned index = rand() % node.opt_func.a_points.size();
      if (node.opt_func.a_distances[index] < best_dist)
      {
        best_dist = node.opt_func.a_distances[index];
        best_pos = node.opt_func.a_points[index];
      }
    }
    return best_pos;
  }

  UPGParametersRaw get_initial_parameters_random(const GRNode &node, const UPGStructure &structure, glm::vec3 initial_pos)
  {
    auto gen = node.opt_func.get_generator(structure);
    auto pd = node.opt_func.get_full_parameters_description(gen.get());

    glm::ivec3 pos_ids(-1,-1,-1);
    std::vector<glm::vec2> borders;
    for (const auto &p : pd.get_block_params())
    {
      for (const auto &param_info : p.second.p)
      {
        if (param_info.type != ParameterType::CONST)
          borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
        if (param_info.name=="move_x")
        {
          assert(pos_ids.x == -1);
          pos_ids.x = borders.size()-1;
        }
        else if (param_info.name=="move_y")
        {
          assert(pos_ids.y == -1);
          pos_ids.y = borders.size()-1;
        }
        else if (param_info.name=="move_z")
        {
          assert(pos_ids.z == -1);
          pos_ids.z = borders.size()-1;
        }
      }
    }
    assert(pos_ids.x>=0 && pos_ids.y>=0 && pos_ids.z>=0);

    UPGParametersRaw res;
    res.p.resize(borders.size());
    float range = 0.5;

    for (int i=0;i<borders.size();i++)
    {
      float sz = borders[i].y - borders[i].x;
      if (borders[i].x < 0)
        res.p[i] = urand(borders[i].x + 0.5*(1-range)*sz, borders[i].y - 0.5*(1-range)*sz);
      else
        res.p[i] = urand(borders[i].x, borders[i].x + 0.25*sz);
    }

    res.p[pos_ids.x] = initial_pos.x;
    res.p[pos_ids.y] = initial_pos.y;
    res.p[pos_ids.z] = initial_pos.z;

    return res;
  }

  GRPrimitive create_start_primitive(const GROptimizationContext &ctx, const GRNode &node)
  {
    UPGStructure structure = get_structure_known_structure(ctx, node);
    glm::vec3 initial_pos = get_position_min_distance(node);
    UPGParametersRaw parameters = get_initial_parameters_random(node, structure, initial_pos);

    return {structure, parameters};
  }

  void set_opt_hyperparameters(const GROptimizationContext &ctx, FieldSdfCompare &opt_func)
  {
    opt_func.set_hyperparameters(ctx.settings.distance_base_thr, 2.0f, 2048, 100);
  }

  void set_fine_opt_hyperparameters(const GROptimizationContext &ctx, FieldSdfCompare &opt_func)
  {
    opt_func.set_hyperparameters(ctx.settings.distance_fine_thr, 2.0f, opt_func.a_points.size(), 100);
  }

  float estimate_positioning_quality(const UPGReconstructionResult &res, const FieldSdfCompare &opt_func, float threshold)
  {
    SdfGenInstance gen(res.structure);
    ProceduralSdf sdf = gen.generate(res.parameters.p);

    int inside_count = 0;
    int apr_count = 0;
    int inside_apr_count = 0;
    for (int i=0;i<opt_func.a_points.size();i++)
    {
      if (opt_func.a_distances[i] < 0)
        inside_count++;
      if (abs(sdf.get_distance(opt_func.a_points[i]) - opt_func.a_distances[i]) < threshold)
      {
        if (opt_func.a_distances[i] < 0)
          inside_apr_count++;
        apr_count++;
      }
    }

    return ((float)inside_apr_count/inside_count);
  }

  void activate_node_on_edge(const GROptimizationContext &ctx, GREdge &edge)
  {
    //fine optimization step 
    UPGReconstructionResult start_p;
    start_p.structure = edge.opt_prim.structure;
    start_p.parameters = edge.opt_prim.parameters;

    set_fine_opt_hyperparameters(ctx, edge.start->opt_func);
    std::shared_ptr<UPGOptimizer> optimizer = get_optimizer_adam(&(edge.start->opt_func), ctx.fine_opt_settings_blk, start_p);
    optimizer->optimize();
    auto partial_result = optimizer->get_best_results()[0];
    logerr("fine res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
                                     partial_result.parameters.p[2], partial_result.parameters.p[3]);
 
    //calculate edge quality for fine optimized element
    edge.fine_prim = GRPrimitive{partial_result.structure, partial_result.parameters};
    edge.fine_quality = estimate_positioning_quality(partial_result, edge.start->opt_func, ctx.settings.distance_fine_thr);

    //create new node with less active points
    SdfGenInstance gen(partial_result.structure);
    ProceduralSdf sdf = gen.generate(partial_result.parameters.p);

    edge.end.reset(new GRNode(edge.start->opt_func));
    edge.end->opt_func.remove_inactive_points(sdf, 0.001f); 
    edge.end->parent = edge.start;
    edge.end->depth = edge.start->depth+1;
    edge.end->quality = -1;//TODO
    edge.end->primitives = edge.start->primitives;
    edge.end->primitives.push_back(edge.fine_prim);
  }

  void opt_step(const GROptimizationContext &ctx, GRNode &node)
  {
    //get start parameters
    GRPrimitive start_primitive = create_start_primitive(ctx, node);

    //optimize 
    UPGReconstructionResult start_p;
    start_p.structure = start_primitive.structure;
    start_p.parameters = start_primitive.parameters;

    set_opt_hyperparameters(ctx, node.opt_func);
    std::shared_ptr<UPGOptimizer> optimizer = get_optimizer_adam(&(node.opt_func), ctx.opt_settings_blk, start_p);
    optimizer->optimize();
    auto partial_result = optimizer->get_best_results()[0];
    logerr("res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
                                partial_result.parameters.p[2], partial_result.parameters.p[3]);
  
    //calculate edge quality
    float result_quality = estimate_positioning_quality(partial_result, node.opt_func, ctx.settings.distance_fine_thr);
    logerr("Q = %f", result_quality);

    //create new edge
    GREdge new_edge;
    new_edge.start = &node;
    new_edge.start_prim = start_primitive;
    new_edge.opt_prim = GRPrimitive{partial_result.structure, partial_result.parameters};
    new_edge.opt_quality = result_quality;
    node.edges.push_back(new_edge);
  }

  GRNode *graph_step_N_greedy(const GROptimizationContext &ctx, GRNode &node, int N_tries = 50)
  {
    //if some tries left, stay in the same node and try again
    if (node.edges.size() < N_tries)
      return &node;

    //else choose best edge and go to it's end node
    int best_pos = 0;
    for (int i=0;i<node.edges.size();i++)
    {
      if (node.edges[i].opt_quality > node.edges[best_pos].opt_quality)
        best_pos = i;
    }
    if (!node.edges[best_pos].end)
      activate_node_on_edge(ctx, node.edges[best_pos]);
    logerr("best Q = %f", node.edges[best_pos].opt_quality);
    return node.edges[best_pos].end.get();
  }

  UPGReconstructionResult merge_primitives(const std::vector<GRPrimitive> &primitives)
  {
    std::vector<uint16_t> structure_inv;
    std::vector<float> params_inv;
    uint32_t num = 1;
    for (auto &pr : primitives)
    {
      for (int i=pr.structure.s.size()-1;i>=0;i--)
        structure_inv.push_back(pr.structure.s[i]);
      for (int i=pr.parameters.p.size()-1;i>=0;i--)
        params_inv.push_back(pr.parameters.p[i]);      
      int32_t S = 1;
      while (num>= S && ((num & S) == 0))
      {
        structure_inv.push_back(3); //merge node id
        S = S << 1;
      }
      num++;
    }
    int32_t S = 1;
    while (num>= S && ((num & S) == 0))
    {
      structure_inv.push_back(3); //merge node id
      S = S << 1;
    }

    UPGReconstructionResult res;
    for (int i=0;i<structure_inv.size();i++)
      res.structure.s.push_back(structure_inv[structure_inv.size()-i-1]);
    for (int i=0;i<params_inv.size();i++)
      res.parameters.p.push_back(params_inv[params_inv.size()-i-1]);
    debug("res structure {");
    for (auto &s : res.structure.s)
      debug("%d ",(int)s);
    debug("}\n");
    return res;
  }

  UPGReconstructionResult GR_agent_main_loop(GROptimizationContext &ctx)
  {
    int total_nodes = 0;
    int steps = 0;
    GRNode *final_node = nullptr;
    GRNode *node = ctx.root.get();

    int max_depth = ctx.target_structure_parts.size();
    int max_steps = 10000;

    while (node->depth < max_depth && steps < max_steps)
    {
      opt_step(ctx, *node);
      GRNode *n = graph_step_N_greedy(ctx, *node);
      if (node != n)
      {
        SdfGenInstance gen(n->primitives.back().structure);
        ProceduralSdf sdf = gen.generate(n->primitives.back().parameters.p);
        CameraSettings camera;
        camera.origin = glm::vec3(0,0,3);
        camera.target = glm::vec3(0,0,0);
        camera.up = glm::vec3(0,1,0);
        Texture t = render_sdf(sdf, camera, 512, 512, 4);
        engine::textureManager->save_png(t, "sdf_node_"+std::to_string(total_nodes));
        total_nodes++;
      }
      node = n;
      steps++;
    }

    return merge_primitives(node->primitives);
  }

  std::vector<UPGReconstructionResult> reconstruction_graph_based(Block *step_blk, const std::vector<glm::vec3> &points,
                                                                  const std::vector<float> &distances)
  {
    GROptimizationContext ctx;

    ctx.target_points = points;
    ctx.target_distances = distances;
    ctx.point_cloud_bbox = get_point_cloud_bbox(points);
    step_blk->get_arr("structure", ctx.target_structure.s);
    if (!ctx.target_structure.s.empty())
      ctx.target_structure_parts = get_sdf_parts(ctx.target_structure);

    ctx.opt_settings_blk.set_int("iterations", 1000);
    ctx.opt_settings_blk.set_bool("verbose", false);
    ctx.opt_settings_blk.set_double("learning_rate", 0.005f);

    ctx.fine_opt_settings_blk.set_int("iterations", 200);
    ctx.fine_opt_settings_blk.set_bool("verbose", false);
    ctx.fine_opt_settings_blk.set_double("learning_rate", 0.0005f);

    ctx.settings.distance_base_thr = 0.01f;
    ctx.settings.distance_fine_thr = 0.002f;

    ctx.root.reset(new GRNode(FieldSdfCompare(points, distances)));
    ctx.root->depth = 0;

    return {GR_agent_main_loop(ctx)};
  }
}