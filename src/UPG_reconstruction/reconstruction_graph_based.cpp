#include "reconstruction_graph_based.h"
#include "sdf_rendering.h"
#include "common_utils/bbox.h"
#include "common_utils/distribution.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include <algorithm>

namespace upg
{
  static int id_counter = 0;
  static float time_base_opt = 0;
  static float time_est_q = 0;
  static float time_fine_opt = 0;
  static float time_all = 0;
  static float time_search_and_select = 0;
  static float time_debug = 0;
  
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
      id = id_counter++;
    }
    FieldSdfCompare opt_func;
    std::vector<GRPrimitive> primitives;
    std::vector<GREdge> edges;
    GRNode *parent = nullptr;

    float quality = 0;
    int depth = 0;
    int id = 0;
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
    int max_number_of_primitives = 16;
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
    UPGParametersRaw target_parameters; //usually empty, is used only for testing purposes by so-called Oracles
    std::vector<UPGPart> target_structure_parts; //can be empty
    AABB point_cloud_bbox;
  };

  //[0.000]--(-0.11)-->[-0.11]
  //   |                
  //   ------(0.093)-->[0.093]--(0.112)-->
  //   |                  |
  //   |                  ------(0.011)-->
  //   |
  //   ------(0.007)-->[0.007]
  std::string get_node_string(const GRNode &node)
  {
    char ns[8];
    char op = node.edges.empty() ? '#' : '[';
    char cl = node.edges.empty() ? '#' : ']';
    if (node.quality >= 0)
      snprintf(ns, 8, "%c%.3f%c", op, MIN(1, node.quality), cl);
    else
      snprintf(ns, 8, "%c%.2f%c", op, MAX(-1, node.quality), cl);
    std::string n_string = node.edges.empty() ? std::to_string(node.id) : "";
    return std::string(ns) + n_string;
  }

  std::string get_edge_string(const GREdge &edge)
  {
    char es[13];
    if (edge.opt_quality >= 0)
      snprintf(es, 13, "--(%.3f)%s", MIN(1, edge.opt_quality), edge.end ? "-->" : "   ");
    else
      snprintf(es, 13, "--(%.2f)%s", MAX(-1, edge.opt_quality), edge.end ? "-->" : "   ");
    return std::string(es);
  }

  std::string get_v_string()
  {
    return "   |               ";
  }

  void print_reconstruction_node_rec(std::vector<std::string> &lines, const GRNode &node, int offset)
  {    
    std::string empty_line(offset, ' ');
    lines.back() += get_node_string(node);

    if (node.edges.empty())
      return;
    
    std::vector<unsigned> indices(node.edges.size());
    for (int i=0;i<node.edges.size();i++)
      indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&node](int a, int b)-> bool {
      if (!node.edges[a].end && node.edges[b].end)
        return false;
      else if (node.edges[a].end && !node.edges[b].end)
        return true;
      else if (node.edges[a].end && node.edges[b].end)
        return node.edges[a].end->id < node.edges[b].end->id;
      else
        return node.edges[a].opt_quality > node.edges[b].opt_quality;
    });

    lines.back() += get_edge_string(node.edges[indices[0]]);
    int prev_lines_size = lines.size();
    int next_row_offset = lines.back().size();
    if (node.edges[indices[0]].end)
      print_reconstruction_node_rec(lines, *(node.edges[indices[0]].end), next_row_offset);
    for (int i=1;i<node.edges.size();i++)
    {
      if (node.edges[indices[i]].end)
      {
        for (int l = prev_lines_size; l<lines.size();l++)
          lines[l].at(offset + 3) = '|';
        lines.emplace_back(); 
        lines.back() = empty_line + get_v_string();
        lines.emplace_back();
        lines.back() = empty_line + "   ----" + get_edge_string(node.edges[indices[i]]);
        prev_lines_size = lines.size();
        print_reconstruction_node_rec(lines, *(node.edges[indices[i]].end), next_row_offset);
      }
      else if (false)
      {
        if (lines.size() == prev_lines_size)
        {
          lines.emplace_back(); 
          lines.back() = empty_line + get_v_string();
          lines.emplace_back();
          lines.back() = empty_line + "   ----" + get_edge_string(node.edges[indices[i]]);
        }
        else
        {
          std::string l1 = get_v_string();
          for (int c=0;c<l1.size();c++)
            lines[prev_lines_size].at(offset + c) = l1.at(c);
          std::string l2 = "   ----" + get_edge_string(node.edges[indices[i]]);
          for (int c=0;c<l2.size();c++)
            lines[prev_lines_size+1].at(offset + c) = l2.at(c);
        }
        prev_lines_size += 2;
      }
    }
  }

  void print_reconstruction_graph(const GROptimizationContext &ctx, const GRNode &root)
  {
    debug("RECONSTRUCTION GRAPH\n");
    for (int i=0;i<80;i++)
      debug("=");
    debug("\n");
    std::vector<std::string> lines = {""};
    print_reconstruction_node_rec(lines, root, 0);
    for (std::string &l : lines)
      debug("%s\n",l.c_str());
    for (int i=0;i<80;i++)
      debug("=");
    debug("\n");
  }

  UPGStructure get_structure_known_structure(const GROptimizationContext &ctx, const GRNode &node)
  {
    assert(ctx.target_structure_parts.size() > 0);
    UPGStructure part_structure;
    auto &part = ctx.target_structure_parts[node.depth];
    for (int j=part.s_range.first; j<part.s_range.second; j++)
      part_structure.s.push_back(ctx.target_structure.s[j]);
    //logerr("depth %d created structure %d %d", node.depth, (int)part_structure.s[0], (int)part_structure.s[1]);
    return part_structure;
  }
  UPGStructure get_random_structure(const GROptimizationContext &ctx, const GRNode &node)
  {
    std::vector<UPGStructure> available_primitives = {{{SdfNode::MOVE, SdfNode::SPHERE}},
                                                      {{SdfNode::MOVE, SdfNode::BOX}},
                                                      {{SdfNode::MOVE, SdfNode::ROUNDED_BOX}},
                                                      {{SdfNode::MOVE, SdfNode::CYLINDER}},
                                                      {{SdfNode::MOVE, SdfNode::CONE}},
                                                      {{SdfNode::MOVE, SdfNode::PRISM}}};
    return available_primitives[urandi(0, available_primitives.size())];
  }
  UPGStructure get_structure_oracle(const GROptimizationContext &ctx, const GRNode &node, float error_chance)
  {
    if (urand() < error_chance)
      return get_random_structure(ctx, node);
    else
      return get_structure_known_structure(ctx, node);
  }

  glm::vec3 get_position_min_distance(const GROptimizationContext &ctx, const GRNode &node, int tries = 25)
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

  glm::ivec3 get_pos_ids(const GRNode &node, const UPGStructure &structure)
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
    return pos_ids;
  }

  std::vector<glm::vec2> get_borders(const GRNode &node, const UPGStructure &structure)
  {
    auto gen = node.opt_func.get_generator(structure);
    auto pd = node.opt_func.get_full_parameters_description(gen.get());

    std::vector<glm::vec2> borders;
    for (const auto &p : pd.get_block_params())
    {
      for (const auto &param_info : p.second.p)
      {
        if (param_info.type != ParameterType::CONST)
          borders.push_back(glm::vec2(param_info.min_val, param_info.max_val));
      }
    }
    return borders;
  }

  UPGParametersRaw get_initial_parameters_known_parameters(const GROptimizationContext &ctx, const GRNode &node, const UPGStructure &structure,
                                                           glm::vec3 initial_pos);
  glm::vec3 get_position_oracle(const GROptimizationContext &ctx, const GRNode &node, float error_std_dev)
  {
    Normal norm(0, error_std_dev);

    UPGStructure ref_s = get_structure_known_structure(ctx, node);
    UPGParametersRaw ref_p = get_initial_parameters_known_parameters(ctx, node, ref_s, {0,0,0});
    glm::ivec3 pos_ids = get_pos_ids(node, ref_s);

    glm::vec3 p(ref_p.p[pos_ids.x] + norm.get(), ref_p.p[pos_ids.y] + norm.get(), ref_p.p[pos_ids.z] + norm.get());
    //logerr("pos = %f %f %f", p.x, p.y, p.z);
    return p;
  }

  UPGParametersRaw get_initial_parameters_random(const GROptimizationContext &ctx, const GRNode &node, const UPGStructure &structure, 
                                                 glm::vec3 initial_pos)
  {
    UPGParametersRaw res;
    glm::ivec3 pos_ids = get_pos_ids(node, structure);
    std::vector<glm::vec2> borders = get_borders(node, structure);
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

  UPGParametersRaw get_initial_parameters_known_parameters(const GROptimizationContext &ctx, const GRNode &node, const UPGStructure &structure,
                                                           glm::vec3 initial_pos)
  {
    UPGStructure ref_s = get_structure_known_structure(ctx, node);
    bool same_structure = structure.s.size() == ref_s.s.size();
    if (same_structure)
    {
      for (int i=0;i<structure.s.size();i++)
      {
        if (structure.s[i] != ref_s.s[i])
        {
          same_structure = false;
          break;
        }
      }
    }

    if (same_structure) //pick reference parameters
    {
      UPGParametersRaw part_params;
      auto &part = ctx.target_structure_parts[node.depth];
      for (int j=part.p_range.first; j<part.p_range.second; j++)
        part_params.p.push_back(ctx.target_parameters.p[j]);
      return part_params;
    }
    else //pick random parameters
      return get_initial_parameters_random(ctx, node, structure, initial_pos);
  }

  void set_opt_hyperparameters(const GROptimizationContext &ctx, FieldSdfCompare &opt_func)
  {
    opt_func.set_hyperparameters(ctx.settings.distance_base_thr, 2.0f, 512, 100);
  }

  void set_fine_opt_hyperparameters(const GROptimizationContext &ctx, FieldSdfCompare &opt_func)
  {
    opt_func.set_hyperparameters(ctx.settings.distance_fine_thr, 2.0f, opt_func.a_points.size(), 100);
  }

  float estimate_positioning_quality(const UPGReconstructionResult &res, const FieldSdfCompare &opt_func, float threshold)
  {
    ProceduralSdf sdf(res.structure);
    sdf.set_parameters(res.parameters.p);

    int inside_count = 0;
    int apr_count = 0;
    int inside_apr_count = 0;
    int swallowed_points = 0;
    for (int i=0;i<opt_func.a_points.size();i++)
    {
      //if (opt_func.a_distances[i] < 0)
      //  inside_count++;
      if (abs(sdf.get_distance(opt_func.a_points[i]) - opt_func.a_distances[i]) < threshold)
      {
        if (opt_func.a_distances[i] < 0)
          inside_apr_count++;
        apr_count++;
      }
      else if (sdf.get_distance(opt_func.a_points[i]) < 0 && opt_func.a_distances[i] > 0)
        swallowed_points++;
    }

    for (int i=0;i<opt_func.points.size();i++)
    {
      if (opt_func.distances[i] < 0)
        inside_count++;
    }

    return ((float)(inside_apr_count - 100*swallowed_points)/inside_count);
  }

  float estimate_solution_quality_MAE(const UPGReconstructionResult &res, const FieldSdfCompare &opt_func)
  {
    ProceduralSdf sdf(res.structure);
    sdf.set_parameters(res.parameters.p);

    double AE = 0;
    for (int i=0;i<opt_func.points.size();i++)
    {
      AE += abs(opt_func.distances[i] - sdf.get_distance(opt_func.points[i]));
    }

    return MAX(0, 1-AE/opt_func.points.size());
  }

  void activate_node_on_edge(const GROptimizationContext &ctx, GREdge &edge)
  {
std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    //fine optimization step 
    UPGReconstructionResult start_p;
    start_p.structure = edge.opt_prim.structure;
    start_p.parameters = edge.opt_prim.parameters;

    set_fine_opt_hyperparameters(ctx, edge.start->opt_func);
    std::shared_ptr<UPGOptimizer> optimizer = get_optimizer_adam(&(edge.start->opt_func), ctx.fine_opt_settings_blk, start_p);
    optimizer->optimize();
    auto partial_result = optimizer->get_best_results()[0];
    //logerr("fine res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
    //                                 partial_result.parameters.p[2], partial_result.parameters.p[3]);
 
    //calculate edge quality for fine optimized element
    edge.fine_prim = GRPrimitive{partial_result.structure, partial_result.parameters};
    edge.fine_quality = estimate_positioning_quality(partial_result, edge.start->opt_func, ctx.settings.distance_fine_thr);

    //create new node with less active points
    ProceduralSdf sdf(partial_result.structure);
    sdf.set_parameters(partial_result.parameters.p);

    edge.end.reset(new GRNode(edge.start->opt_func));
    edge.end->opt_func.remove_inactive_points(sdf, 0.001f); 
    edge.end->parent = edge.start;
    edge.end->depth = edge.start->depth+1;
    edge.end->quality = edge.fine_quality;
    edge.end->primitives = edge.start->primitives;
    edge.end->primitives.push_back(edge.fine_prim);
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      time_fine_opt += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  }

  void opt_step(const GROptimizationContext &ctx, GRNode &node, int tries = 1)
  {
    set_opt_hyperparameters(ctx, node.opt_func);
    int offset = node.edges.size();
    node.edges.resize(offset+tries);
    
    std::vector<FieldSdfCompare> opt_funcs;
    std::vector<GRPrimitive> start_primitives;
    for (int i=0;i<tries;i++)
    {
      opt_funcs.emplace_back(node.opt_func);
      UPGStructure structure = get_structure_oracle(ctx, node, 0.0);
      glm::vec3 initial_pos = get_position_oracle(ctx, node, 0.0);
      UPGParametersRaw parameters = get_initial_parameters_random(ctx, node, structure, initial_pos);
      start_primitives.push_back({structure, parameters});
    }

    #pragma omp parallel for 
    for (int i=0;i<tries;i++)
    {
      //optimize 
      UPGReconstructionResult start_p;
      start_p.structure = start_primitives[i].structure;
      start_p.parameters = start_primitives[i].parameters;

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
      std::shared_ptr<UPGOptimizer> optimizer = get_optimizer_multistep_adam(&(opt_funcs[i]), ctx.opt_settings_blk, start_p);
      optimizer->optimize();
      auto partial_result = optimizer->get_best_results()[0];
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      time_base_opt += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      //logerr("res = %f %f %f %f", partial_result.parameters.p[0], partial_result.parameters.p[1],
      //                            partial_result.parameters.p[2], partial_result.parameters.p[3]);

std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();  
      //calculate edge quality
      float result_quality = estimate_positioning_quality(partial_result, opt_funcs[i], ctx.settings.distance_fine_thr);
      //logerr("%d Q = %f", i, result_quality);

      //create new edge
      GREdge new_edge;
      new_edge.start = &node;
      new_edge.start_prim = start_primitives[i];
      new_edge.opt_prim = GRPrimitive{partial_result.structure, partial_result.parameters};
      new_edge.opt_quality = result_quality;
      node.edges[offset+i] = new_edge;
std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
      time_est_q += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
    }
  }

  UPGReconstructionResult merge_primitives(const std::vector<GRPrimitive> &primitives)
  {
    assert(primitives.size() > 0);
    if (primitives.size() == 1)
    {
      UPGReconstructionResult res;
      res.structure = primitives[0].structure;
      res.parameters = primitives[0].parameters;
      return res;
    }
    UPGReconstructionResult res;
    for (int i=0;i<primitives.size();i++)
    {
      if (i != primitives.size()-1)
        res.structure.s.push_back(SdfNode::OR);
      for (auto &s : primitives[i].structure.s)
        res.structure.s.push_back(s);
      for (auto &p : primitives[i].parameters.p)
        res.parameters.p.push_back(p);
    }
    debug("res structure {");
    for (auto &s : res.structure.s)
      debug("%d ",(int)s);
    debug("}\n");
    return res;
  }

  void debug_render_node(const GRNode *n, int id)
  {
std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    ProceduralSdf sdf(n->primitives.back().structure);
    sdf.set_parameters(n->primitives.back().parameters.p);
    CameraSettings camera;
    camera.origin = glm::vec3(0, 0, 3);
    camera.target = glm::vec3(0, 0, 0);
    camera.up = glm::vec3(0, 1, 0);
    Texture t = render_sdf(sdf, camera, 512, 512, 4);
    engine::textureManager->save_png(t, "sdf_node_" + std::to_string(id));
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      time_debug += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  }

  UPGReconstructionResult GR_agent_N_greedy(GROptimizationContext &ctx, int N_tries = 25)
  {
    int total_nodes = 0;
    int steps = 0;
    GRNode *final_node = nullptr;
    GRNode *node = ctx.root.get();

    int max_depth = ctx.settings.max_number_of_primitives;
    int max_steps = 10000;

    while (node->depth < max_depth && steps < max_steps)
    {
      opt_step(ctx, *node);
      GRNode *new_node = node;

      //choose best edge and go to the next stop with it
      if (node->edges.size() >= N_tries)
      {
        int best_pos = 0;
        for (int i=0;i<node->edges.size();i++)
        {
          if (node->edges[i].opt_quality > node->edges[best_pos].opt_quality)
            best_pos = i;
        }
        if (!node->edges[best_pos].end)
          activate_node_on_edge(ctx, node->edges[best_pos]);
        //logerr("best Q = %f", node->edges[best_pos].opt_quality);
        new_node = node->edges[best_pos].end.get();
      }
      if (node != new_node)
      {
        debug_render_node(new_node, total_nodes);
        total_nodes++;
      }
      node = new_node;
      steps++;
    }

    return merge_primitives(node->primitives);
  }

  UPGReconstructionResult GR_agent_stochastic_returns(const GROptimizationContext &ctx, int path_limit = 25, int tries_per_step = 3)
  {
    float best_quality = -1e9;
    UPGReconstructionResult best_result;
    GRNode *best_end_node = nullptr;
    GRNode *node = ctx.root.get();
    int max_depth = ctx.settings.max_number_of_primitives;
    int total_nodes = 0;
    int paths = 0;

    while (paths < path_limit /*&& best_quality < (1-1e-5)*/)
    {
      GRNode *new_node = node;

      //pick best and go deeper
      if (node->depth < max_depth)
      {
        opt_step(ctx, *node, tries_per_step);
        int best_pos = 0;
        for (int i=0;i<node->edges.size();i++)
        {
          if (node->edges[i].opt_quality > node->edges[best_pos].opt_quality)
            best_pos = i;
        }
        if (!node->edges[best_pos].end)
          activate_node_on_edge(ctx, node->edges[best_pos]);
        //logerr("best Q = %f", node->edges[best_pos].opt_quality);
        new_node = node->edges[best_pos].end.get();
      }
      else //we found some solution. Remember it if it's the best and go to some previous node to improve later
      {
std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        UPGReconstructionResult res = merge_primitives(node->primitives);
        float solution_quality = estimate_solution_quality_MAE(res, node->opt_func);
        node->quality = solution_quality;
        if (solution_quality > best_quality)
        {
          best_quality = solution_quality;
          best_end_node = node;
          best_result = res;
        }
        paths++;

        std::vector<float> improvement_chances;
        std::vector<float> qualities;
        std::vector<GRNode *> path;
        GRNode *n = best_end_node;
        while (n->parent != nullptr)
        {
          //for every node chance of improvement is inverse proportional to it's best edge quality
          path.push_back(n->parent);
          float prev_chance = improvement_chances.empty() ? 0 : improvement_chances.back();
          improvement_chances.push_back(1.0f/MAX(1e-6, n->quality) + prev_chance);
          qualities.push_back(n->quality);
          n = n->parent;
        }

        float rnd = urand(0, improvement_chances.back());
        int i_node_id = 0;
        for (int i=0;i<improvement_chances.size();i++)
        {
          if (rnd < improvement_chances[i])
          {
            i_node_id = i;
            break;
          } 
        }

        logerr("improving node %d(path %d) with chance %f", i_node_id, best_end_node->id, qualities[i_node_id]);
        //print_reconstruction_graph(ctx, *(ctx.root));
        new_node = path[i_node_id];
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      time_search_and_select += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      if (node != new_node && new_node->depth > 0)
      {
        debug_render_node(new_node, total_nodes);
        total_nodes++;
      }
      node = new_node;
    }

    return best_result;
  }

  UPGReconstructionResult GR_agent_stochastic_changes(const GROptimizationContext &ctx, int path_limit = 25, int tries_per_step = 3)
  {
    float best_quality = -1e9;
    UPGReconstructionResult best_result;
    GRNode *best_end_node = nullptr;
    GRNode *base_path_node = nullptr;
    GRNode *node = ctx.root.get();
    int max_depth = ctx.settings.max_number_of_primitives;
    int total_nodes = 0;
    int paths = 0;

    while (paths < path_limit /*&& best_quality < (1-1e-5)*/)
    {
      GRNode *new_node = node;
      GRNode *new_base_path_node = nullptr;
      if (base_path_node)
      {
        int best_pos = 0;
        for (int i=0;i<base_path_node->edges.size();i++)
        {
          if (base_path_node->edges[i].opt_quality > base_path_node->edges[best_pos].opt_quality)
            best_pos = i;
        }
        if (base_path_node->edges.size() > 0 && base_path_node->edges[best_pos].end)
          new_base_path_node = base_path_node->edges[best_pos].end.get();
        else
          new_base_path_node = nullptr;
      }

      //pick best and go deeper
      if (node->depth < max_depth)
      {
        if (base_path_node == nullptr || base_path_node == node)
        {
          opt_step(ctx, *node, tries_per_step);
        }
        else
        {
          set_opt_hyperparameters(ctx, node->opt_func);
          GRPrimitive start_primitive = base_path_node->primitives.back();
          
          //optimize 
          UPGReconstructionResult start_p;
          start_p.structure = start_primitive.structure;
          start_p.parameters = start_primitive.parameters;

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
          std::shared_ptr<UPGOptimizer> optimizer = get_optimizer_multistep_adam(&(node->opt_func), ctx.opt_settings_blk, start_p);
          optimizer->optimize();
          auto partial_result = optimizer->get_best_results()[0];
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
          time_base_opt += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();  
          //calculate edge quality
          float result_quality = estimate_positioning_quality(partial_result, node->opt_func, ctx.settings.distance_fine_thr);
          logerr("Qa = %f", result_quality);

          //create new edge
          GREdge new_edge;
          new_edge.start = node;
          new_edge.start_prim = start_primitive;
          new_edge.opt_prim = GRPrimitive{partial_result.structure, partial_result.parameters};
          new_edge.opt_quality = result_quality;
          node->edges.push_back(new_edge);
std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
          time_est_q += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
        }

        base_path_node = new_base_path_node;

        int best_pos = 0;
        for (int i=0;i<node->edges.size();i++)
        {
          if (node->edges[i].opt_quality > node->edges[best_pos].opt_quality)
            best_pos = i;
        }
        if (!node->edges[best_pos].end)
          activate_node_on_edge(ctx, node->edges[best_pos]);
        //logerr("best Q = %f", node->edges[best_pos].opt_quality);
        new_node = node->edges[best_pos].end.get();
      }
      else //we found some solution. Remember it if it's the best and go to some previous node to improve later
      {
std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        UPGReconstructionResult res = merge_primitives(node->primitives);
        float solution_quality = estimate_solution_quality_MAE(res, node->opt_func);
        node->quality = solution_quality;
        if (solution_quality > best_quality)
        {
          best_quality = solution_quality;
          best_end_node = node;
          best_result = res;
        }

        std::vector<float> improvement_chances;
        std::vector<float> qualities;
        std::vector<GRNode *> path;
        GRNode *n = best_end_node;
        while (n->parent != nullptr)
        {
          //for every node chance of improvement is inverse proportional to it's best edge quality
          path.push_back(n->parent);
          float prev_chance = improvement_chances.empty() ? 0 : improvement_chances.back();
          improvement_chances.push_back(1.0f/MAX(1e-6, n->quality) + prev_chance);
          qualities.push_back(n->quality);
          n = n->parent;
        }

        float rnd = urand(0, improvement_chances.back());
        int i_node_id = 0;
        for (int i=0;i<improvement_chances.size();i++)
        {
          if (rnd < improvement_chances[i])
          {
            i_node_id = i;
            break;
          } 
        }

        logerr("improving node %d(path %d) with chance %f", i_node_id, best_end_node->id, qualities[i_node_id]);
        print_reconstruction_graph(ctx, *(ctx.root));
        new_node = path[i_node_id];
        base_path_node = new_node;
        paths++;
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      time_search_and_select += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      }

      if (node != new_node && new_node->depth > 0)
      {
        debug_render_node(new_node, total_nodes);
        total_nodes++;
      }
      node = new_node;
    }

    return best_result;
  }

  std::vector<UPGReconstructionResult> reconstruction_graph_based(Block *step_blk, const std::vector<glm::vec3> &points,
                                                                  const std::vector<float> &distances)
  {
    GROptimizationContext ctx;

    ctx.target_points = points;
    ctx.target_distances = distances;
    ctx.point_cloud_bbox = get_point_cloud_bbox(points);
    step_blk->get_arr("structure", ctx.target_structure.s);
    step_blk->get_arr("params", ctx.target_parameters.p);
    if (!ctx.target_structure.s.empty())
      ctx.target_structure_parts = get_sdf_parts(ctx.target_structure);

    ctx.opt_settings_blk.set_int("iterations", 1000);
    ctx.opt_settings_blk.set_bool("verbose", false);
    ctx.opt_settings_blk.set_double("learning_rate", 0.01f);
    ctx.opt_settings_blk.set_vec3("step_params_0", glm::vec3( 100,  0.01,   1e9));
    ctx.opt_settings_blk.set_vec3("step_params_1", glm::vec3( 500,  0.01,  0.00));
    ctx.opt_settings_blk.set_vec3("step_params_2", glm::vec3(1000, 0.005, -0.01));

    ctx.fine_opt_settings_blk.set_int("iterations", 200);
    ctx.fine_opt_settings_blk.set_bool("verbose", false);
    ctx.fine_opt_settings_blk.set_double("learning_rate", 0.0005f);

    ctx.settings.distance_base_thr = 0.01f;
    ctx.settings.distance_fine_thr = 0.002f;
    ctx.settings.max_number_of_primitives = ctx.target_structure.s.empty() ? 16 : ctx.target_structure_parts.size();

    ctx.root.reset(new GRNode(FieldSdfCompare(points, distances)));
    ctx.root->depth = 0;

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    auto res = GR_agent_stochastic_changes(ctx, 30, 5);
std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    time_all += 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    print_reconstruction_graph(ctx, *(ctx.root));

    /*static float time_base_opt = 0;
    static float time_est_q = 0;
    static float time_fine_opt = 0;
    static float time_all = 0;
    static float time_search_and_select = 0;
    static float time_debug = 0*/
    logerr("TIME %f %f %f %f %f %f", time_base_opt, time_est_q, time_fine_opt, time_search_and_select, time_debug, time_all);

    return {res};
  }
}