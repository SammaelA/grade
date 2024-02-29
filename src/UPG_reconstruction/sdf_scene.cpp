#include "sdf_scene.h"
#include "sdf_node.h"

namespace upg
{
  struct LogicNode
  {
    unsigned type;
    unsigned left_id;
    unsigned right_id;
    bool is_literal;
    int literal_id = -1;
  };
  void get_transform_rec(const SdfNode *node, 
                         std::vector<LogicNode> &nodes, 
                         std::vector<float> &parameters,
                         std::vector<SdfObject> &objects,
                         glm::mat4 transform, float distance_mult, float distance_add)
  {
    logerr("GTR %lu %d",(uint64_t)node, node->child_cnt());
    SdfNodeClass cl = get_sdf_node_properties(node->type).node_class;
    int param_cnt = get_sdf_node_properties(node->type).param_count;
    switch (cl)
    {
    case SdfNodeClass::PRIMITIVE:
    {
      objects.emplace_back();
      objects.back().type = node->type;
      objects.back().bbox = node->get_bbox();
      objects.back().distance_add = distance_add;
      objects.back().distance_mult = distance_mult;
      objects.back().params_count = param_cnt;
      objects.back().params_offset = parameters.size();
      objects.back().transform = transform;

      parameters.insert(parameters.end(), node->p.begin(), node->p.end());

      nodes.push_back({(unsigned)node->type, 0u, 0u, true, (int)(objects.size()-1)});
    }
      break;
    case SdfNodeClass::TRANSFORM:
    {
      assert(node->child_cnt() == 1);
      switch (node->type)
      {
      case SdfNodeType::MOVE:
        transform = glm::translate(transform, glm::vec3(-node->p[0], -node->p[1], -node->p[2]));
        break;
      case SdfNodeType::ROTATE:
        {
          float x = cosf(node->p[0]) * cosf(node->p[1]);
          float y = sinf(node->p[0]) * cosf(node->p[1]);
          float z = sinf(node->p[1]);
          glm::vec3 axis = normalize(glm::vec3{x,y,z});
          transform = glm::rotate(transform, -node->p[2], axis);
        }
        break;
      case SdfNodeType::SCALE:
        transform = glm::scale(transform, glm::vec3(1.0f/node->p[0], 1.0f/node->p[0], 1.0f/node->p[0]));
        distance_mult *= node->p[0];
        break;
      case SdfNodeType::ROUND:
        distance_add -= distance_mult*node->p[0];
        ///distance
        break;      
      default:
        logerr("unknown node %s from TRANSFORM class\n", get_sdf_node_properties(node->type).name.c_str());
        break;
      }
      get_transform_rec(node->get_children()[0], nodes, parameters, objects, transform, distance_mult, distance_add);
    }
      break;    
    case SdfNodeClass::COMBINE:
    {
      assert(node->child_cnt() == 2);
      nodes.push_back({(unsigned)node->type, 0u, 0u, false});
      unsigned cur_node_id = nodes.size()-1;
      
      nodes[cur_node_id].left_id = nodes.size();
      get_transform_rec(node->get_children()[0], nodes, parameters, objects, transform, distance_mult, distance_add);
      
      nodes[cur_node_id].right_id = nodes.size();
      get_transform_rec(node->get_children()[1], nodes, parameters, objects, transform, distance_mult, distance_add);
    }
      break;
    default:
      logerr("Sdf node %s belongs to the unknown or unsupported class %d", get_sdf_node_properties(node->type).name.c_str(), (int)cl);
      break;
    }
  }

  unsigned duplicate_subtree_rec(std::vector<LogicNode> &nodes, unsigned orig_id)
  {
    unsigned new_id = nodes.size();
    nodes.push_back({nodes[orig_id].type, 0, 0, nodes[orig_id].is_literal, nodes[orig_id].literal_id});  
    if (!nodes[orig_id].is_literal)
    {
      nodes[new_id].left_id = duplicate_subtree_rec(nodes, nodes[orig_id].left_id);
      nodes[new_id].right_id = duplicate_subtree_rec(nodes, nodes[orig_id].right_id);
    }
    logerr("created node %d", (int)new_id);
    return new_id;
  }

  void normalize_node_tree(std::vector<LogicNode> &nodes, unsigned cur_id)
  {
    if (nodes[cur_id].is_literal)
      return;
    constexpr unsigned SUB = SdfNodeType::SUBTRACT;
    constexpr unsigned AND = SdfNodeType::AND;
    constexpr unsigned OR  = SdfNodeType::OR;
    bool node_finished = false;
    while (!node_finished)
    {
      bool match = true;
      while (match)
      {
        #define T (nodes[cur_id])
        #define L (nodes[nodes[cur_id].left_id])
        #define R (nodes[nodes[cur_id].right_id])
        #define  L_id (nodes[cur_id].left_id)
        #define  R_id (nodes[cur_id].right_id)
        match = false;

debug("START ([%u]%u %u %u)([%u]%u %u %u)([%u]%u %u %u) [\n",
nodes[cur_id].type, cur_id, nodes[cur_id].left_id, nodes[cur_id].right_id, 
nodes[nodes[cur_id].left_id].type, nodes[cur_id].left_id, nodes[nodes[cur_id].left_id].left_id, nodes[nodes[cur_id].left_id].right_id,
nodes[nodes[cur_id].right_id].type, nodes[cur_id].right_id, nodes[nodes[cur_id].right_id].left_id, nodes[nodes[cur_id].right_id].right_id);

        if (T.type == SUB && R.type == OR)
        {
          unsigned x_id = L_id; 
          unsigned y_id = R.left_id;
          unsigned z_id = R.right_id;

          R = {SUB, x_id, y_id, false};
          T = {SUB, R_id, z_id, false};

          match = true;
        }
        else if (T.type == AND && R.type == OR)
        {
          unsigned x_id = L_id; 
          unsigned y_id = R.left_id;
          unsigned z_id = R.right_id;    

          nodes.push_back({AND, duplicate_subtree_rec(nodes, x_id), y_id, false});   
          unsigned new_id = nodes.size()-1;
          logerr("1created new node %u",new_id);

          R = {AND, x_id, z_id, false};
          T = {OR, new_id, R_id, false};
          match = true;

        }
        else if (T.type == SUB && R.type == AND)
        {
          unsigned x_id = L_id; 
          unsigned y_id = R.left_id;
          unsigned z_id = R.right_id;    

          nodes.push_back({SUB, duplicate_subtree_rec(nodes, x_id), y_id, false});  
          unsigned new_id = nodes.size()-1;

          R = {SUB, x_id, z_id, false};
          T = {OR, new_id, R_id, false};
          match = true;
        }
        else if (T.type == AND && R.type == AND)
        {
          unsigned x_id = L_id; 
          unsigned y_id = R.left_id;
          unsigned z_id = R.right_id; 
          
          R = {AND, x_id, y_id, false};
          T = {AND, R_id, z_id, false};
          match = true;
        }
        else if (T.type == SUB && R.type == SUB)
        {
          unsigned x_id = L_id; 
          unsigned y_id = R.left_id;
          unsigned z_id = R.right_id;    

          nodes.push_back({SUB, duplicate_subtree_rec(nodes, x_id), y_id, false}); 
          unsigned new_id = nodes.size()-1;

          R = {AND, x_id, z_id, false};
          T = {OR, new_id, R_id, false};  
          match = true;     
        }
        else if (T.type == AND && R.type == SUB)
        {
          unsigned x_id = L_id; 
          unsigned y_id = R.left_id;
          unsigned z_id = R.right_id;
          
          R = {AND, x_id, y_id, false};
          T = {SUB, R_id, z_id, false}; 
          match = true;              
        }
        else if (T.type == SUB && L.type == OR)
        {
          unsigned x_id = L.left_id; 
          unsigned y_id = L.right_id;
          unsigned z_id = R_id;   

          nodes.push_back({SUB, y_id, duplicate_subtree_rec(nodes, z_id), false});  
          unsigned new_id = nodes.size()-1;

          L = {SUB, x_id, z_id, false};
          T = {OR, L_id, new_id, false}; 
          match = true;         
        }
        else if (T.type == AND && L.type == OR)
        {
          unsigned x_id = L.left_id; 
          unsigned y_id = L.right_id;
          unsigned z_id = R_id;   

          nodes.push_back({AND, y_id, duplicate_subtree_rec(nodes, z_id), false});  
          unsigned new_id = nodes.size()-1;
          logerr("2creaeted new node %u",new_id);

          L = {AND, x_id, z_id, false};
          T = {OR, L_id, new_id, false};  

          match = true;        
        }
        else if (T.type == AND && L.type == SUB)
        {
          unsigned x_id = L.left_id; 
          unsigned y_id = L.right_id;
          unsigned z_id = R_id;   

          L = {AND, x_id, z_id, false};
          T = {SUB, L_id, y_id, false}; 
          match = true;         
        }
debug("] ([%u]%u %u %u)([%u]%u %u %u)([%u]%u %u %u) END\n",
nodes[cur_id].type, cur_id, nodes[cur_id].left_id, nodes[cur_id].right_id, 
nodes[nodes[cur_id].left_id].type, nodes[cur_id].left_id, nodes[nodes[cur_id].left_id].left_id, nodes[nodes[cur_id].left_id].right_id,
nodes[nodes[cur_id].right_id].type, nodes[cur_id].right_id, nodes[nodes[cur_id].right_id].left_id, nodes[nodes[cur_id].right_id].right_id);
      }

      normalize_node_tree(nodes, nodes[cur_id].left_id);

      node_finished = (T.type == OR) || (R.is_literal && L.type != OR);
    }
    normalize_node_tree(nodes, nodes[cur_id].right_id);
  }

  void get_conjunction_nodes(std::vector<LogicNode> &nodes, std::vector<unsigned> &conj_node_ids, unsigned cur_id)
  {
    if (nodes[cur_id].type == SdfNodeType::OR)
    {
      get_conjunction_nodes(nodes, conj_node_ids, nodes[cur_id].left_id);
      get_conjunction_nodes(nodes, conj_node_ids, nodes[cur_id].right_id);
    }
    else
      conj_node_ids.push_back(cur_id);
  }

  void get_literal_nodes(std::vector<LogicNode> &nodes, std::vector<std::pair<unsigned, bool>> &literal_node_ids, unsigned cur_id, bool complement)
  {
    if (nodes[cur_id].is_literal)
    {
      literal_node_ids.push_back({cur_id, complement});
    }
    else if (nodes[cur_id].type == SdfNodeType::AND)
    {
      get_literal_nodes(nodes, literal_node_ids, nodes[cur_id].left_id, complement);
      get_literal_nodes(nodes, literal_node_ids, nodes[cur_id].right_id, complement);
    } 
    else if (nodes[cur_id].type == SdfNodeType::SUBTRACT)
    {
      get_literal_nodes(nodes, literal_node_ids, nodes[cur_id].left_id, complement);
      get_literal_nodes(nodes, literal_node_ids, nodes[cur_id].right_id, !complement);      
    }
    else
      assert(false);
  }

  AABB transform_bbox(AABB bbox, glm::mat4 transform)
  {
    glm::mat4 itr = glm::inverse(transform);
    glm::vec3 p_min(1e9,1e9,1e9);
    glm::vec3 p_max(-1e9,-1e9,-1e9);
    for (unsigned i=0;i<8;i++)
    {
      glm::vec3 idx = glm::vec3(i&4>>2,i&2>>1,i&1);
      glm::vec3 corner = itr*glm::vec4(idx*bbox.min_pos + (glm::vec3(1,1,1)-idx)*bbox.max_pos, 1.0f);
      logerr("corner %f %f %f", corner.x, corner.y, corner.z);
      //corner = itr*glm::vec4(1,1,1,1);
      p_min = min(p_min, corner);
      p_max = max(p_max, corner);
    }
    return AABB(p_min, p_max);
  }

  SdfScene create_sdf_scene(const UPGStructure &structure, const UPGParametersRaw &params)
  {
    ProceduralSdf reference_sdf(structure);
    reference_sdf.set_parameters(params.p);

    std::vector<LogicNode> nodes;
    std::vector<float> all_parameters;
    std::vector<SdfObject> basic_objects;
    get_transform_rec(reference_sdf.root, nodes, all_parameters, basic_objects, glm::mat4(1.0f), 1.0f, 0.0f);
    int i=0;
    //for (auto &n : nodes)
    //  logerr("%d node t%u (%u %u) %d",i++, n.type, n.left_id, n.right_id, (int)n.literal_id);    
    normalize_node_tree(nodes, 0);
    i = 0;
    for (auto &n : nodes)
      logerr("%d node t%u (%u %u) %d",i++, n.type, n.left_id, n.right_id, (int)n.literal_id);    
    
    SdfScene scene;
    scene.parameters = all_parameters;

    std::vector<unsigned> conj_node_ids;
    get_conjunction_nodes(nodes, conj_node_ids, 0);
    for (auto &id : conj_node_ids)
    {
      logerr("Conjuction %u",id);
      std::vector<std::pair<unsigned, bool>> literal_node_ids;
      get_literal_nodes(nodes, literal_node_ids, id, false);
      scene.conjunctions.push_back({(unsigned)scene.objects.size(), (unsigned)literal_node_ids.size()});
      for (auto &p : literal_node_ids)
      {
        SdfObject &base_obj = basic_objects[nodes[p.first].literal_id];
        scene.objects.emplace_back();
        scene.objects.back().type = base_obj.type;
        scene.objects.back().params_offset = base_obj.params_offset;
        scene.objects.back().params_count = base_obj.params_count;
        scene.objects.back().distance_mult = base_obj.distance_mult;
        scene.objects.back().distance_add = base_obj.distance_add;
        scene.objects.back().bbox = transform_bbox(base_obj.bbox, base_obj.transform);
        scene.objects.back().transform = base_obj.transform;
        scene.objects.back().complement = p.second;
        logerr("%u (prim %d) %s", p.first, nodes[p.first].literal_id, p.second ? "COMPLEMENT" : "");
        logerr("transform (%f %f %f)(%f %f %f)", 
               scene.objects.back().bbox.min_pos.x, scene.objects.back().bbox.min_pos.y, scene.objects.back().bbox.min_pos.z,
               scene.objects.back().bbox.max_pos.x, scene.objects.back().bbox.max_pos.y, scene.objects.back().bbox.max_pos.z);
      }
      logerr("Conjuction %u END",id);
    }

    for (auto &c : scene.conjunctions)
    {
      printf("(%u %u)[", c.offset, c.offset+c.size);
      for (int i=0;i<c.size;i++)
      printf("%u ", scene.objects[c.offset+i].type);
      printf("]\n");
    }
    for (auto &v : scene.parameters)
      printf("%f, ", v);
    printf("\n");

    return scene;
  }
}