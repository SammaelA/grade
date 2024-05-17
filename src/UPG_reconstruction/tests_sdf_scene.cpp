#include "upg.h"
#include "tinyEngine/engine.h"
#include "sdf_scene_convert.h"
#include "sdf_octree.h"
#include "common_sdf_scenes.h"
#include "sdf_grid.h"
#include "sdf_bboxes.h"
#include <functional>



#include "LiteRT/utils/watertight_mesh.h"
#include "LiteScene/cmesh4.h"



namespace upg
{

  //TEST 1 SDF SCENE DISTANCE EVALUATION
  //checks if distances evaluated for ProceduralSdf and
  //SdfScene match each other
  void scene_test_1()
  {
    SceneDesc s = scene_chair();

    // SDF represented with tree of SdfNode
    ProceduralSdf reference_sdf(s.first);
    reference_sdf.set_parameters(s.second.p);

    // SDF represented in CGS manner, optimized for render
    SdfScene scene = create_sdf_scene(s.first, s.second);
    auto scene_sdf = get_SdfSceneFunction(scene);

    AABB box = reference_sdf.get_bbox();

    unsigned tries = 50000;
    unsigned errors = 0;
    for (int i=0;i<tries;i++)
    {
      float3 p = box.min_pos + box.size()*float3(urand(), urand(), urand());
      float d1 = reference_sdf.get_distance(p);
      float d2 = scene_sdf->eval_distance(p);
      if (abs(d1 - d2) > 1e-6)
      {
        errors++;
        logerr("wrong distance in (%f %f %f): d1=%f d2=%f", p.x, p.y, p.z, d1, d2);
      }
    }

    debug("TEST 1. SDF SCENE DISTANCE EVALUATION\n");
    debug("  1.1. %-64s", "Chair scene ");
    if (errors == 0)
      debug("passed\n");
    else
      debug("FAILED n_errors = %d\n", errors);
  }

  //TEST 2 SDF GRID
  //testing grid sdf, both the quality of grid representation
  //and that identity of representation with ISdfGridFunction
  //and SdfGridNode
  void scene_test_2()
  {
    auto scene = scene_chair();
    ProceduralSdf reference_sdf(scene.first);
    reference_sdf.set_parameters(scene.second.p);
    AABB bbox = AABB({-1,-1,-1},{1,1,1});

    unsigned grid_size = 256;
    std::vector<float> grid_data(grid_size*grid_size*grid_size, 0);

    GridSdfNode::primitive_SDF_to_grid(reference_sdf, bbox, grid_data.data(), grid_size);
    ProceduralSdf sdf({{SdfNodeType::GRID_256}});
    sdf.set_parameters(grid_data);

    auto grid_sdf = get_SdfGridFunction(SdfGridView(uint3(grid_size,grid_size,grid_size), grid_data));

    unsigned tries = 50000;
    unsigned errors = 0;
    double sum_diff = 0;
    for (int i=0;i<tries;i++)
    {
      float3 p = bbox.min_pos + bbox.size()*float3(urand(), urand(), urand());
      float d1 = reference_sdf.get_distance(p);
      float d2 = grid_sdf->eval_distance(p);//sdf.get_distance(p);
      float d3 = grid_sdf->eval_distance(p);
      if (abs(d2 - d3) > 1e-6)
      {
        errors++;
        logerr("wrong distance in (%f %f %f): d1=%f d2=%f d3=%f", p.x, p.y, p.z, d1, d2, d3);
      }
      sum_diff += abs(d1 - d3);
    }
    float diff = sum_diff/tries;

    debug("TEST 2. SDF GRID\n");
    debug("  2.1. %-64s", "Grid distance evaluation ");
    if (errors == 0)
      debug("passed\n");
    else
      debug("FAILED n_errors = %d\n", errors);
    debug("  2.2. %-64s", "Grid reconstruction ");
    if (diff <= 0.0001)
      debug("passed\n");
    else
      debug("FAILED diff = %f\n", diff);
  }

  //TEST 3 SDF OCTREE
  //testing grid sdf, both the quality of grid representation
  //and that identity of representation with ISdfGridFunction
  //and SdfGridNode
  void scene_test_3()
  {
    auto scene = scene_chair();
    ProceduralSdf reference_sdf(scene.first);
    reference_sdf.set_parameters(scene.second.p);
    AABB bbox = AABB({-1,-1,-1},{1,1,1});

    std::vector<float> data = {0};
    ProceduralSdf g_sdf({{SdfNodeType::OCTREE}});
    g_sdf.set_parameters(data);
    dynamic_cast<OctreeSdfNode*>(g_sdf.root)->construct( [&reference_sdf](const float3 &p) {
      /*logerr("sample %f %f %f",p.x,p.y,p.z);*/ return reference_sdf.get_distance(p); });

    SparseOctreeBuilder &octree = dynamic_cast<OctreeSdfNode*>(g_sdf.root)->octree;
    auto octree_sdf = get_SdfOctreeFunction({(unsigned)(octree.get_nodes().size()), octree.get_nodes().data()});

    unsigned tries = 50000;
    unsigned errors = 0;
    double sum_diff = 0;
    for (int i=0;i<tries;i++)
    {
      float3 p = bbox.min_pos + bbox.size()*float3(urand(), urand(), urand());
      float d1 = reference_sdf.get_distance(p);
      float d2 = g_sdf.get_distance(p);
      float d3 = octree_sdf->eval_distance(p);
      if (abs(d2 - d3) > 1e-6)
      {
        errors++;
        logerr("wrong distance in (%f %f %f): d1=%f d2=%f d3=%f", p.x, p.y, p.z, d1, d2, d3);
      }
      sum_diff += abs(d1 - d3);
    }
    float diff = sum_diff/tries;

    debug("TEST 3. SDF OCTREE\n");
    debug("  3.1. %-64s", "Octree distance evaluation ");
    if (errors == 0)
      debug("passed\n");
    else
      debug("FAILED n_errors = %d\n", errors);
    debug("  3.2. %-64s", "Octree reconstruction ");
    if (diff <= 0.01)
      debug("passed\n");
    else
      debug("FAILED diff = %f\n", diff);
  }

  void scene_test_4()
  {
    auto scene = scene_chair();
    ProceduralSdf reference_sdf(scene.first);
    reference_sdf.set_parameters(scene.second.p);
    AABB bbox = AABB({-1,-1,-1},{1,1,1});
    //reference_sdf.get_distance(float3(0,0,0));

    std::vector<AABB> bbox_list = get_bbox_list([&reference_sdf](const float3 &pos) { return reference_sdf.get_distance(pos); },
                                                bbox, 0, 3);
    
    unsigned tries = 50000;
    unsigned errors = 0;
    for (int i=0;i<tries;i++)
    {
      float3 p = bbox.min_pos + bbox.size()*float3(urand(), urand(), urand());
      float d_ref = reference_sdf.get_distance(p);
      if (d_ref > 0)
        continue;
      int box_id = -1;
      for (int j=0;j<bbox_list.size();j++)
      {
        if (bbox_list[j].contains(p))
        {
          if (box_id == -1)
            box_id = j;
          else
          {
            errors++;
            logerr("boxes %d and %d both contain (%f %f %f)", box_id, j, p.x, p.y, p.z);
          }
        }
      }
      if (box_id == -1)
      {
        errors++;
        logerr("(%f %f %f) does not belog to any box", p.x, p.y, p.z);
      }
    }

    debug("TEST 4. SDF BBOX LIST\n");
    debug("  4.1. %-64s", "Bbox list is correct ");
    if (errors == 0)
      debug("passed\n");
    else
      debug("FAILED n_errors = %d\n", errors);
  }


// -------------------------- watertight -----------------------------

  // Checking the teapot for watertight
  void scene_test_5(){
    debug("TEST 5. TEAPOT WATERTIGHT CHECK\n");
    auto mesh = cmesh4::LoadMeshFromVSGF("modules/LiteRT/scenes/01_simple_scenes/data/bunny.vsgf");
    
    bool res = watertight_mesh(mesh);
    if(res == 1)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 1 the cube for watertight
  void scene_test_6(){
    debug("TEST 6. CUBE WATERTIGHT CHECK\n");
    cmesh4::SimpleMesh mesh;

    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 1)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 2 the cube for watertight:
  // Instead of triangle (0, 4, 7), consider triangle (0, 7, 8) triangle (8, 4, 7)
  void scene_test_7(){
    debug("TEST 7. CUBE WATERTIGHT CHECK - TWO AND ONE EDGES\n");
    cmesh4::SimpleMesh mesh;
    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));
    mesh_vertices.push_back(float4(0.5, 0, 0, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 7, 8, 8, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 1)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 3 the cube for watertight: Add 3 equal triangles
  void scene_test_8(){
    debug("TEST 8. CUBE WATERTIGHT CHECK - REDUNDANT VERICES\n");
    cmesh4::SimpleMesh mesh;

    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0,
      2, 4, 7
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 1)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 4 the cube for watertight: Triangle in a triangle
  void scene_test_9(){
    debug("TEST 9. CUBE WATERTIGHT CHECK - TRIANGLE IN TRIANGLE\n");
    cmesh4::SimpleMesh mesh;
    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));

    mesh_vertices.push_back(float4(1, 0.7, 0.7, 8));
    mesh_vertices.push_back(float4(1, 0.3, 0.3, 8));
    mesh_vertices.push_back(float4(1, 0.3, 0.7, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0,
      8, 9, 10
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 1)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 5 the cube for watertight
  // Two equal elements in mesh_vertices that are pointed to by
  // different mesh_indices indexes
  void scene_test_10(){
    debug("TEST 10. EQUAL VERTICES\n");
    cmesh4::SimpleMesh mesh;
    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));
    mesh_vertices.push_back(float4(0, 0, 0, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 8, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 1)
      debug("passed\n");
    else
      debug("error\n");
  }



  // -------------------------- non- watertight -----------------------------


  // Checking 1 the cube for non-watertight: Indexes are not multiples of three
  void scene_test_11(){
    debug("TEST 11. CUBE - BROKEN INDEX BUFFER\n");
    cmesh4::SimpleMesh mesh;
    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 0)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 2 the cube for non-watertight:
  // Two equal vertices in the same triangle
  void scene_test_12(){
    debug("TEST 12. CUBE - ZERO AREA TRIANGLE\n");
    cmesh4::SimpleMesh mesh;
    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 6, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 0)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 3 the cube for non-watertight: There are edge(s) that occur only once
  void scene_test_13(){
    debug("TEST 13. CUBE WITH HOLE\n");
    cmesh4::SimpleMesh mesh;

    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 2, 2, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 8, 0,
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 0)
      debug("passed\n");
    else
      debug("error\n");
  }

  // Checking 4 the cube for non-watertight:
  // There are edge(s) that occur more than twice
  void scene_test_14(){
    debug("TEST 14. CUBE - NOT TWO-MANIFOLD\n");
    cmesh4::SimpleMesh mesh;
    std::vector<LiteMath::float4> mesh_vertices;

    mesh_vertices.push_back(float4(0, 0, 0, 8));
    mesh_vertices.push_back(float4(0, 0, 1, 8));
    mesh_vertices.push_back(float4(1, 1, 1, 8));
    mesh_vertices.push_back(float4(0, 1, 1, 8));
    mesh_vertices.push_back(float4(1, 0, 0, 8));
    mesh_vertices.push_back(float4(1, 1, 0, 8));
    mesh_vertices.push_back(float4(0, 1, 0, 8));
    mesh_vertices.push_back(float4(1, 0, 1, 8));
    mesh_vertices.push_back(float4(2, 0, 1, 8));

    mesh.vPos4f = mesh_vertices;

    std::vector<unsigned int> mesh_indices = {
      0, 4, 7, 1, 0, 7, 6, 2, 5, 3, 6, 2, 1, 6, 0, 3, 1, 6,
      4, 2, 5, 4, 2, 7, 3, 7, 1, 7, 3, 2, 5, 4, 0, 6, 5, 0,
      0, 4, 8, 0, 8, 4
    };

    mesh.indices = mesh_indices;

    bool res = watertight_mesh(mesh);
    if(res == 0)
      debug("passed\n");
    else
      debug("error\n");
  }

  void perform_tests_sdf_scene(const std::vector<int> &test_ids)
  {
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      scene_test_1, scene_test_2, scene_test_3, scene_test_4,

      scene_test_5, scene_test_6, scene_test_7, scene_test_8, scene_test_9,
      scene_test_10, scene_test_11, scene_test_12, scene_test_13, scene_test_14
    };

    if (tests.empty())
    {
      tests.resize(test_functions.size());
      for (int i=0;i<test_functions.size();i++)
        tests[i] = i+1;
    }

    for (int i=0;i<80;i++)
      debug("#");
    debug("\nSDF SCENE TESTS\n");
    for (int i=0;i<80;i++)
      debug("#");
    debug("\n");
    
    for (int i : tests)
    {
      assert(i > 0 && i <= test_functions.size());
      test_functions[i-1]();
    }
  }
}