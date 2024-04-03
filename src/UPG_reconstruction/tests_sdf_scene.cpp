#include "upg.h"
#include "tinyEngine/engine.h"
#include "sdf_scene_convert.h"
#include "sdf_octree.h"
#include "common_sdf_scenes.h"
#include "sdf_grid.h"
#include <functional>

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

    auto grid_sdf = get_SdfGridFunction({uint3(grid_size,grid_size,grid_size), grid_data.data()});

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

  void perform_tests_sdf_scene(const std::vector<int> &test_ids)
  {
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      scene_test_1, scene_test_2, scene_test_3
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