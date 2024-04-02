#include "upg.h"
#include "tinyEngine/engine.h"
#include "sdf_scene_convert.h"
#include "sdf_octree.h"
#include "common_sdf_scenes.h"
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

  void perform_tests_sdf_scene(const std::vector<int> &test_ids)
  {
    std::vector<int> tests = test_ids;

    std::vector<std::function<void(void)>> test_functions = {
      scene_test_1
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