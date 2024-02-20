#include "sdf_reconstruction_common.h"
#include "sdf_node.h"
#include "optimization.h"
#include "tinyEngine/camera.h"
#include "tinyEngine/engine.h"
#include "common_utils/bbox.h"
#include "sdf_rendering.h"

namespace upg
{
  float normal_pdf(float x, float x0, float sigma)
  {
    return (1 / sqrt(2 * PI * sigma * sigma)) * exp(-SQR(x - x0) / (2 * sigma * sigma));
  }

  PointCloudReference get_point_cloud_reference(const Block &input_blk)
  {
    PointCloudReference reference;
    Block *synthetic_reference = input_blk.get_block("synthetic_reference"); //reference is parameters for our own generator
    Block *model_reference = input_blk.get_block("model_reference"); //reference is an .obj (or other) file with 3D model
    assert(!(synthetic_reference && model_reference));
    if (synthetic_reference)
    {
      reference.is_synthetic = true;
      synthetic_reference->get_arr("structure", reference.structure.s);
      synthetic_reference->get_arr("params", reference.parameters.p);
      ProceduralSdf sdf(reference.structure);
      sdf.set_parameters(reference.parameters.p);

      int points = synthetic_reference->get_int("points_count", 10000);
      sdf_to_point_cloud(sdf, points, &(reference.points), &(reference.outside_points));
      sdf_to_point_cloud_with_dist(sdf, 2*points, &(reference.d_points), &(reference.d_distances));

      CameraSettings camera;
      camera.origin = glm::vec3(0,0,3);
      camera.target = glm::vec3(0,0,0);
      camera.up = glm::vec3(0,1,0);
      Texture t = render_sdf(sdf, camera, 512, 512, 16);
      engine::textureManager->save_png(t, "reference_sdf");
    }
    else if (model_reference)
    {
      //TODO
    }

    return reference;
  };
}