#pragma once

namespace voxelization
{
  float max(float x, float y);
  float min(float x, float y);
  float clamp(float x, float x_min, float x_max);

  void software_render_test_3d();
};