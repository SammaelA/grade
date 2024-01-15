#pragma once
#include "public_camera.h"
#include "tinyEngine/camera.h"
namespace nsdf
{
  Camera convert(const CameraSettings &cam)
  {
    Camera c{cam.origin.x, cam.origin.y, cam.origin.z,
             cam.target.x, cam.target.y, cam.target.z,
             cam.up.x, cam.up.y, cam.up.z,
             cam.fov_rad,
             cam.z_near,
             cam.z_far};
    return c;
  }
  CameraSettings convert(const Camera &cam)
  {
    CameraSettings c{ glm::vec3(cam.pos_x, cam.pos_y, cam.pos_z),
                      glm::vec3(cam.target_x, cam.target_y, cam.target_z),
                      glm::vec3(cam.up_x, cam.up_y, cam.up_z),
                      cam.fov_rad,
                      cam.z_near,
                      cam.z_far};
    return c;
  }
}