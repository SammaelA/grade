#pragma once
#include <glm/gtc/matrix_transform.hpp>
#include <glm/glm.hpp>

struct Camera
{
    glm::vec3 pos = glm::vec3(0, 10, 10);
    glm::vec3 front = glm::vec3(0, 0, -10);
    glm::vec3 up = glm::vec3(0, 1, 0);
    glm::mat4 camera_mat = glm::mat4(1.0f);
    float yaw = 0;
    float pitch = 0;
    float roll = 0;
    glm::mat4 &camera() { camera_mat = glm::lookAt(pos, pos + front, up); return camera_mat;}
};
struct DirectedLight
{
    glm::vec3 dir;
    float intensity;
    glm::vec3 color;
    float ambient_q;
    float diffuse_q;
    float specular_q;
    bool has_shadow_map;
    glm::vec2 shadow_map_size;
};