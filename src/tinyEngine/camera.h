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
    glm::vec3 dir = glm::vec3(0,1,0);
    float intensity = 0;
    glm::vec3 color = glm::vec3(1,1,1);
    float ambient_q = 1;
    float diffuse_q = 0;
    float specular_q = 0;
    bool has_shadow_map = false;
    glm::vec2 shadow_map_size = glm::vec2(4096, 4096);
};