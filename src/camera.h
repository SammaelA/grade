#pragma once
#include <glm/gtc/matrix_transform.hpp>
#include <glm/glm.hpp>

struct Camera
{
    glm::vec3 pos = glm::vec3(0, 10, 10);
    glm::vec3 front = glm::vec3(0, 0, -10);
    glm::vec3 up = glm::vec3(0,1,0);
    float yaw = 0;
    float pitch = 0;
    float roll = 0;
    glm::mat4 camera() {return glm::lookAt(pos, pos + front, up);}
};