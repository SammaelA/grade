#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

class CHTurtle
{
public:
    glm::vec3 dir = glm::vec3(0,0,1);
    glm::vec3 pos = glm::vec3(0,0,0);
    glm::vec3 right = glm::vec3(1,0,0);
    float width = 0;
    CHTurtle() {};
    CHTurtle(const CHTurtle &other) = default;

    void turn_right(float angle)
    {
        glm::vec3 axis = glm::normalize(glm::cross(dir, right));
        auto rot_quat = glm::angleAxis(glm::radians(angle), axis); 

        dir = glm::normalize(glm::rotate(rot_quat, dir));
        right = glm::normalize(glm::rotate(rot_quat, right));
    }
    void turn_left(float angle)
    {
        glm::vec3 axis = glm::normalize(glm::cross(dir, right));
        auto rot_quat = glm::angleAxis(glm::radians(-angle), axis); 

        dir = glm::normalize(glm::rotate(rot_quat, dir));
        right = glm::normalize(glm::rotate(rot_quat, right));
    }
    void pitch_up(float angle)
    {
        dir = glm::normalize(glm::rotate(glm::angleAxis(glm::radians(angle), right), dir));
    }
    void pitch_down(float angle)
    {
        dir = glm::normalize(glm::rotate(glm::angleAxis(glm::radians(-angle), right), dir));
    }
    void roll_right(float angle)
    {
        right = glm::normalize(glm::rotate(glm::angleAxis(glm::radians(angle), dir), right));
    }
    void roll_left(float angle)
    {
        right = glm::normalize(glm::rotate(glm::angleAxis(glm::radians(-angle), dir), right));
    }
    void move(float distance)
    {
        pos += distance*dir;
    }
    void set_width(float width)
    {
        this->width = width;
    }
};