#pragma once
#include "common_utils/LiteMath_ext.h"
#include "common_utils/quaternion.h"

class CHTurtle
{
public:
    float3 dir = float3(0,0,1);
    float3 pos = float3(0,0,0);
    float3 right = float3(1,0,0);
    float width = 0;
    CHTurtle() {};

    void turn_right(float angle)
    {
        float3 axis = normalize(cross(dir, right));
        auto rot_quat = LiteMath::angleAxis(LiteMath::to_radians(angle), axis); 

        dir = normalize(LiteMath::rotate(rot_quat, dir));
        right = normalize(LiteMath::rotate(rot_quat, right));
    }
    void turn_left(float angle)
    {
        float3 axis = normalize(cross(dir, right));
        auto rot_quat = LiteMath::angleAxis(LiteMath::to_radians(-angle), axis); 

        dir = normalize(LiteMath::rotate(rot_quat, dir));
        right = normalize(LiteMath::rotate(rot_quat, right));
    }
    void pitch_up(float angle)
    {
        dir = normalize(LiteMath::rotate(LiteMath::angleAxis(LiteMath::to_radians(angle), right), dir));
    }
    void pitch_down(float angle)
    {
        dir = normalize(LiteMath::rotate(LiteMath::angleAxis(LiteMath::to_radians(-angle), right), dir));
    }
    void roll_right(float angle)
    {
        right = normalize(LiteMath::rotate(LiteMath::angleAxis(LiteMath::to_radians(angle), dir), right));
    }
    void roll_left(float angle)
    {
        right = normalize(LiteMath::rotate(LiteMath::angleAxis(LiteMath::to_radians(-angle), dir), right));
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