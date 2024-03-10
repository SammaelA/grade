#pragma once
#include "common_utils/LiteMath_ext.h"

struct EnvironmentParameters
{
    int day = 1;
    int month = 1;
    int year = 2021;

    float hours = 12;
    float minutes = 0;
    float seconds = 0;

    float latitude_deg = 60;
    float longitude_deg = 15;
};
class Sun
{
public:
    static glm::vec3 sun_direction(EnvironmentParameters &params);
};
