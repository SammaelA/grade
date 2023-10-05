#pragma once
#include "cities_generator/global.h"
#include <functional>

void biomeBreadthSearch(
    const std::vector<vec2Int>& startPoints,
    std::function<std::vector<vec2Int>(vec2Int)> neighbsFunc,
    std::function<void(vec2Int, float)> applyFunc,
    std::function<float(vec2Int, vec2Int, float, vec2Int)> distanceFunc
);