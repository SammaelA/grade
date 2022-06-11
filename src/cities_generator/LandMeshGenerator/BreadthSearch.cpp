#include "BreadthSearch.h"
#include <memory>
#include <unordered_map>
#include <queue>

void biomeBreadthSearch(
    const std::vector<vec2Int>& startPoints,
    std::function<std::vector<vec2Int>(vec2Int)> neighbsFunc,
    std::function<void(vec2Int, float)> applyFunc,
    std::function<float(vec2Int, vec2Int, float, vec2Int)> distanceFunc)
{
    struct graphNode
    {
        vec2Int position;
        float distance;
        vec2Int origin;
        graphNode(vec2Int p, float dist, vec2Int _origin) 
        {
            position = p;
            distance = dist;
            origin = _origin;
        }
    };
    auto cmp = [](const graphNode& a, const graphNode& b) 
    { 
        return a.distance > b.distance; 
    };
    std::priority_queue<graphNode, std::vector<graphNode>, decltype(cmp)> border(cmp);
    std::unordered_set<vec2Int> reachedPoints{};
    for (vec2Int point: startPoints)
        border.emplace(point, 0, point);
        
    while (border.size() > 0)
    {
        graphNode current = border.top();

        border.pop();
        if (reachedPoints.find(current.position) != reachedPoints.end())
            continue;

        applyFunc(current.position, current.distance);
        reachedPoints.emplace(current.position);

        std::vector<vec2Int> neighbs = neighbsFunc(current.position);
        for (vec2Int n: neighbs)
        {
            if (reachedPoints.find(n) != reachedPoints.end())
                continue;

            border.emplace(
                n, 
                distanceFunc(n, current.position, current.distance, current.origin),
                current.origin
            );
        }
    }
}
