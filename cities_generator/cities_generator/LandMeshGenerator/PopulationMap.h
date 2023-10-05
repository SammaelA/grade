#pragma once
#include "DefaultMap.h"
#include "BiomeMap.h"
#include "cities_generator/global.h"
#include <array>

enum class Population
{
    UNSCANNED = 0,
    UNUSABLE = 1,
    EMPTY = 2,
    POPULATED = 3
};

struct PopulationPoint
{
    Population population;
    float intensity;
    PopulationPoint(Population population = Population::UNSCANNED, float intensity = 0);
};

class PopulationMap : public DefaultMap<PopulationPoint>
{
    friend class Landscape;
    protected:
        vec2Int realCenterPos, shoreCenterPos;

        PopulationMap(CGenBiomeMap& biomemap);
        PopulationMap(const PopulationMap& p);
        virtual ~PopulationMap();
        // virtual void Set(int x, int y, float val);
        virtual unsigned PointToPixel(PopulationPoint point);
        virtual glm::vec2 PointIntToNormalized(vec2Int point, bool safe = false);
        vec2Int PointNormalizedToInt(glm::vec2 pos);
        std::vector<vec2Int> populationFindNeighbs_FindPopulationVec(vec2Int, std::vector<Population>);
        std::vector<vec2Int> populationFindNeighbs_FindPopulation(vec2Int, Population);
};