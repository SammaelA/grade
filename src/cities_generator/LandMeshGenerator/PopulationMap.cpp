#include "PopulationMap.h"

PopulationPoint::PopulationPoint(Population p, float f)
{
    population = p;
    intensity = f;
}

PopulationMap::PopulationMap(CGenBiomeMap& biomemap) :
    DefaultMap(biomemap.GetSize().x, biomemap.GetSize().y) 
{}

PopulationMap::~PopulationMap() {}

unsigned PopulationMap::PointToPixel(PopulationPoint point)
{
    static const unsigned POPULATION_LOW = (227u << 16) | (202u << 8) | (189u);
    static const unsigned POPULATION_HIGH = (142u << 16) | (88u << 8) | (59u);
    static const unsigned UNPOPULATED = (246u << 16) | (243u << 8) | (225u);
    static const unsigned UNUSABLE = (180u << 16) | (180u << 8) | (180u);
    static const unsigned UNSCANNED = (208u << 16) | (75u << 8) | (230u);

    static const float POPULATED_COLOR_STEP = 0.1f;

    unsigned color = 0;
    switch (point.population)
    {
        case (Population::UNSCANNED):
        {
            color = UNSCANNED;
            break;
        }
        case (Population::EMPTY):
        {
            color = UNPOPULATED;
            break;
        }
        case (Population::UNUSABLE):
        {
            color = UNUSABLE;
            break;
        }
        case (Population::POPULATED):
        {
            float drawnIntensity = point.intensity;
            drawnIntensity = floorf(drawnIntensity / POPULATED_COLOR_STEP) * POPULATED_COLOR_STEP;
            color = ColorLerp(POPULATION_LOW, POPULATION_HIGH, drawnIntensity);
            break;
        }
    }
    return (color << 8) | 255;        
}

glm::vec2 PopulationMap::PointIntToNormalized(vec2Int pos, bool safe)
{
    if (!IsPointValid(pos.x, pos.y) && !safe)
    {
        debug("BAD INPUT! [PointIntToNormalized]");
        throw std::exception{};
    }
    return glm::vec2 {
        (float) (pos.x + m_width + 0.5) / (2 * m_width + 1),
        (float) (pos.y + m_height + 0.5) / (2 * m_height + 1)
    };
}

vec2Int PopulationMap::PointNormalizedToInt(glm::vec2 pos)
{
    if (pos.x < 0 || pos.y < 0 || pos.x > 1 || pos.y > 1)
    {
        debug("BAD INPUT! [PointNormalizedToInt]", pos.x, pos.y);
        throw std::exception{};
    }
    pos.x = clamp(pos.x, 0.0f, 0.999999f);
    pos.y = clamp(pos.y, 0.0f, 0.999999f);
    vec2Int res;
    res.x = int(floorf(pos.x * (2 * m_width + 1) - m_width));
    res.y = int(floorf(pos.y * (2 * m_height + 1) - m_height));
    return res;
}

std::vector<vec2Int> PopulationMap::populationFindNeighbs_FindPopulationVec(
        vec2Int coords, 
        std::vector<Population> populations
    )
{
    std::vector<vec2Int> neighbs;
    for (int i = 0; i < 4; i++)
    {
        vec2Int current = coords + vec2Int{(i % 2) * (i - 2), ((i + 1) % 2) * (1 - i)};
        if (IsPointValid(current))
        {
            PopulationPoint pp = Get(current);
            for (Population p: populations)
            {
                if (Get(current).population == p)
                {
                    neighbs.push_back(current);
                }
            }
        }
    }
    for (int i = -1; i <= 1; i+=2)
    {
        for (int j = -1; j <= 1; j+=2)
        {
            vec2Int current = coords + vec2Int{i, j};
            if (IsPointValid(current))
            {
                PopulationPoint pp = Get(current);
                for (Population p: populations)
                {
                    if (pp.population == p)
                    {
                        if (Get(coords + vec2Int{0, j}).population == p ||
                            Get(coords + vec2Int{i, 0}).population == p)
                        {
                            neighbs.push_back(current);
                        }
                    }
                }
            }
        }
    }
    return neighbs;
}

std::vector<vec2Int> PopulationMap::populationFindNeighbs_FindPopulation(
        vec2Int coords, 
        Population population
    )
{
    std::vector<Population> p{population};
    return populationFindNeighbs_FindPopulationVec(coords, p);
}

PopulationMap::PopulationMap(const PopulationMap& p)
    :DefaultMap(p)
{
    realCenterPos = p.realCenterPos;
    shoreCenterPos = p.shoreCenterPos;
}