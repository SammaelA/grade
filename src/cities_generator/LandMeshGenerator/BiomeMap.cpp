#include "BiomeMap.h"
#include "BreadthSearch.h"
#include <unordered_map>
#include <vector>
#include <memory>

bool operator<= (const Biome a, const Biome b)
{
    if (a == b)
        return true;
    if (a == Biome::PLAIN)
    {
        if (b == Biome::GROUND)
            return true;
    }
    if (a == Biome::MOUNTAINS)
    {
        if (b == Biome::GROUND)
            return true;
    }
    if (a == Biome::CITY)
    {
        if (b == Biome::GROUND)
            return true;
    }
    if (a == Biome::SHORELINE)
    {
        if (b == Biome::GROUND)
            return true;
    }
    return false;
}

BiomePoint::BiomePoint(Biome b, float f)
{
    biome = b;
    intensity = f;
}

BiomePoint::operator Biome() const
{
    return biome;
}

CGenBiomeMap::CGenBiomeMap(int width, int height) :
    DefaultMap(width, height) 
{}

CGenBiomeMap::~CGenBiomeMap() {}

unsigned CGenBiomeMap::PointToPixel(BiomePoint point)
{
    static const unsigned WATER_LOW = (165u << 16) | (217u << 8) | (255u);
    static const unsigned WATER_HIGH = (0u << 16) | (145u << 8) | (255u);
    static const unsigned GROUND_LOW = (229u << 16) | (149u << 8) | (95u);
    static const unsigned GROUND_HIGH = (135u << 16) | (54u << 8) | (0u);
    static const unsigned PLAIN_LOW = (79u << 16) | (199u << 8) | (78u);
    static const unsigned PLAIN_HIGH = (20u << 16) | (63u << 8) | (20u);
    static const unsigned MOUNTAIN_LOW = (190u << 16) | (175u << 8) | (168u);
    static const unsigned MOUNTAIN_HIGH = (72u << 16) | (59u << 8) | (54u);
    static const unsigned CITY_LOW = (255 << 16) | (217 << 8) | (194);
    static const unsigned CITY_HIGH = (110 << 16) | (93 << 8) | (83);


    unsigned color;
    switch (point.biome)
    {
        case (Biome::WATER):
        {
            // float grade = point.intensity * 2 / std::max(m_width, m_height);
            // grade = clamp(grade, 0.0f, 1.0f);
            color = ColorLerp(WATER_LOW, WATER_HIGH, 0.5f);
            break;
        }
        case (Biome::GROUND):
        {
            // float grade = point.intensity / 2 / std::max(m_width, m_height);
            // grade = clamp(grade, 0.0f, 1.0f);
            color = ColorLerp(GROUND_LOW, GROUND_HIGH, 0.5f);
            break;
        }
        case (Biome::MOUNTAINS):
        {
            float grade = point.intensity * 2.0 / std::max(m_width, m_height);
            grade = clamp(grade, 0.0f, 1.0f);
            color = ColorLerp(MOUNTAIN_LOW, MOUNTAIN_HIGH, grade);
            break;
        }
        case (Biome::PLAIN):
        {
            // float grade = point.intensity / std::max(m_width, m_height);
            // grade = clamp(grade, 0.0f, 1.0f);
            // color = ColorLerp(PLAIN_LOW, PLAIN_HIGH, grade);
            color = ColorLerp(PLAIN_LOW, PLAIN_HIGH, 0.5f);
            break;
        }
        case (Biome::CITY):
        {
            float grade = point.intensity / std::max(m_width, m_height);
            grade = clamp(grade, 0.0f, 1.0f);
            color = ColorLerp(CITY_LOW, CITY_HIGH, grade);
            break;
        }
        case (Biome::SHORELINE):
            color = (252u << 16) | (221u << 8) | (118u);
            break;

        case (Biome::EMPTY):
            color = (150u << 16) | (150u << 8) | (150);
            break;
    }
    return (color << 8) | 255;        
}

vec2Int CGenBiomeMap::PointNormalizedToInt(glm::vec2 pos)
{
    if (pos.x < 0 || pos.y < 0 || pos.x > 1 || pos.y > 1)
    {
        debug("BAD INPUT! [BiomeMap_PointNormalizedToInt]", pos.x, pos.y);
        throw std::exception{};
    }
    pos.x = clamp(pos.x, 0.0f, 0.999999f);
    pos.y = clamp(pos.y, 0.0f, 0.999999f);
    vec2Int res;
    res.x = int(floorf(pos.x * (2 * m_width + 1) - m_width));
    res.y = int(floorf(pos.y * (2 * m_height + 1) - m_height));
    return res;
}

bool CGenBiomeMap::IsPointValid_IntContinuous(glm::vec2 v)
{
    float x = v.x;
    float y = v.y;
    return (x >= -m_width - 0.5f && x <= m_width + 0.5f && y >= -m_height - 0.5f && y <= m_height + 0.5f);
}

glm::vec2 CGenBiomeMap::PointNormalizedToIntContinuous(glm::vec2 pos, bool safe)
{
    if (pos.x < 0 || pos.y < 0 || pos.x > 1 || pos.y > 1)
    {
        if (!safe)
        {
            debug("BAD INPUT! [BiomeMap_PointNormalizedToIntContinuous]", pos.x, pos.y);
            throw std::exception{};
        }
    }
    glm::vec2 res;
    res.x = pos.x * (2 * m_width + 1) - m_width - 0.5f;
    res.y = pos.y * (2 * m_height + 1) - m_height - 0.5f;
    return res;
}

std::vector<BiomePoint> CGenBiomeMap::GetInterpolated_Normalized(glm::vec2 pos)
{
    std::vector<BiomePoint> points;
    std::array<vec2Int, 4> surroundingRect;
    vec2Int coords = PointNormalizedToInt(pos);
    glm::vec2 gridPos = PointIntToNormalized(coords);
    vec2Int dir = vec2Int(sign_no_zero(pos.x - gridPos.x), sign_no_zero(pos.y - gridPos.y));
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            surroundingRect[i * 2 + j] = coords + vec2Int(dir.x * i, dir.y * j);
        }
    }
    
    glm::vec2 cellSize = PointIntToNormalized(vec2Int(1,1)) - PointIntToNormalized(vec2Int(0,0));
    cellSize = glm::abs(cellSize);
    for (int i = 0; i < 4; i++)
    {
        if (!IsPointValid(surroundingRect[i]))
            continue;
        BiomePoint currentPoint = Get(surroundingRect[i]);
        glm::vec2 currentPointPos = PointIntToNormalized(surroundingRect[i]);
        
        int j = 0;
        for (; j < points.size(); j++)
        {
            if (points[j].biome == currentPoint.biome)
                break;
        }
        if (j == points.size())
            points.push_back(currentPoint.biome);

        float intensity = currentPoint.intensity * 
                          (1 - fabsf(pos.x - currentPointPos.x) / cellSize.x) *
                          (1 - fabsf(pos.y - currentPointPos.y) / cellSize.y);  
        points[j].intensity += intensity;
    }

    return points;
}

std::vector<vec2Int> CGenBiomeMap::biomeFindNeighbs_FindBiome(vec2Int coords, Biome BIOME)
{
    return biomeFindNeighbs_FindBiomeVec(coords, std::vector<Biome> {BIOME});
}

std::vector<vec2Int> CGenBiomeMap::biomeFindNeighbs_FindBiomeVec(vec2Int coords, std::vector<Biome> biomes)
{
    std::vector<vec2Int> neighbs;
    for (int i = 0; i < 4; i++)
    {
        vec2Int current = coords + vec2Int{(i % 2) * (i - 2), ((i + 1) % 2) * (1 - i)};
        if (IsPointValid(current))
        {
            BiomePoint bp = Get(current);
            for (Biome b: biomes)
            {
                if (Get(current) == b)
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
                BiomePoint bp = Get(current);
                for (Biome b: biomes)
                {
                    if (bp == b)
                    {
                        if (Get(coords + vec2Int{0, j}) == b ||
                            Get(coords + vec2Int{i, 0}) == b)
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

float CGenBiomeMap::biomeFindDistance_FromOrigin(vec2Int current, vec2Int previous, float prevDist, vec2Int origin)
{
    return (current - origin).length();
};


glm::vec2 CGenBiomeMap::PointIntToNormalized(vec2Int pos, bool safe)
{
    return PointIntContinuousToNormalized(glm::vec2(pos), safe);
}

glm::vec2 CGenBiomeMap::PointIntContinuousToNormalized(glm::vec2 pos, bool safe)
{
    if (!IsPointValid_IntContinuous(pos) && !safe)
    {
        debug("BAD INPUT! [PointIntContinuousToNormalized]");
        throw std::exception{};
    }
    return glm::vec2 {
        (float) (pos.x + m_width + 0.5) / (2 * m_width + 1),
        (float) (pos.y + m_height + 0.5) / (2 * m_height + 1)
    };
}

void CGenBiomeMap::CalibrateBiomemap()
{
    std::unordered_map<Biome, std::vector<vec2Int>> areas;
    for (int i = 0; i < (int)Biome::ENUM_LEN; i++)
    {
        Biome b = (Biome)i;
        areas.emplace(b, std::vector<vec2Int>{});
    }

    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            vec2Int coords{i, j};
            BiomePoint bp = Get(coords);
            areas[bp.biome].push_back(coords);
            Set(coords, BiomePoint{bp.biome, 1000000});
        }
    }

    for (int i = 0; i < (int)Biome::ENUM_LEN; i++)
    {
        Biome biome = (Biome)i;

        if (areas[biome].size() == 0)
            continue;

        std::function biomeFindNeighbs_DifferentFromCurrent
        {
            [this, biome]
            (vec2Int coords) 
            {
                std::vector<vec2Int> result;
                for (int i = 0; i < (int)Biome::ENUM_LEN; i++)
                {
                    Biome b = (Biome)i;
                    if (b == biome)
                        continue;

                    auto cur = biomeFindNeighbs_FindBiome(coords, b);
                    result.insert(result.end(), cur.begin(), cur.end());
                }
                return result;
            }
        };
        
        std::function biomeApplyToReachedPoints_UpdateIntensity
        {
            [this, biome](vec2Int coords, float distance)
            {
                BiomePoint bp = Get(coords);

                if (bp.biome == biome)
                    return;

                bp.intensity = std::min(distance, bp.intensity);
                Set(coords, bp);
            }
        };
        
        biomeBreadthSearch(
            areas[biome],
            biomeFindNeighbs_DifferentFromCurrent,
            biomeApplyToReachedPoints_UpdateIntensity,
            biomeFindDistance_FromOrigin
        );
    }

    SaveAsTexture();
}