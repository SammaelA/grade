#include "SurfaceMap.h"
#include "BreadthSearch.h"

unsigned SurfaceMap::PointToPixel(BiomePoint point)
{
    static const unsigned WATER_LOW = (165u << 16) | (217u << 8) | (255u);
    static const unsigned WATER_HIGH = (0u << 16) | (145u << 8) | (255u);
    static const unsigned GROUND_LOW = (229u << 16) | (149u << 8) | (95u);
    static const unsigned GROUND_HIGH = (135u << 16) | (54u << 8) | (0u);

    unsigned color;
    switch (point.biome)
    {
        case (Biome::WATER):
        {
            float grade = point.intensity * 2 / std::max(m_width, m_height);
            grade = clamp(grade, 0.0f, 1.0f);
            color = ColorLerp(WATER_LOW, WATER_HIGH, grade);
            break;
        }
        case (Biome::GROUND):
        {
            float grade = point.intensity / 2 / std::max(m_width, m_height);
            grade = clamp(grade, 0.0f, 1.0f);
            color = ColorLerp(GROUND_LOW, GROUND_HIGH, grade);
            break;
        }
        default:
        {
            debug("WRONG SURFACE MAP!");
            throw std::exception{};
        }
    }
    return (color << 8) | 255;        
}

SurfaceMap::SurfaceMap(int w, int h) :
    CGenBiomeMap(w, h) {};

SurfaceMap* SurfaceMap::ScanBiomemap(const CGenBiomeMap& map)
{
    SurfaceMap* surface = new SurfaceMap(map.m_width, map.m_height);
    for (int i = -surface->m_width; i <= surface->m_width; i++)
    {
        for (int j = -surface->m_height; j <= surface->m_height; j++)
        {
            if (map.Get(vec2Int(i,j)).biome == Biome::WATER)
            {
                surface->Set(vec2Int(i,j), Biome::WATER);
            } 
            else if (map.Get(vec2Int(i,j)).biome <= Biome::GROUND)
            {
                surface->Set(vec2Int(i,j), Biome::GROUND);
            }          
        }
    }


    //FIND SHORELINE
    vec2Int justPoint{0, 0};
    std::function biomeFindNeighbs_FindAny
    {
        [surface](vec2Int coords) 
        {
            std::vector<Biome> biomes;
            for (int i = 0; i != static_cast<int>(Biome::ENUM_LEN); i++)
            {
                biomes.push_back(static_cast<Biome>(i));
            }
            return surface->biomeFindNeighbs_FindBiomeVec(coords, biomes);
        }
    };
    std::vector<vec2Int> shoreline;
    std::function biomeApplyToReachedPoints_SaveShoreline
    {
        [surface, &shoreline](vec2Int coords, float distance)
        {
            if (surface->Get(coords).biome == Biome::GROUND)
            {
                std::vector<vec2Int> neighbs;
                for (int i = 0; i < 4; i++)
                {
                    vec2Int current = coords + vec2Int{(i % 2) * (i - 2), ((i + 1) % 2) * (1 - i)};
                    if (surface->IsPointValid(current))
                    {
                        if (surface->Get(current).biome == Biome::WATER)
                            shoreline.push_back(current);
                    }
                }
            }
        }
    };

    biomeBreadthSearch(
        std::vector<vec2Int>{justPoint},
        biomeFindNeighbs_FindAny,
        biomeApplyToReachedPoints_SaveShoreline,
        surface->biomeFindDistance_FromOrigin
    );

    //CALIBRATE WATER
    std::function biomeFindNeighbs_FindWater
    {
        [surface](vec2Int coords) 
        {
            return surface->biomeFindNeighbs_FindBiome(coords, Biome::WATER);
        }
    };
    std::function biomeApplyToReachedPoints_SetWater_DistanceWise
    {
        [surface](vec2Int coords, float distance)
        {
            if (surface->Get(coords) == Biome::WATER)
            {
                surface->Set(coords, BiomePoint{Biome::WATER, distance});
            }
        }
    };
    biomeBreadthSearch(
        shoreline,
        biomeFindNeighbs_FindWater,
        biomeApplyToReachedPoints_SetWater_DistanceWise,
        surface->biomeFindDistance_FromOrigin
    );

    //CALIBRATE GROUND
    std::function biomeFindNeighbs_FindGround
    {
        [surface](vec2Int coords) 
        {
            return surface->biomeFindNeighbs_FindBiome(coords, Biome::GROUND);
        }
    };
    std::function biomeApplyToReachedPoints_SetGround_DistanceWise
    {
        [surface](vec2Int coords, float distance)
        {
            if (surface->Get(coords) == Biome::GROUND)
            {
                surface->Set(coords, BiomePoint{Biome::GROUND, distance});
            }
        }
    };
    biomeBreadthSearch(
        shoreline,
        biomeFindNeighbs_FindGround,
        biomeApplyToReachedPoints_SetGround_DistanceWise,
        surface->biomeFindDistance_FromOrigin
    );
    
    surface->SaveAsTexture();
    return surface;
}

SurfacePoint::SurfacePoint(Surface surf, float _intensity)
{
    surface = surf;
    intensity = _intensity;
}

SurfacePoint::SurfacePoint(const BiomePoint& bp)
{
    intensity = bp.intensity;
    if (bp.biome == Biome::GROUND)
        surface = Surface::EARTH;
    else if (bp.biome == Biome::WATER)
        surface = Surface::LIQUID;
    else
    {
        debug("INCORRECT TRANSFORMATION FROM BIOMEPOINT!");
        throw std::exception{};
    }
}

std::vector<SurfacePoint> SurfaceMap::GetSurfaceInterpolated_Normalized(glm::vec2 pos)
{
    std::vector<BiomePoint> temp = GetInterpolated_Normalized(pos);
    std::vector<SurfacePoint> res;
    for (BiomePoint bp: temp)
        res.push_back(bp);
    return res;
}

SurfacePoint SurfaceMap::GetSurface(vec2Int v)
{
    return Get(v);
}

float SurfaceMap::GetMaxIntensity(Surface surf)
{
    float maxSurfaceIntensity = 0.f;
    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            SurfacePoint sp = Get(vec2Int{i,j});
            if (sp.surface != surf)
                continue;
            float cur = Get(vec2Int{i,j}).intensity;
            if (cur > maxSurfaceIntensity)
                maxSurfaceIntensity = cur;
        }
    }
    return maxSurfaceIntensity;
}
