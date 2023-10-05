#pragma once
#include "DefaultMap.h"
#include "BiomeMap.h"
#include "cities_generator/global.h"

enum Surface
{
    LIQUID = 0,
    EARTH = 1,
};

class SurfacePoint
{
    public:
        Surface surface;
        float intensity;
        SurfacePoint(Surface surf,float intensity);
        SurfacePoint(const BiomePoint& bp);
};

class SurfaceMap : public CGenBiomeMap
{
    virtual unsigned PointToPixel(BiomePoint point);
    SurfaceMap(int w, int h);

    private:
        using CGenBiomeMap::GetInterpolated_Normalized;
        using DefaultMap::Get;
        using DefaultMap::Set;
        
    public:
        static SurfaceMap* ScanBiomemap(const CGenBiomeMap& map);
        std::vector<SurfacePoint> GetSurfaceInterpolated_Normalized(glm::vec2 pos);
        SurfacePoint GetSurface(vec2Int v);
        float GetMaxIntensity(Surface s);
};