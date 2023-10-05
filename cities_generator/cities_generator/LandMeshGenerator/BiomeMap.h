#pragma once
#include "DefaultMap.h"
#include "cities_generator/global.h"
#include <array>

enum class Biome
{
    WATER = 0,
    GROUND = 1,
    SHORELINE = 2,
    PLAIN = 3,
    CITY = 4,
    MOUNTAINS = 5,
    EMPTY = 6,
    ENUM_LEN = 7
};

bool operator<= (Biome a, Biome b);

struct BiomePoint
{
    Biome biome;
    float intensity;
    BiomePoint(Biome biome = Biome::EMPTY, float intensity = 0);
    operator Biome() const;
};

class CGenBiomeMap : public DefaultMap<BiomePoint>
{
    friend class Landscape;
    protected:
        CGenBiomeMap(int width, int height);
        virtual ~CGenBiomeMap();
        // virtual void Set(int x, int y, float val);
        virtual unsigned PointToPixel(BiomePoint point);
        virtual glm::vec2 PointIntToNormalized(vec2Int point, bool safe = false);
        virtual glm::vec2 PointIntContinuousToNormalized(glm::vec2 pos, bool safe = false);
        virtual glm::vec2 PointNormalizedToIntContinuous(glm::vec2 pos, bool safe = false);
        virtual bool IsPointValid_IntContinuous(glm::vec2 v);
        vec2Int PointNormalizedToInt(glm::vec2 pos);
        std::vector<BiomePoint> GetInterpolated_Normalized(glm::vec2 pos);

        //Breadth search functions
        std::vector<vec2Int> biomeFindNeighbs_FindBiome(vec2Int coords, Biome BIOME);
        std::vector<vec2Int> biomeFindNeighbs_FindBiomeVec(vec2Int coords, std::vector<Biome> biomes);
        static float biomeFindDistance_FromOrigin(vec2Int current, vec2Int previous, float prevDist, vec2Int origin);
        
        void CalibrateBiomemap();

    public:
        static constexpr float DENSITY_PER_REAL = 0.33f;
};