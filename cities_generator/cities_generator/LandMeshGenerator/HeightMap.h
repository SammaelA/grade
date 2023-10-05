#pragma once
#include "DefaultMap.h"
#include "cities_generator/global.h"
#include <array>

class HeightMap : public DefaultMap<float>
{
    friend class Landscape;

    static constexpr float DENSITY_PER_REAL = 0.8f; //FINAL 0.8

    float m_maxHeight, m_minHeight;
    bool m_isValidBorders;

    HeightMap(int wigth, int height);
    HeightMap(const HeightMap&);
    virtual ~HeightMap();
    float GetSafe_EdgeStrectched(int x, int y);
    float GetSafe_EdgeStrectched(vec2Int v);
    float GetSafe_EdgeRepeated(int x, int y);
    float GetSafe_EdgeMirrored(vec2Int v);
    float GetSafe_EdgeMirrored(int x, int y);
    float GetInterpolated_IntContinuous(glm::vec2 coords);
    virtual void Set(int x, int y, float val);
    virtual void Set(vec2Int, float val);
    void FillByPerlin(float min, float max, glm::ivec2 shift = glm::ivec2(0,0));
    void UpdateHeightBorders(float height);
    void RecalculateHeightBorders();
    void operator-=(const HeightMap& other);
    void operator+=(const HeightMap& other);
    void ComponentWiseOperation(const HeightMap& other, const char operation);
    void BlurHeightmap();
    virtual unsigned PointToPixel(float height);
    virtual glm::vec2 PointIntToNormalized(vec2Int pos, bool safe = false);
    virtual glm::vec2 PointIntContinuousToNormalized(glm::vec2 pos, bool safe = false);
    virtual glm::vec2 PointNormalizedToIntContinuous(glm::vec2 pos, bool safe = false);
    // virtual vec2Int PointNormalizedToInt(glm::vec2 pos);
    
    static std::array<std::array<int, 2>, 2> GetSurroundingRect(glm::vec2 TC);

    public:
};