#pragma once
#include "cities_generator/global.h"
#include "RenderableAsMap.h"
#include "third_party/stb_image.h"

template<class T>
class DefaultMap : RenderableAsMap
{
    friend class Landscape;
    friend class SurfaceMap;

    protected:
        T* m_data;
        int m_width, m_height;
        bool m_isValidBorders;
        int dataSize;

        DefaultMap(int wigth, int height);
        DefaultMap() = default;
        DefaultMap(const DefaultMap& m);
        virtual ~DefaultMap();
        T Get(int x, int y, bool safe = false, T baseValue = T{}) const;
        T Get(vec2Int pos, bool safe = false, T baseValue = T{}) const;
        virtual void Set(int x, int y, T val);
        virtual void Set(vec2Int pos, T val);
        virtual bool IsTextureReady();
        virtual void SaveAsTexture();
        bool IsPointValid(int x, int y);
        bool IsPointValid(vec2Int v);
        bool IsPointValid_Normalized(glm::vec2 v);
        virtual bool IsPointValid_IntContinuous(glm::vec2 v);
        virtual unsigned PointToPixel(T val);
        virtual glm::vec2 PointIntToNormalized(vec2Int, bool = false) {};
        virtual glm::vec2 PointIntContinuousToNormalized(glm::vec2, bool = false) {debug("NO METHOD");throw std::exception{};};
        virtual glm::vec2 PointNormalizedToIntContinuous(glm::vec2, bool = false) {debug("NO METHOD");throw std::exception{};};
        // virtual vec2Int PointNormalizedToInt(glm::vec2 pos) = 0;
        unsigned ColorLerp(unsigned col1, unsigned col2, float t);
    
    public:
        virtual GLuint GetTexture();
        glm::vec2 GetSize();
};

#include "DefaultMap.tpp"