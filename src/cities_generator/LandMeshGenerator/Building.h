#pragma once
#include "cities_generator/global.h"
#include "../ShaderDataHolders/Archetypes/BasicGLTFModel.h"
#include "../ShaderDataHolders/RenderableObjectDataHolder.h"
#include <memory>
#include <map>

struct Rect;

struct Polygon
{
    std::vector<glm::vec2> points;
    bool IsPointInside(glm::vec2 p);
    glm::vec2 GetPoint(int index);
    Polygon() = default;
    explicit Polygon(Rect);
    bool CheckRectIntersect(Rect& r);
    std::unique_ptr<Rect> Insribe();
    Rect Outscribe();
};

struct Rect : public Polygon
{
    glm::vec2 RandomPoint(float shift);
    glm::vec2 center;
    glm::vec2 size;
    float angleRad; //counter-clockwise
    Rect(glm::vec2 p, glm::vec2 s, float a);
    std::vector<glm::vec2> GetPoints();
    float Area();
    
    protected:
        using Polygon::points;
};

struct Basement : public Rect
{
    float height;
    Basement() = delete;
    Basement(const Rect& p);
};

enum class BUILDING_TYPE
{
    SCYSCRAPER,
    HIGH_RISE
};

struct BuildingInfo
{
    BuildingInfo() = delete;
    BuildingInfo(std::string path);

    std::string path, roofPath, specularPath;
    float realHeight, roofRealSize;

    glm::vec2 GetRealSize();

    private:
        float WidthToHeight;
};

class Building : public BasicGLTFModel
{
    int sunPos_Uniform = MISSING_UNIFORM_POS;
    int sunColor_Uniform = MISSING_UNIFORM_POS;
    int cameraPos_Uniform = MISSING_UNIFORM_POS;

    int buildingTex_Uniform = MISSING_UNIFORM_POS;
    int buildingSpecTex_Uniform = MISSING_UNIFORM_POS;
    int roofTex_Uniform = MISSING_UNIFORM_POS;
    // int asphaltTexScale_Uniform = MISSING_UNIFORM_POS;
    // int markingTex_Uniform = MISSING_UNIFORM_POS;
    // int crossroadMarkingTex_Uniform = MISSING_UNIFORM_POS;

    GLuint buildingTex, roofTex, buildingSpecTex;
    // GLfloat asphaltTexScale = 0.2;
    // GLuint markingTex, crossroadMarkingTex;
    float elevation = 0.f;
    float height = 1.f;
    
    public:
        static std::map<BUILDING_TYPE, BuildingInfo> buildingTypes;
        static void InitBuildingTypes();
        
        Basement basement;
        Building(Basement, BUILDING_TYPE type, float height);
        // Building(Polygon p, float h);
        virtual void Render(RENDER_MODE mode);
};