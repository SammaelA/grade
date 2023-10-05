#pragma once
#include "cities_generator/global.h"
#include "../ShaderDataHolders/Archetypes/BasicGLTFModel.h"
#include "../ShaderDataHolders/RenderableObjectDataHolder.h"

class RoadModel : public BasicGLTFModel
{
    int asphaltTex_Uniform = MISSING_UNIFORM_POS;
    int asphaltTexScale_Uniform = MISSING_UNIFORM_POS;
    int markingTex_Uniform = MISSING_UNIFORM_POS;
    int crossroadMarkingTex_Uniform = MISSING_UNIFORM_POS;

    GLuint asphaltTex;
    GLfloat asphaltTexScale = 0.2;
    GLuint markingTex, crossroadMarkingTex;

    public:
        RoadModel();
        virtual void Render(RENDER_MODE mode);
};

struct RoadSlice
{
    static const std::vector<glm::vec3> borderPattern;
    static const int POSITION_SIZE_IN_FLOAT = 3;
    static const int NORMAL_SIZE_IN_FLOAT = 3;
    static const int ASPHALT_TEX_COORDS_SIZE_IN_FLOAT = 2;
    static const int MARKING_TEX_COORDS_SIZE_IN_FLOAT = 1;
    static const int BEZIER_HANDLES_SIZE_IN_FLOAT = 4;
    static const int BEZIER_ANCHORS_SIZE_IN_FLOAT = 4;
    static const int VERTEX_SIZE_IN_FLOAT = 
        ASPHALT_TEX_COORDS_SIZE_IN_FLOAT + POSITION_SIZE_IN_FLOAT + MARKING_TEX_COORDS_SIZE_IN_FLOAT +
        BEZIER_ANCHORS_SIZE_IN_FLOAT + BEZIER_HANDLES_SIZE_IN_FLOAT + NORMAL_SIZE_IN_FLOAT;
    static const int TRIANGLE_SIZE_IN_VERTEX = 3;
    static const glm::vec3 foundationDefault;
    static float modelScale;

    std::vector<glm::vec3> border;
    std::vector<glm::vec3> normals;
    std::vector<float> distances;
    glm::vec3 position, foundation;
    float sliceHorizontalScale;

    RoadSlice() = default;
    RoadSlice(glm::vec3 _position, glm::vec3 _foundation = RoadSlice::foundationDefault, float s = 1.f);

    static float GetModelScale();
    static float GetDefaultFullSlicewidth();
    static void AddRoadConnectionToBuffer(const RoadSlice&, const RoadSlice&, GLfloat*, glm::vec2 = glm::vec2{0});
    static uint GetSliceVertexesCount();
    static uint GetConnectionVertexesCount();
    static std::array<int, 2> GetPlanarBorderIndexes();
    static float GetDefaultPlanarSlicewidth();
    void CalculateDistances();
    void CalculateNormals(RoadSlice* prevPtr, RoadSlice* nextPtr);
    void CalculateNormals(glm::vec3 groundNormal, glm::vec3 realFoundation = glm::vec3{0});
    std::array<glm::vec3, 2> GetBoundRealPoints();
    void AddFacesToBuffer(GLfloat* buffer);
};

void UpdateBufferCapacity(int vertexesAmount, int extraVertexes, int& availableAmount, GLfloat*& generalBuffer);