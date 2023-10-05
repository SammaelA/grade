#pragma once
#include "cities_generator/global.h"
#include "RenderableAsMap.h"
#include "Building.h"

using glm::vec2;

struct RoadSection;

struct RoadNode
{
    vec2 pos;
    float data;
    bool blocked = false;
    std::vector<RoadSection*> edges;

    RoadNode(vec2 p);

    bool IsCrossroad();
    bool IsEndNode();
    bool IsTurn();
    bool operator== (const RoadNode& a);
    bool operator!= (const RoadNode& a);
};


struct RoadSection
{
    RoadNode *start, *end;
    glm::vec2 startDir, endDir;
    bool blocked = false;

    bool IsStraight();
    bool IsRoundTurn();
    glm::vec2 StraightDir(RoadNode* from = nullptr);
    RoadSection(RoadNode *, RoadNode *, glm::vec2 _startDir = vec2{0}, glm::vec2 _endDir = vec2{0});
    RoadSection(RoadNode *, RoadNode *, float circularAngleDegrees);
    std::array<RoadNode*, 2> GetEdges();
    bool operator== (const RoadSection& a);
    float GetLength();
    static bool AreParrallel(RoadSection*, RoadSection*);
    RoadNode* Other(RoadNode* n);
    static RoadSection* SectionByEdges(RoadNode* , RoadNode*);
    glm::vec2 RoundTurnCenter();
};

struct RoadGraph
{
    RoadNode* root;
    RoadNode* AddNode(glm::vec2 pos);
    RoadSection* AddSection(RoadNode*, RoadNode*, glm::vec2, glm::vec2);
    RoadSection* AddSection(RoadNode*, RoadNode*, float);
    void RemoveSection(RoadNode*, RoadNode*);
    void RemoveNode(RoadNode*);
    void MergeNodes(RoadNode* from, RoadNode* to);
    std::vector<RoadSection*> GetSections();
    std::vector<RoadNode*> GetNodes();
    RoadGraph();
    RoadNode* DivideSection(RoadNode* a, RoadNode* b, float proportion);
};

class RoadScheme : public RenderableAsMap
{
    friend class Landscape;

    static constexpr int MAX_SIZE = 512 * 512;
    float widthToHeightRatio;

    protected:
        RoadGraph graph;

        RoadScheme(float ratio);
        virtual bool IsTextureReady();
        virtual void SaveAsTexture();
        
    public:
        virtual GLuint GetTexture();
        GLuint& GetTextureRef();
};