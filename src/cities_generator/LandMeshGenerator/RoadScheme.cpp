#include "RoadScheme.h"
#include "third_party/stb_image.h"
#include <algorithm>
#include <functional>

bool RoadScheme::IsTextureReady()
{
    return m_textureId != BAD_TEXTURE;
}

void RoadScheme::SaveAsTexture()
{
    //Preparing
    if (widthToHeightRatio <= 0)
    {
        debug("CANT SAVE TEXTURE FOR ROAD");
        throw std::exception{};
    }

    int TEX_WIDTH = sqrtf(RoadScheme::MAX_SIZE * widthToHeightRatio);
    int TEX_HEIGHT = sqrtf(RoadScheme::MAX_SIZE / widthToHeightRatio);

    GLubyte* image = (GLubyte*)malloc(TEX_WIDTH * TEX_HEIGHT * 4 * sizeof(GLubyte));

    std::function<void(vec2Int, unsigned)> SetColor =
        [&TEX_WIDTH, &TEX_HEIGHT, &image](vec2Int pos, unsigned col)
    {
        image[0 + (pos.x + pos.y*(TEX_WIDTH))*4] = (GLubyte) (col >> 24) ;
        image[1 + (pos.x + pos.y*(TEX_WIDTH))*4] = (GLubyte) (col << 8 >> 24);
        image[2 + (pos.x + pos.y*(TEX_WIDTH))*4] = (GLubyte) (col << 16 >> 24);
        image[3 + (pos.x + pos.y*(TEX_WIDTH))*4] = (GLubyte) (col << 24 >> 24);
    };
    
    for (unsigned j = 0; j < TEX_HEIGHT; j++) {
        for (unsigned i = 0; i < TEX_WIDTH; i++) {

            unsigned color = (((204u << 16) | (202u << 8) | (102)) << 8) | 255;
            SetColor(vec2Int{(int)i, (int)j}, color);
        }
    }

    //Drawning road lines

    std::function<vec2Int(vec2)> PointNormalizedToInt = [&TEX_WIDTH, &TEX_HEIGHT](vec2 pos)
    {
        if (pos.x < 0 || pos.y < 0 || pos.x > 1 || pos.y > 1)
        {
            debug("BAD INPUT! [PointNormalizedToInt RoadScheme]", pos.x, pos.y);
            throw std::exception{};
        }
        pos.x = clamp(pos.x, 0.0f, 0.999999f);
        pos.y = clamp(pos.y, 0.0f, 0.999999f);
        vec2Int res;
        res.x = int(floorf(pos.x * TEX_WIDTH));
        res.y = int(floorf(pos.y * TEX_HEIGHT));
        return res;
    };

    float stepLength = 1.0f / (std::max(TEX_WIDTH, TEX_HEIGHT) + 1); 
    for (RoadSection* section : graph.GetSections())
    {
        glm::vec2 start = section->start->pos;
        glm::vec2 end = section->end->pos;
        glm::vec2 curPoint = start;
        while (glm::length(curPoint - start) - stepLength < glm::length(end - start))
        {
            glm::vec2 realPoint = curPoint;
            if (glm::length(curPoint - start) > glm::length(end - start))
                realPoint = end;
            vec2Int curCoords = PointNormalizedToInt(realPoint);
            SetColor(curCoords, (((0 << 16) | (0 << 8) | (0)) | 255));
            curPoint += stepLength * glm::normalize(end - start);
        }
    }

    //Saving
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &m_textureId);
    glBindTexture(GL_TEXTURE_2D, m_textureId);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, 
                    GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
                    GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEX_WIDTH, 
                TEX_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, 
                image);
    stbi_image_free(image);
    glBindTexture(GL_TEXTURE_2D, 0);
}

RoadScheme::RoadScheme(float ratio)
{
    widthToHeightRatio = ratio;
    m_textureId = BAD_TEXTURE;
}

GLuint RoadScheme::GetTexture()
{
    if (IsTextureReady())
        return m_textureId;
    else
    {
        debug("ROAD TEXTURE IS NOT READY");
        throw std::exception{};
    }
}

GLuint& RoadScheme::GetTextureRef()
{
    if (IsTextureReady())
        return m_textureId;
    else
    {
        debug("ROAD TEXTURE IS NOT READY");
        throw std::exception{};
    }
}

RoadNode::RoadNode(vec2 _pos)
{
    pos = _pos;
}

RoadSection::RoadSection(RoadNode *n1, RoadNode *n2, glm::vec2 startD, glm::vec2 endD)
{
    start = n1;
    end = n2;
    startDir = glm::normalize((glm::length(startD) < 0.00001) ? StraightDir() : startD);
    endDir = glm::normalize((glm::length(endD) < 0.00001) ? StraightDir() : endD);
}

RoadSection::RoadSection(RoadNode *n1, RoadNode *n2, float circularAngleDegrees)
{
    glm::vec2 forward = glm::normalize(n2->pos - n1->pos);
    glm::vec2 left = LeftNormal(forward);
    float angle = glm::radians(circularAngleDegrees);
    glm::vec2 startDir = cosf(angle) * forward + sinf(angle) * left;
    glm::vec2 endDir = cosf(angle) * forward - sinf(angle) * left;
    *this = RoadSection(n1, n2, startDir, endDir);
}

bool RoadSection::IsStraight()
{
    return glm::length(startDir - endDir) < 0.00001 && glm::length(startDir - StraightDir()) < 0.0001;
}

bool RoadSection::IsRoundTurn()
{
    if (!Approx(glm::dot(StraightDir(), startDir), glm::dot(StraightDir(), endDir)))
        return false;
    if (IsStraight())
        return false;
    return !(fabs(glm::dot(startDir, endDir) - 1.f) < 0.00001);
}

glm::vec2 RoadSection::StraightDir(RoadNode* from)
{
    if (from == nullptr || from == start)
        return glm::normalize(end->pos - start->pos);
    else if (from == end)
        return glm::normalize(start->pos - end->pos);
    else
    {
        debug("INCORRECT FROM NODE!");
        throw std::exception{};
    }
    
}

RoadNode* RoadGraph::AddNode(glm::vec2 pos)
{
    RoadNode* res = new RoadNode(pos);
    if (root == nullptr)
        root = res;
    return res;
}

RoadSection* RoadGraph::AddSection(RoadNode* start, RoadNode* end, glm::vec2 startDir, glm::vec2 endDir)
{
    RoadSection* section = new RoadSection(start, end, startDir, endDir);
    start->edges.push_back(section);
    end->edges.push_back(section);
    return section;
}

RoadSection* RoadGraph::AddSection(RoadNode* start, RoadNode* end, float circularAngleDegrees)
{
    RoadSection* section = new RoadSection(start, end, circularAngleDegrees);
    start->edges.push_back(section);
    end->edges.push_back(section);
    return section;
}

void RoadGraph::RemoveSection(RoadNode* n1, RoadNode* n2)
{
    if (n1 == root || n2 == root)
    {
        if (n1 == root && n1->edges.size() == 1)
        {
            if (n2->edges.size() > 1)
                root = n2;
            else
            {
                debug("UNABLE TO REDEFINE ROOT");
                throw std::exception{};
            }
        }
        if (n2 == root && n2->edges.size() == 1)
        {
            if (n1->edges.size() > 1)
                root = n1;
            else
            {
                debug("UNABLE TO REDEFINE ROOT");
                throw std::exception{};
            }
        }
    }
    
    for (int i = 0; i < n1->edges.size(); i++)
    {
        RoadSection* section = n1->edges[i];
        if (section->end == n2 || section->start == n2)
        {
            n2->edges.erase(std::find(n2->edges.begin(), n2->edges.end(), section));
            n1->edges.erase(n1->edges.begin() + i);
            delete section;
            return;
        }
    }
    debug("THERE IS NO SUCH SECTION!");
    throw std::exception{};
}

RoadGraph::RoadGraph()
{
    root = nullptr;
}

std::vector<RoadSection*> RoadGraph::GetSections()
{
    std::unordered_set<RoadNode*> visited;
    std::vector<RoadSection*> result;
    std::unordered_set<RoadNode*> reached;
    
    if (root == nullptr)
    {
        debug("ROOT IS NULL");
        throw std::exception{};
    }
    reached.emplace(root);

    while(reached.size() != 0)
    {
        auto currentIt = reached.begin();
        RoadNode* current = *currentIt;
        reached.extract(currentIt);
        visited.emplace(current);
        for (RoadSection* section : current->edges)
        {
            RoadNode* otherNode = (section->start == current) ? section->end : section->start;
            if (visited.find(otherNode) == visited.end())
            {
                result.push_back(section);
                reached.emplace(otherNode);
            }
        }
    }
    return result;
}

std::vector<RoadNode*> RoadGraph::GetNodes()
{
    std::unordered_set<RoadNode*> visited;
    std::unordered_set<RoadNode*> reached;
    
    if (root == nullptr)
    {
        debug("ROOT IS NULL");
        throw std::exception{};
    }
    reached.emplace(root);

    while(reached.size() != 0)
    {
        auto currentIt = reached.begin();
        RoadNode* current = *currentIt;
        reached.extract(currentIt);
        visited.emplace(current);
        for (RoadSection* section : current->edges)
        {
            RoadNode* otherNode = (section->start == current) ? section->end : section->start;
            if (visited.find(otherNode) == visited.end())
            {
                reached.emplace(otherNode);
            }
        }
    }
    
    std::vector<RoadNode*> result;
    result.insert(result.end(), visited.begin(), visited.end());
    return result;
}

bool RoadNode::IsCrossroad()
{
    return edges.size() > 2;
}

bool RoadSection::operator== (const RoadSection& a)
{
    return (*(a.start) == *start) && (*(a.end) == *end);
}

bool RoadNode::operator== (const RoadNode& a)
{
    return glm::length(a.pos - pos) < 0.00001;
}

bool RoadNode::operator!= (const RoadNode& a)
{
    return !(*this == a);
}


std::array<RoadNode*, 2> RoadSection::GetEdges()
{
    return std::array<RoadNode*, 2>{start, end};
}

RoadNode* RoadGraph::DivideSection(RoadNode* a, RoadNode* b, float proportion)
{
    RoadNode* newNode = AddNode(LERP(a->pos, b->pos, proportion));
    AddSection(a, newNode, 0.f);
    AddSection(newNode, b, 0.f);
    RemoveSection(a, b);
    return newNode;
}

bool RoadNode::IsEndNode()
{
    return edges.size() == 1;
}

void RoadGraph::MergeNodes(RoadNode* from, RoadNode* to)
{
    for (RoadSection* sec: from->edges)
    {
        RoadNode* s = (sec->start == from) ? sec->end : sec->start;
        RemoveSection(sec->start, sec->end);
        if (s != to)
            auto a = AddSection(s, to, 0.f);
    }
}

void RoadGraph::RemoveNode(RoadNode* n)
{
    bool ok = !(n == root);
    for (RoadSection* s: n->edges)
    {
        RoadNode* other = (s->start == n) ? s->end : s->start;
        RemoveSection(n, other);
        if (other->edges.size() == 0)
        {
            delete other;
        }
        else
        {
            if (!ok)
            {
                root = other;
                ok = true;
            }
        }
    }
    delete n; 
    if (!ok)
    {
        debug("UNABLE TO SAVE ROOT");
        throw std::exception{};
    }
}

float RoadSection::GetLength()
{
    return glm::length(start->pos - end->pos);
}

bool RoadNode::IsTurn()
{
    return edges.size() == 2;
}

bool RoadSection::AreParrallel(RoadSection* s1, RoadSection* s2)
{
    if (s1->IsStraight() && s2->IsStraight())
    {
        return Approx(fabsf(glm::dot(s1->StraightDir(), s2->StraightDir())), 1.f);
    }
    return false;
}

RoadNode* RoadSection::Other(RoadNode* n)
{
    if (n == start)
        return end;
    else if (n == end)
        return start;
    else
    {
        debug("IT IS NOT SECTION'S NODE!");
        throw std::exception{};
    }
}

RoadSection* RoadSection::SectionByEdges(RoadNode* n1, RoadNode* n2)
{
    for (int i = 0; i < n1->edges.size(); i++)
    {
        RoadSection* section = n1->edges[i];
        if (section->end == n2 || section->start == n2)
        {
            return section;
        }
    }
    debug("THERE IS NO SECTION BETWEEN!");
    throw std::exception{};
}

glm::vec2 RoadSection::RoundTurnCenter()
{
    if (!IsRoundTurn())
    {
        debug("IT IS NOT ROUND TURN!");
        throw std::exception{};
    }
    if (Approx(glm::normalize(startDir), -glm::normalize(endDir), 0.0001))
    {
        return (start->pos + end->pos) * 0.5f;
    }
    float angleRadians = SafeAcosf(glm::dot(glm::normalize(startDir), glm::normalize(-endDir)));
    angleRadians /= 2;
    float dist = tan(angleRadians) / sinf(angleRadians) * (glm::length(start->pos - end->pos) / 2);
    
    glm::vec2 _startDir = (glm::dot(startDir, StraightDir()) > 0) ? startDir : -startDir;
    glm::vec2 _endDir = (glm::dot(endDir, StraightDir()) > 0) ? endDir : -endDir;

    glm::vec2 dir = LeftNormal(_startDir);
    dir = (glm::dot(dir, _endDir) > 0) ? dir : -dir;
    
    return start->pos + dir * dist;
}