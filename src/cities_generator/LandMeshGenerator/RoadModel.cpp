#include "RoadModel.h"
#include "../ECSClass.h"
#include "third_party/stb_image.h"
#include "cities_generator/ImageLoader.h"
#include "LandGeneratorClass.h"

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;

RoadModel::RoadModel()
{
    m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer("roadShader");
    // asphaltTex = ImageLoader::StbLoadImage(
    //     std::string(ASSETS_FOLDER) + "Textures/chessTexture.jpg");
    asphaltTex = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/asphaltTexture.jpg");
    markingTex = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/markingTexture.png");
    crossroadMarkingTex = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/markingTexture.png");

    glBindTexture(GL_TEXTURE_2D, markingTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, crossroadMarkingTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    UseShader();

    Scene3DObjectDataHolder::VerifyUnoformPos(viewProjMatrix_Uniform, "view_proj_matrix");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunPos_Uniform, "sun_pos");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunColor_Uniform, "sun_color");
    Scene3DObjectDataHolder::VerifyUnoformPos(asphaltTex_Uniform, "asphalt_texture");
    Scene3DObjectDataHolder::VerifyUnoformPos(asphaltTexScale_Uniform, "asphalt_scale");
    Scene3DObjectDataHolder::VerifyUnoformPos(markingTex_Uniform, "marking_texture");
    Scene3DObjectDataHolder::VerifyUnoformPos(crossroadMarkingTex_Uniform, "cross_marking_texture");
}

void RoadModel::Render(RENDER_MODE mode)
{
    //SET STATE
    glProvokingVertex(GL_LAST_VERTEX_CONVENTION);
    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);

    //UNIFORMS
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    glm::mat4 viewProj = cam->GetViewProjectionMatrix();

    Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
            2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];

    TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
    SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();

    UseShader();
    
    LoadInShader_Scene3DObject();
    glUniformMatrix4fv(viewProjMatrix_Uniform, 1, GL_FALSE, glm::value_ptr(viewProj));
    glUniform3fv(sunPos_Uniform, 1, glm::value_ptr(sunTransform->position));
    glUniform3fv(sunColor_Uniform, 1, glm::value_ptr(sunComp->color));
    glUniform1f(asphaltTexScale_Uniform, asphaltTexScale);
    
    //TEXTURES
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, asphaltTex);
    glUniform1i(asphaltTex_Uniform, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, markingTex);
    glUniform1i(markingTex_Uniform, 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, crossroadMarkingTex);
    glUniform1i(crossroadMarkingTex_Uniform, 2);

    //DRAW CALL
    RenderableObjectDataHolder::Render(mode);
}

float RoadSlice::GetModelScale()
{
    if (modelScale < 0)
    {
        float sliceWidth = 0;
        float _min = borderPattern[0].x;
        float _max = borderPattern[0].x;
        for (int i = 0; i < borderPattern.size(); i++)
        {
            _min = std::min(borderPattern[i].x, _min);
            _max = std::max(borderPattern[i].x, _max);
        }
        sliceWidth = _max - _min;
        modelScale = GetDefaultFullSlicewidth() / sliceWidth;
    }
    return modelScale;
}

void RoadSlice::CalculateDistances()
{
    distances = std::vector<float>{};
    float currentLen = 0;
    for (int i = 0; i < border.size(); i++)
    {
        if (i != 0)
            currentLen += glm::length(border[i] - border[i - 1]);
        
        distances.push_back(currentLen);
    }
}

void RoadSlice::CalculateNormals(RoadSlice* prevPtr, RoadSlice* nextPtr)
{
    normals.clear();
    for (int i = 0; i < border.size(); i++)
    {
        glm::vec3 center = border[i];
        std::array<glm::vec3, 4> points;
        points[0] = (nextPtr != nullptr) ?  nextPtr->border[i] : center;
        points[1] = (i != 0) ?  border[i - 1] : center;
        points[2] = (prevPtr != nullptr) ?  prevPtr->border[i] : center;
        points[3] = (i != border.size() - 1) ?  border[i + 1] : center;

        glm::vec3 result{0};
        for (int j = 0; j < 4; j++)
        {
            glm::vec3 e1 = points[j] - center;
            glm::vec3 e2 = points[(j + 1) % 4] - center;
            float weight = sinf(SafeAcosf(glm::dot(glm::normalize(e1), glm::normalize(e2))));
            weight *= glm::length(e1) * glm::length(e2) * 0.5;
            
            if (weight > 0.000001)
            {
                glm::vec3 n = glm::normalize(glm::cross(e1, e2));
                if (glm::dot(n, glm::vec3(0, 1, 0)) < 0.f)
                    n = -n;
                result += n * weight;
            }
        }
        result = glm::normalize(result);
        normals.push_back(result);
    }
}

void RoadSlice::CalculateNormals(glm::vec3 groundNormal, glm::vec3 realFoundation)
{
    normals.clear();
    RoadSlice temp(glm::vec3{0});
    glm::vec3 verticalDirection = glm::normalize(groundNormal);
    glm::vec3 horizontalDirection = Approx(realFoundation, glm::vec3{0}) ?
        glm::cross(verticalDirection, glm::cross(verticalDirection, foundation)) :
        glm::cross(verticalDirection, glm::cross(verticalDirection, realFoundation));
    horizontalDirection = glm::normalize(horizontalDirection);
    if (glm::dot(horizontalDirection, foundation) < 0)
        horizontalDirection *= -1;

    auto planarInds = RoadSlice::GetPlanarBorderIndexes();

    for (int i = 0; i < temp.border.size(); i++)
    {
        glm::vec2 resultingNormal2{0};
        for (int j = -1; j <= 1; j+=2)
        {
            if (i + j < 0 || i + j >= temp.border.size())
                continue;
            float weight = fabsf(temp.distances[i] - temp.distances[i + j]);
            glm::vec2 curNormal = XY(temp.border[i]) - XY(temp.border[i + j]);
            curNormal = LeftNormal(curNormal);
            if (glm::dot(curNormal, glm::vec2(0, 1)) < 0)
                curNormal = -curNormal;
            resultingNormal2 += curNormal * weight;
        }
        resultingNormal2 = glm::normalize(resultingNormal2);
        glm::vec3 resultingNormal3 = resultingNormal2.x * horizontalDirection + resultingNormal2.y * verticalDirection;
        resultingNormal3 = glm::normalize(resultingNormal3);
        normals.push_back(resultingNormal3);
    }
}

float RoadSlice::GetDefaultFullSlicewidth()
{
    return Landscape::RoadModelGenerator::ROAD_FULL_REAL_WIDTH;
}

std::array<glm::vec3, 2> RoadSlice::GetBoundRealPoints()
{
    return std::array<glm::vec3, 2> {border[0], border[border.size() - 1]};
}

RoadSlice::RoadSlice(glm::vec3 _position, glm::vec3 _foundation, float s)
{
    sliceHorizontalScale = s;
    position = _position;
    foundation = glm::normalize(_foundation);
    
    float rollAngle = asinf(foundation.y);
    float yawAngle = -atan2f(foundation.z, foundation.x);

    for (int i = 0; i < RoadSlice::borderPattern.size(); i++)
    {
        glm::vec3 curPos = RoadSlice::borderPattern[i];
        curPos *= RoadSlice::GetModelScale();
        curPos.x *= sliceHorizontalScale;
        curPos = glm::rotateZ(curPos, rollAngle);
        curPos = glm::rotateY(curPos, yawAngle);
        curPos += position;
        border.push_back(curPos);
    }
    CalculateDistances();
}

void RoadSlice::AddFacesToBuffer(GLfloat* buffer)
{
    debug("NOT SUPPORTED BEZIER CURVES! AND NORMALS!");
    throw std::exception{};

    GLfloat* ptr = buffer;
    for (int i = 0; i < border.size() - 2; i++)
    {
        for (int j = 0; j < TRIANGLE_SIZE_IN_VERTEX; j++)
        {
            int currentBorderInd = 0;
            if (j != 0)
                currentBorderInd = i + j;
            memcpy(
                ptr, 
                glm::value_ptr(border[currentBorderInd]), 
                sizeof(border[currentBorderInd])
            );
            ptr += sizeof(border[currentBorderInd]) / sizeof(border[currentBorderInd].x);
            glm::vec2 asphaltTexCoords{0};
            memcpy(ptr, glm::value_ptr(asphaltTexCoords), sizeof(asphaltTexCoords));
            ptr += (sizeof(asphaltTexCoords) / sizeof(asphaltTexCoords.x));
            
            float markingTexCoords = 0.f;
            memcpy(ptr, (void*)(&markingTexCoords), sizeof(markingTexCoords));
            ptr += (sizeof(markingTexCoords) / sizeof(markingTexCoords));
        }
    }
}
void RoadSlice::AddRoadConnectionToBuffer(
    const RoadSlice& originSlice, 
    const RoadSlice& targetSlice,    
    GLfloat* buffer, 
    glm::vec2 roadDir)
{
    auto planarInds = RoadSlice::GetPlanarBorderIndexes();
    bool isRoundTurn = (glm::length(roadDir) < 0.0001);
    glm::vec4 RoundCenterSmallRadiuses{0};
    glm::vec4 RoundBigRadiusesEmpty{0};
    if (isRoundTurn)
    {
        std::array<glm::vec2, 2> s1Positions{XZ(targetSlice.border[planarInds[0]]), XZ(targetSlice.border[planarInds[1]])};
        std::array<glm::vec2, 2> s2Positions{XZ(originSlice.border[planarInds[0]]), XZ(originSlice.border[planarInds[1]])};
        glm::vec2 s1Center = (s1Positions[0] + s1Positions[1]) * 0.5f;
        glm::vec2 s2Center = (s2Positions[0] + s2Positions[1]) * 0.5f;
        glm::vec2 middle = (s1Center + s2Center) * 0.5f;
        float angle = glm::dot(
            glm::normalize(XZ((targetSlice.border[planarInds[0]] - targetSlice.border[planarInds[1]]) * 0.5f)),
            glm::normalize(XZ((originSlice.border[planarInds[0]] - originSlice.border[planarInds[1]]) * 0.5f))
        );
        angle = fabsf(angle);
        angle = SafeAcosf(angle);
        angle *= 0.5;
        float normalLength = glm::length(s2Center - s1Center) * 0.5 / glm::tan(angle);
        bool sideFlag = glm::length(s1Positions[0] - s2Positions[0]) < glm::length(s1Positions[1] - s2Positions[1]);
        glm::vec2 direction = (sideFlag) ?
            (s1Positions[0] + s2Positions[0]) * 0.5f - (s1Positions[1] + s2Positions[1]) * 0.5f :
            (s1Positions[1] + s2Positions[1]) * 0.5f - (s1Positions[0] + s2Positions[0]) * 0.5f;
        direction = glm::normalize(direction);
        glm::vec2 roundCenter = middle + direction * normalLength;
        glm::vec2 radiusesOne = glm::vec2 {glm::length(roundCenter - s1Positions[0]), glm::length(roundCenter - s2Positions[0])};
        glm::vec2 radiusesTwo = glm::vec2 {glm::length(roundCenter - s1Positions[1]), glm::length(roundCenter - s2Positions[1])};
        RoundCenterSmallRadiuses = (sideFlag) ? 
            glm::vec4{roundCenter, radiusesOne.x, radiusesOne.y} :
            glm::vec4{roundCenter, radiusesTwo.x, radiusesTwo.y};
        RoundBigRadiusesEmpty = (sideFlag) ? 
            glm::vec4{radiusesTwo.x, radiusesTwo.y, 0, 0} :
            glm::vec4{radiusesOne.x, radiusesOne.y, 0, 0};
    }
        
    GLfloat* ptr = buffer;
    for (int i = 0; i < targetSlice.border.size() - 1; i++)
    {
        for (int k = 0; k < 2; k++)
        {
            const RoadSlice& from = (k == 0) ? targetSlice : originSlice;
            const RoadSlice& to = (k == 0) ? originSlice : targetSlice;
            for (int j = 0; j < TRIANGLE_SIZE_IN_VERTEX; j++)
            {
                vec2Int currentBorderInd; // x - ind for [from, to], y - ind for border
                if (j != 2)
                {
                    currentBorderInd.x = 0;
                    currentBorderInd.y = i + j;
                }
                else
                {
                    currentBorderInd.x = 1;
                    currentBorderInd.y = i + ((k == 0) ? 0 : 1);
                }

                const RoadSlice& cur = (currentBorderInd.x == 0) ? from : to;
                memcpy(ptr, glm::value_ptr(cur.border[currentBorderInd.y]), sizeof(cur.border[0]));
                ptr += (sizeof(cur.border[0]) / sizeof(cur.border[0].x));

                glm::vec3 normal = cur.normals[currentBorderInd.y];
                memcpy(ptr, glm::value_ptr(normal), sizeof(normal));
                ptr += (sizeof(normal) / sizeof(normal.x));

                glm::vec2 asphaltTexCoords;
                if (!isRoundTurn)
                {
                    glm::vec2 left = LeftNormal(roadDir);
                    if (currentBorderInd.y < planarInds[0])
                    {
                        asphaltTexCoords = XZ(cur.border[planarInds[0]]) + left * 
                            (cur.distances[planarInds[0]] - cur.distances[currentBorderInd.y]);
                    }
                    else if (currentBorderInd.y > planarInds[1])
                    {
                        asphaltTexCoords = XZ(cur.border[planarInds[1]]) - left * 
                            (cur.distances[currentBorderInd.y] - cur.distances[planarInds[1]]);
                    }
                    else
                    {
                        asphaltTexCoords = XZ(cur.border[currentBorderInd.y]);
                    }
                }
                else
                {
                    glm::vec2 curSliceCenter = XZ((cur.border[planarInds[0]] + cur.border[planarInds[1]]) * 0.5f);
                    glm::vec2 dirStartToCenter = glm::normalize(XZ(cur.foundation));
                    float fullDist = cur.distances[cur.distances.size() - 1];
                    asphaltTexCoords = curSliceCenter + dirStartToCenter * (cur.distances[currentBorderInd.y] - fullDist / 2);
                }
                memcpy(ptr, glm::value_ptr(asphaltTexCoords), sizeof(asphaltTexCoords));
                ptr += (sizeof(asphaltTexCoords) / sizeof(asphaltTexCoords.x));
                
                float markingTexCoords = 0.f;
                if (!isRoundTurn)
                {
                    float planarDist = fabsf(cur.distances[planarInds[0]] - cur.distances[planarInds[1]]);
                    markingTexCoords = clamp(
                        (cur.distances[currentBorderInd.y] - cur.distances[planarInds[0]]) / planarDist, 0.f, 1.f);
                }
                else
                {
                    //2 means targetSlice points, 3 means originSlice points
                    markingTexCoords = Approx(cur.border[planarInds[0]], originSlice.border[planarInds[0]])
                        ? 2.f
                        : 3.f;
                }
                memcpy(ptr, (void*)(&markingTexCoords), sizeof(markingTexCoords));
                ptr += (sizeof(markingTexCoords) / sizeof(markingTexCoords));

                memcpy(ptr, glm::value_ptr(RoundCenterSmallRadiuses), sizeof(RoundCenterSmallRadiuses));
                ptr += (sizeof(RoundCenterSmallRadiuses) / sizeof(RoundCenterSmallRadiuses.x));

                memcpy(ptr, glm::value_ptr(RoundBigRadiusesEmpty), sizeof(RoundBigRadiusesEmpty));
                ptr += (sizeof(RoundBigRadiusesEmpty) / sizeof(RoundBigRadiusesEmpty.x));
            }
        } 
    }
}
uint RoadSlice::GetSliceVertexesCount()
{
    return TRIANGLE_SIZE_IN_VERTEX * (borderPattern.size() - 2);
}
uint RoadSlice::GetConnectionVertexesCount()
{
    return TRIANGLE_SIZE_IN_VERTEX * (borderPattern.size() - 1) * 2;
}
std::array<int, 2> RoadSlice::GetPlanarBorderIndexes()
{
    return std::array<int, 2>{2, 3};
}
float RoadSlice::GetDefaultPlanarSlicewidth()
{
    glm::vec3 p1 = borderPattern[GetPlanarBorderIndexes()[0]];
    glm::vec3 p2 = borderPattern[GetPlanarBorderIndexes()[1]];
    return glm::length(p1 - p2) * GetModelScale();
}

void UpdateBufferCapacity(int vertexesAmount, int extraVertexes, int& availableAmount, GLfloat*& generalBuffer)
{
    const int FLOAT_SIZE_IN_BYTE = sizeof(generalBuffer[0]);
    while(vertexesAmount + extraVertexes > availableAmount)
    {
        int save = availableAmount;
        availableAmount *= 2;
        GLfloat* temp = new GLfloat[availableAmount * RoadSlice::VERTEX_SIZE_IN_FLOAT];
        memcpy(temp, generalBuffer, save * RoadSlice::VERTEX_SIZE_IN_FLOAT * FLOAT_SIZE_IN_BYTE);
        delete [] generalBuffer;
        generalBuffer = temp;
    }
}

//IF YOU CHANGE IT, CHANGE (GetPlanarBorderIndexes, CalculateNormals(glm::vec3))
const std::vector<glm::vec3> RoadSlice::borderPattern{
    {-9.7, -2.234332, 0},
    {-8.23936, -0.56, 0},
    {-7.30317, 0, 0},
    {7.30317, 0, 0},
    {8.23936, -0.56, 0},
    {9.7, -2.234332, 0}
};
const glm::vec3 RoadSlice::foundationDefault{1,0,0};
float RoadSlice::modelScale = -1;