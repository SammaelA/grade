#include "Building.h"
#include "../ECSClass.h"
#include "LandGeneratorClass.h"
#include "cities_generator/ImageLoader.h"

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;

std::map<BUILDING_TYPE, BuildingInfo> Building::buildingTypes{};

Building::Building(Basement _basement, BUILDING_TYPE type, float _height) : basement(_basement)
{
    InitBuildingTypes();

    BuildingInfo curInfo = buildingTypes.at(type);
    // BuildingInfo curInfo = buildingTypes.at(BUILDING_TYPE::HIGH_RISE);

    basement.size.x = floorf(basement.size.x / curInfo.GetRealSize().x) * curInfo.GetRealSize().x;
    basement.size.y = floorf(basement.size.y / curInfo.GetRealSize().x) * curInfo.GetRealSize().x;

    height = floorf(_height / curInfo.GetRealSize().y) * curInfo.GetRealSize().y;
    elevation = basement.height;
    position = X0Z(basement.center);
    position.y = elevation;
    rotation = glm::yawPitchRoll(-basement.angleRad, 0.f ,0.f);

    //preload
    m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer("buildingShader");

    UseShader();

    Scene3DObjectDataHolder::VerifyUnoformPos(viewProjMatrix_Uniform, "view_proj_matrix");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunPos_Uniform, "sun_pos");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunColor_Uniform, "sun_color");
    Scene3DObjectDataHolder::VerifyUnoformPos(cameraPos_Uniform, "camera_position");
    Scene3DObjectDataHolder::VerifyUnoformPos(buildingTex_Uniform, "building_texture");
    Scene3DObjectDataHolder::VerifyUnoformPos(roofTex_Uniform, "roof_texture");
    Scene3DObjectDataHolder::VerifyUnoformPos(buildingSpecTex_Uniform, "specular_texture");


    //Vertex data
    const int POLYGON_SIZE_IN_VERTEX = 3;
    const int POSITION_SIZE_IN_FLOAT = 3;
    const int NORMAL_SIZE_IN_FLOAT = 3;
    const int TC_SIZE_IN_FLOAT = 3;
    const int VERTEX_SIZE_IN_FLOAT = POSITION_SIZE_IN_FLOAT + NORMAL_SIZE_IN_FLOAT + TC_SIZE_IN_FLOAT;

    int SIDE_SIZE_IN_POLYGON = (basement.GetPoints().size() - 2);
    int SIDE_COUNT = 6;
    int POLYGON_COUNT = SIDE_SIZE_IN_POLYGON * SIDE_COUNT;

    int _vertexNum = POLYGON_COUNT * POLYGON_SIZE_IN_VERTEX;
    GLfloat* generalBuffer = new GLfloat[_vertexNum * VERTEX_SIZE_IN_FLOAT];
    const int FLOAT_SIZE_IN_BYTE = sizeof(generalBuffer[0]);


    std::array<glm::vec2, 4> centeredPoints;
    for (int i = 0; i < 4; i++)
    {
        bool isLeft = (i < 2);
        bool isTop = (i % 3 != 0);
        centeredPoints[i] = glm::vec2{basement.size.x / 2 * (isLeft ? -1 : 1), basement.size.y / 2 * (isTop ? 1 : -1)};
    }
    
    enum {
        LEFT = 0,
        RIGHT = 1,
        FORWARD = 2,
        BACK = 3,
        TOP = 4,
        DOWN = 5
    } side;
    
    GLfloat* ptr = generalBuffer;
    for (int side = 0; side < SIDE_COUNT; side++)
    {
        std::array<glm::vec3, 4> curPoints;
        for (int i = 0; i < 4; i++)
        {
            switch (side)
            {
                case LEFT:
                    curPoints[i] = (i < 2) ? X0Z(centeredPoints[i]) : X0Z(centeredPoints[3 - i]) + _0Y0(height);
                    break;
                case RIGHT:
                    curPoints[i] = (i >= 2) ? X0Z(centeredPoints[i]) : X0Z(centeredPoints[3 - i]) + _0Y0(height);
                    break;
                case FORWARD:
                    curPoints[i] = (i % 3 != 0) 
                        ? X0Z(centeredPoints[i])
                        : ((i == 0) ? X0Z(centeredPoints[1]) + _0Y0(height) : X0Z(centeredPoints[2]) + _0Y0(height));
                    break;
                case BACK:
                    curPoints[i] = (i % 3 == 0) 
                        ? X0Z(centeredPoints[i])
                        : ((i == 1) ? X0Z(centeredPoints[0]) + _0Y0(height) : X0Z(centeredPoints[3]) + _0Y0(height));
                    break;
                case TOP:
                    curPoints[i] = X0Z(centeredPoints[i]) + _0Y0(height);
                    break;
                case DOWN:
                    curPoints[i] = X0Z(centeredPoints[i]);
                    break;
            }
        }

        std::array<glm::vec3, 4> curTCs;
        for (int i = 0; i < 4; i++)
        {
            switch (side)
            {
                case LEFT: case RIGHT: case FORWARD: case BACK:
                    curTCs[i] = glm::vec3{glm::length(XZ(curPoints[i]) - XZ(curPoints[0])), curPoints[i].y, -1.f};
                    curTCs[i].x /= curInfo.GetRealSize().x;
                    curTCs[i].y /= curInfo.GetRealSize().y;
                    break;
                case TOP:
                    curTCs[i] = XY0(centeredPoints[i]) - XY0(centeredPoints[0]);
                    curTCs[i] *= 1.f / curInfo.roofRealSize;
                    curTCs[i].z = 1;
                    break;
                case DOWN:
                    curTCs[i] = glm::vec3{0.5f, 0.5f, -1};
                    break;
            }
        }

        glm::vec3 normal;
        switch (side)
        {
            case LEFT:
                normal = glm::vec3{-1, 0, 0};
                break;
            case RIGHT:
                normal = glm::vec3{1, 0, 0};
                break;
            case FORWARD:
                normal = glm::vec3{0, 0, 1};
                break;
            case BACK:
                normal = glm::vec3{0, 0, -1};
                break;
            case TOP:
                normal = glm::vec3{0, 1, 0};
                break;
            case DOWN:
                normal = glm::vec3{0, -1, 0};
                break;
        }
        
        for (int i = 1; i < centeredPoints.size() - 1; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                glm::vec3 pos = (j == 2)
                    ? curPoints.at(0)
                    : curPoints.at(i + j);
                glm::vec3 tc = (j == 2)
                    ? curTCs.at(0)
                    : curTCs.at(i + j);

                memcpy(ptr, glm::value_ptr(pos), sizeof(pos));
                ptr += (sizeof(pos) / sizeof(pos.x));

                memcpy(ptr, glm::value_ptr(normal), sizeof(normal));
                ptr += (sizeof(normal) / sizeof(normal.x));

                memcpy(ptr, glm::value_ptr(tc), sizeof(tc));
                ptr += (sizeof(tc) / sizeof(tc.x));
            }   
        }
    }
    

    
    vertexesNum = _vertexNum;
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(
        GL_ARRAY_BUFFER, 
        vertexesNum * VERTEX_SIZE_IN_FLOAT * FLOAT_SIZE_IN_BYTE, 
        generalBuffer, 
        GL_STATIC_DRAW);

    glVertexAttribPointer(
        0, 
        POSITION_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)0
    );
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(
        1, 
        NORMAL_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)(POSITION_SIZE_IN_FLOAT * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(
        2, 
        TC_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)((POSITION_SIZE_IN_FLOAT + NORMAL_SIZE_IN_FLOAT) * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    delete [] generalBuffer;

    // Texture
    buildingTex = ImageLoader::StbLoadImage(std::string(ASSETS_FOLDER) + curInfo.path);
    buildingSpecTex = ImageLoader::StbLoadImage(std::string(ASSETS_FOLDER) + curInfo.specularPath);
    roofTex = ImageLoader::StbLoadImage(std::string(ASSETS_FOLDER) + curInfo.roofPath);

    glBindTexture(GL_TEXTURE_2D, buildingTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, roofTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindTexture(GL_TEXTURE_2D, buildingSpecTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    checkForGlErrors("Building creation");
}

// Building::Building(Polygon polygon, float h)
// {
//     InitBuildingTypes();

//     BuildingInfo curInfo = buildingTypes.at(BUILDING_TYPE::SCYSCRAPER);

//     //preload
//     m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer("buildingShader");

//     UseShader();

//     Scene3DObjectDataHolder::VerifyUnoformPos(viewProjMatrix_Uniform, "view_proj_matrix");
//     Scene3DObjectDataHolder::VerifyUnoformPos(sunPos_Uniform, "sun_pos");
//     Scene3DObjectDataHolder::VerifyUnoformPos(sunColor_Uniform, "sun_color");
//     Scene3DObjectDataHolder::VerifyUnoformPos(cameraPos_Uniform, "camera_position");
//     Scene3DObjectDataHolder::VerifyUnoformPos(buildingTex_Uniform, "building_texture");
//     Scene3DObjectDataHolder::VerifyUnoformPos(roofTex_Uniform, "roof_texture");
//     Scene3DObjectDataHolder::VerifyUnoformPos(buildingSpecTex_Uniform, "specular_texture");


//     //Vertex data
//     const int POLYGON_SIZE_IN_VERTEX = 3;
//     const int POSITION_SIZE_IN_FLOAT = 3;
//     const int NORMAL_SIZE_IN_FLOAT = 3;
//     const int TC_SIZE_IN_FLOAT = 3;
//     const int VERTEX_SIZE_IN_FLOAT = POSITION_SIZE_IN_FLOAT + NORMAL_SIZE_IN_FLOAT + TC_SIZE_IN_FLOAT;

//     int SIDE_SIZE_IN_POLYGON = polygon.points.size() - 2;
//     int SIDE_COUNT = 1;
//     int POLYGON_COUNT = SIDE_SIZE_IN_POLYGON * SIDE_COUNT;

//     int _vertexNum = POLYGON_COUNT * POLYGON_SIZE_IN_VERTEX;
//     GLfloat* generalBuffer = new GLfloat[_vertexNum * VERTEX_SIZE_IN_FLOAT];
//     const int FLOAT_SIZE_IN_BYTE = sizeof(generalBuffer[0]);

    
//     GLfloat* ptr = generalBuffer;
//     for (int i = 1; i < polygon.points.size() - 1; i++)
//     {
//         for (int j = 0; j < 3; j++)
//         {
//             glm::vec3 pos = X0Z((j == 2)
//                 ? polygon.points.at(0)
//                 : polygon.points.at(i + j)) + _0Y0(h);
//             glm::vec3 tc = glm::vec3{fabs(sinf(h * 123.23)), fabs(cosf(h * 123.23)), -1};
//             glm::vec3 normal{0, 1, 0};

//             memcpy(ptr, glm::value_ptr(pos), sizeof(pos));
//             ptr += (sizeof(pos) / sizeof(pos.x));

//             memcpy(ptr, glm::value_ptr(normal), sizeof(normal));
//             ptr += (sizeof(normal) / sizeof(normal.x));

//             memcpy(ptr, glm::value_ptr(tc), sizeof(tc));
//             ptr += (sizeof(tc) / sizeof(tc.x));
//         }   
//     }
    

//     vertexesNum = _vertexNum;
//     glBindVertexArray(VAO);
//     glBindBuffer(GL_ARRAY_BUFFER, VBO);
//     glBufferData(
//         GL_ARRAY_BUFFER, 
//         vertexesNum * VERTEX_SIZE_IN_FLOAT * FLOAT_SIZE_IN_BYTE, 
//         generalBuffer, 
//         GL_STATIC_DRAW);

//     glVertexAttribPointer(
//         0, 
//         POSITION_SIZE_IN_FLOAT, 
//         GL_FLOAT, 
//         GL_FALSE, 
//         VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
//         (GLvoid*)0
//     );
//     glEnableVertexAttribArray(0);

//     glVertexAttribPointer(
//         1, 
//         NORMAL_SIZE_IN_FLOAT, 
//         GL_FLOAT, 
//         GL_FALSE, 
//         VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
//         (GLvoid*)(POSITION_SIZE_IN_FLOAT * sizeof(GLfloat))
//     );
//     glEnableVertexAttribArray(1);

//     glVertexAttribPointer(
//         2, 
//         TC_SIZE_IN_FLOAT, 
//         GL_FLOAT, 
//         GL_FALSE, 
//         VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
//         (GLvoid*)((POSITION_SIZE_IN_FLOAT + NORMAL_SIZE_IN_FLOAT) * sizeof(GLfloat))
//     );
//     glEnableVertexAttribArray(2);

//     glBindBuffer(GL_ARRAY_BUFFER, 0);
//     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
//     glBindVertexArray(0);

//     delete [] generalBuffer;

//     // Texture
//     buildingTex = ImageLoader::StbLoadImage(std::string(ASSETS_FOLDER) + curInfo.path);
//     buildingSpecTex = ImageLoader::StbLoadImage(std::string(ASSETS_FOLDER) + curInfo.specularPath);
//     roofTex = ImageLoader::StbLoadImage(std::string(ASSETS_FOLDER) + curInfo.roofPath);

//     glBindTexture(GL_TEXTURE_2D, buildingTex);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//     glBindTexture(GL_TEXTURE_2D, 0);

//     glBindTexture(GL_TEXTURE_2D, roofTex);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//     glBindTexture(GL_TEXTURE_2D, 0);

//     glBindTexture(GL_TEXTURE_2D, buildingSpecTex);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//     glBindTexture(GL_TEXTURE_2D, 0);

//     checkForGlErrors("Building creation");
// }

void Building::Render(RENDER_MODE mode)
{
    //SET STATE
    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);

    //UNIFORMS
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    glm::mat4 viewProj = cam->GetViewProjectionMatrix();

    UseShader();
    
    LoadInShader_Scene3DObject();
    glUniformMatrix4fv(viewProjMatrix_Uniform, 1, GL_FALSE, glm::value_ptr(viewProj));

    Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
            2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];
    TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
    SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();
    glUniform3fv(sunPos_Uniform, 1, glm::value_ptr(sunTransform->position));
    glUniform3fv(sunColor_Uniform, 1, glm::value_ptr(sunComp->color));

    glUniform3fv(cameraPos_Uniform, 1, glm::value_ptr(cam->GetPosition()));

    //Textures
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, buildingTex);
    glUniform1i(buildingTex_Uniform, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, roofTex);
    glUniform1i(roofTex_Uniform, 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, buildingSpecTex);
    glUniform1i(buildingSpecTex_Uniform, 2);

    //DRAW CALL
    RenderableObjectDataHolder::Render(mode);
}

glm::vec2 Polygon::GetPoint(int index)
{
    if (index == -1)
        return points[points.size() - 1];
    else if (index >= 0 && index < points.size())
        return points[index];
    else if (index == points.size())
        return points[0];
    else
    {
        debug ("WRONG INDEX!", index, points.size());
        throw std::exception{};
    }
}

Rect Polygon::Outscribe()
{
    std::array<std::array<float, 2>, 2> bounds;
    for (int i = 0; i < points.size(); i++)
    {
        if (i == 0)
        {
            bounds[0][0] = bounds[0][1] = points[i].x;
            bounds[1][0] = bounds[1][1] = points[i].y;
        }
        else
        {
            bounds[0][0] = std::min(bounds[0][0], points[i].x);
            bounds[1][0] = std::min(bounds[1][0], points[i].y);
            bounds[0][1] = std::max(bounds[0][1], points[i].x);
            bounds[1][1] = std::max(bounds[1][1], points[i].y);
        }
    }

    glm::vec2 rectCenter{0.5f * (bounds[0][0] + bounds[0][1]), 0.5f * (bounds[1][0] + bounds[1][1])};
    glm::vec2 size{fabsf(bounds[0][0] - bounds[0][1]),fabsf(bounds[1][0] - bounds[1][1])};
    Rect res(rectCenter, size, 0.f);
    return res;
}

std::unique_ptr<Rect> Polygon::Insribe()
{
    const int FIND_SQUARE_POINTS = Landscape::BuildingsGenerator::INSCRIBE_FIND_SQUARE_POINTS;

    const int SIZE_ATTEMPTS = Landscape::BuildingsGenerator::INSCRIBE_SIZE_ATTEMPTS;
    const int OVERAL_ATTEMPTS_AMOUNT = Landscape::BuildingsGenerator::INSCRIBE_OVERAL_ATTEMPTS_AMOUNT;
    const float MAX_SIDE_RATIO = Landscape::BuildingsGenerator::INSCRIBE_MAX_SIDE_RATIO;
    const float MIN_AREA_REAL = Landscape::BuildingsGenerator::INSCRIBE_MIN_AREA_REAL;
    const float MAX_AREA_START_PART = Landscape::BuildingsGenerator::INSCRIBE_MAX_AREA_START_PART;

    //derived
    const int ATTEMPTS_PER_SIZE = OVERAL_ATTEMPTS_AMOUNT / SIZE_ATTEMPTS;

    if (SIZE_ATTEMPTS < 2)
    {
        debug("TOO LESS SIZES!");
        throw std::exception{};
    }

    Rect bBox = Outscribe();

    float areaPart = 0;
    for (int i = 0; i < FIND_SQUARE_POINTS; i++)
    {
        if (IsPointInside(bBox.RandomPoint(0)))
            areaPart += 1.f / FIND_SQUARE_POINTS;
    }

    auto TryInscribeRect =

        [this, &bBox,
            ATTEMPTS_PER_SIZE, MAX_SIDE_RATIO
        ]
        (float _area)
        {
            std::unique_ptr<Rect> result(nullptr);
            for (int i = 0; i < ATTEMPTS_PER_SIZE; i++)
            {
                float division = LERP(MAX_SIDE_RATIO, 1.0f - MAX_SIDE_RATIO, rnd());
                float b = sqrtf((_area / division) - _area);
                glm::vec2 size{_area / b, b};

                float shift = std::min(size.x, size.y) / 2;
                glm::vec2 center = bBox.RandomPoint(shift);
                if (!IsPointInside(center))
                    continue;
                
                float angle = rnd() * M_PI;
                
                Rect cur(center, size, angle);
                if (!CheckRectIntersect(cur))
                {
                    result = std::make_unique<Rect>(cur);
                    break;
                }
            }
            return result;
        };

    std::array<float, 2> areas{MIN_AREA_REAL, MAX_AREA_START_PART * areaPart * bBox.Area()};
    std::unique_ptr<Rect> minRect = TryInscribeRect(areas[0]);
    std::unique_ptr<Rect> maxRect = TryInscribeRect(areas[1]);

    if (maxRect)
    {
        return maxRect;
    }

    if (!minRect)
        return minRect;

    for (int i = 2; i < SIZE_ATTEMPTS; i++)
    {
        float curArea = 0.5f * (areas[0] + areas[1]);
        std::unique_ptr<Rect> curRect = TryInscribeRect(curArea);
        if (curRect)
        {
            minRect = std::move(curRect);
            areas[0] = curArea;
        }
        else
        {
            areas[1] = curArea;
        }
    }
    
    return minRect;
}

bool Polygon::IsPointInside(glm::vec2 p)
{
    int recheck = 3;
    for (int j = 0; j < recheck; j++)
    {
        int counter = 0;
        float angleRad = rnd() * M_PI * 2;
        glm::vec2 dir = cosf(angleRad) * glm::vec2(1, 0) + sinf(angleRad) * glm::vec2(0, 1);
        glm::vec2 newP = p + dir * 100000.f;
        for (int i = 0; i < points.size(); i++)
        {
            if (IntersectSegments(p, newP, GetPoint(i), GetPoint(i + 1)).second)
                counter++;
        }

        if (counter % 2 == 0)
            return false;
    }
    return true;
}

glm::vec2 Rect::RandomPoint(float shift)
{
    glm::vec2 forward = glm::vec2{cosf(angleRad), sinf(angleRad)};
    glm::vec2 left = LeftNormal(forward);
    if (size.x / 2 < shift || size.y / 2 < shift)
        return center;

    glm::vec2 innerSize = size;
    innerSize.x -= shift * 2;
    innerSize.y -= shift * 2;
    
    glm::vec2 s{LERP(-innerSize.x * 0.5f, innerSize.x * 0.5f, rnd()), LERP(-innerSize.y * 0.5f, innerSize.y * 0.5f, rnd())};
    return (center + forward * s.x + left * s.y);
}

Rect::Rect(glm::vec2 p, glm::vec2 s, float a)
{
    center = p;
    size = s;
    angleRad = a; 

    glm::vec2 forward = glm::vec2{cosf(angleRad), sinf(angleRad)};
    glm::vec2 left = LeftNormal(forward);

    for (int i = 0; i < 4; i++)
    {
        float xFactor = (i < 2) ? -1 : 1;
        float yFactor = (i % 3 == 0) ? -1 : 1;
        points.push_back(glm::vec2{});
        points[i] = center;
        points[i] += xFactor * forward * size[0] * 0.5f;
        points[i] += yFactor * left * size[1] * 0.5f;
    }
}

std::vector<glm::vec2> Rect::GetPoints()
{
    return points;
}

Polygon::Polygon(Rect r)
{
    points = r.GetPoints();
}

bool Polygon::CheckRectIntersect(Rect& r)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < points.size(); j++)
        {
            if (IntersectSegments(r.GetPoint(i), r.GetPoint(i + 1), GetPoint(j), GetPoint(j+1)).second)
                return true;
        }
    }
    return false;
}

float Rect::Area()
{
    return size.x * size.y;
}

Basement::Basement(const Rect& r) : Rect(r)
{
    height = 0.f;
}


void Building::InitBuildingTypes()
{
    if (buildingTypes.size() != 0)
        return;

    BuildingInfo scyscraper("Textures/Buildings/building_one.png");
    scyscraper.realHeight = 3.f;
    scyscraper.roofPath = "Textures/Buildings/building_one_roof.jpg";
    scyscraper.specularPath = "Textures/Buildings/building_one_specular.png";
    scyscraper.roofRealSize = 1.5f;

    BuildingInfo panel("Textures/Buildings/building_two.png");
    panel.realHeight = 3.f;
    panel.roofPath = "Textures/Buildings/building_two_roof.png";
    panel.specularPath = "Textures/Buildings/building_two_specular.png";
    panel.roofRealSize = 8.f;

    buildingTypes = std::map<BUILDING_TYPE, BuildingInfo>
    {
        {BUILDING_TYPE::SCYSCRAPER, scyscraper},
        {BUILDING_TYPE::HIGH_RISE, panel}
    };
}

BuildingInfo::BuildingInfo(std::string _path) : path{_path}
{
    int imageW, imageH;
    unsigned char* image = stbi_load((std::string(ASSETS_FOLDER) + path).c_str(), &imageW, &imageH, 0, 4);

    WidthToHeight = (float)imageW / imageH;

    stbi_image_free(image);
}

glm::vec2 BuildingInfo::GetRealSize()
{
    return glm::vec2{realHeight * WidthToHeight, realHeight};
}