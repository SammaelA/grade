#include "LandGeneratorClass.h"
#include "third_party/stb_image.h"
#include "cities_generator/ImageLoader.h"
#include "cities_generator/ECSClass.h"
#include "common_utils/perlin.h"
#include "cities_generator/shaderClass.h"
#include "cities_generator/Renderer.h"
#include "cities_generator/GameManager.h"
#include <cmath>
#include "../EventClasses.h"
#include <algorithm>
//#include "../../Libraries/eigen/Eigen/Householder"
//#include "../../Libraries/eigen/Eigen/Dense"
#include "WaterPlate.h"
#include "BreadthSearch.h"
#include <unordered_map>
#include <memory>

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;

Landscape::HillsGenerator::HillsGenerator()
{
    m_fromDegrees = 25.0;
    m_toDegrees = 40.0;
    m_minRockRatio = 0.15;
    m_maxRockRatio = 1.f;

    m_bordersCos[0] = cosf(glm::radians(m_toDegrees));
    m_bordersCos[1] = cosf(glm::radians(m_fromDegrees));
}

float Landscape::HillsGenerator::GetGrassRatioByCos(float flatAngleCos)
{
    if (flatAngleCos < m_bordersCos[0])
        return 1 - m_maxRockRatio;
    if (flatAngleCos > m_bordersCos[1])
        return 1.0;
    return 1 - LERP(
        m_minRockRatio, 
        m_maxRockRatio, 
        (m_bordersCos[1] - flatAngleCos) / (m_bordersCos[1] - m_bordersCos[0])
    );
}

void Landscape::ReadNormal(int x, int y, GLfloat output[])
{
    glm::vec2 shifts[] = {
        {0, 1},
        {1, 0},
        {0, -1},
        {-1, 0}
    };
    glm::vec3 points[4];
    float defaultHeight = m_heightmap->Get(x, y);
    for (int i = 0; i < 4; i++)
    {
        glm::vec2 realCoords = PointMapIntToReal(m_heightmap, vec2Int(shifts[i][0], shifts[i][1]));
        points[i] = {
            realCoords.x, 
            m_heightmap->Get(x + shifts[i][0], y + shifts[i][1], true, defaultHeight), 
            realCoords.y
        };
    }

    glm::vec3 res;
    glm::vec3 currentPoint {0, defaultHeight, 0};
    for (int i = 0; i < 4; i++)
    {
        glm::vec3 vec1 = points[i] - currentPoint;
        glm::vec3 vec2 = points[(i + 1) % 4] - currentPoint;
        res += glm::normalize(glm::cross(vec1, vec2)) * glm::vec3(0.25);
    }

    res = glm::normalize(res);
    for (int i = 0; i < 3; i++)
    {
        output[i] = res[i];
    }
}

glm::vec3 Landscape::GetNormal(int x, int y)
{
    glm::vec3 normal;
    GLfloat temp[3];
    ReadNormal(x, y, temp);
    for (int i = 0; i < 3; i++)
    {
        normal[i] = temp[i];
    }
    return normal;
}

template<class T>
glm::vec2 Landscape::PointMapIntToReal(DefaultMap<T>* map, vec2Int coords)
{
    return GetMapCellRealSize(map) * glm::vec2(coords);
}

template<class T>
glm::vec2 Landscape::PointRealToMapNormalized(DefaultMap<T>* map, glm::vec2 coords, bool safe)
{
    return map->PointIntContinuousToNormalized(coords / GetMapCellRealSize(map), safe);
}

template<class T>
glm::vec2 Landscape::PointMapNormalizedToReal(DefaultMap<T>* map, glm::vec2 coords)
{
    return map->PointNormalizedToIntContinuous(coords) * GetMapCellRealSize(map);
}

float Landscape::GetHillsGrassRatio(glm::vec3 pos, glm::vec3 normal)
{
    glm::vec3 upVec {0.0, 1.0, 0.0};
    float flatCos = glm::dot(upVec, normal);
    float slopeRatio = hillsGenerator.GetGrassRatioByCos(flatCos);

    float biomeRatio = 1.f;
    float BIOME_TO_REAL = GetMapCellRealSize(m_biomemap);
    BiomePoint bp = m_biomemap->GetInterpolated_Normalized(
            PointRealToMapNormalized(m_heightmap, XZ(pos))
        )[0];
    biomeRatio = (bp.biome == Biome::MOUNTAINS) ?
        (bp.intensity / HillsGenerator::SEMI_MOUNTAIN_AREA_REAL * BIOME_TO_REAL) : 0.f;
    biomeRatio = clamp(biomeRatio, 0.f, 1.f);
    biomeRatio = LERP(1.f, HillsGenerator::MAX_MOUNTAIN_GRASS_TEXTURE_ALPHA, biomeRatio);

    return slopeRatio * biomeRatio;
}

float Landscape::GetHillsSlopeRatio(glm::vec3 pos, glm::vec3 normal)
{
    return std::max(0.f, 1 - GetHillsGrassRatio(pos, normal));
}

float Landscape::GetWaterTextureRatio(glm::vec3 pos)
{
    static const float SAND_MIN_CELLS = WaterGenerator::SAND_TEXTURE_MIN_TO_SHORE *
        WaterGenerator::SHORE_REAL_SIZE / GetMapCellRealSize(m_biomemap);
    static const float SAND_MAX_CELLS = WaterGenerator::SAND_TEXTURE_MAX_TO_SHORE *
        WaterGenerator::SHORE_REAL_SIZE / GetMapCellRealSize(m_biomemap);
    
    if (!IsRenderingWater())
        return 0;

    if (pos.y < GetWaterLevel())
        return 1;
    else
    {
        glm::vec2 tc = PointRealToMapNormalized(m_heightmap, glm::vec2{pos.x, pos.z});
        SurfacePoint sp = m_surfacemap->GetSurfaceInterpolated_Normalized(tc)[0];
        float ratio = (SAND_MAX_CELLS - sp.intensity) / (SAND_MAX_CELLS - SAND_MIN_CELLS);
        return clamp(ratio, 0.f, 1.f);
    }
}

Landscape::Landscape()
{
    m_realSize[0] = 1000; //FINAL 1000
    m_realSize[1] = 1000; //FINAL 1000   
    // m_realSize[0] = 12.5;
    // m_realSize[1] = 12.5;
    // m_realSize[0] = 2;
    // m_realSize[1] = 2;

    m_heightmap = nullptr;

    m_grassTexture = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/grassTexture.png"
    );

    m_soilTexture = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/soilTexture.jpg"
    );

    m_sandTexture = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/sandTexture.jpg"
    );
}

Landscape::~Landscape()
{
    delete m_heightmap;
    delete m_savedHeightmap;
    delete m_biomemap;
    delete m_surfacemap;
    delete m_roadScheme;

    glDeleteTextures(1, &m_grassTexture);
    glDeleteTextures(1, &m_soilTexture);
    glDeleteTextures(1, &m_sandTexture);
}

void Landscape::Init()
{
    checkForGlErrors("LANDSCAPE PRE-INIT");
    
    // GenerateHeightmap_DIAMONDSQUARE(true);
    GenerateBiomemap();
    GenerateHeightmap_PERLIN();
    GeneratePopulationMap();
    GenerateHeightmap_PERLIN();
    GenerateRoad();
    GenerateBuildings();
    
    ApplyShoreImmerse();
    ApplyErosionToHeightmap();
    RecreateMesh();
}

void Landscape::RecreateMesh()
{
    debug("Creating landscape mesh");

    if (!m_heightmap)
    {
        debug("NO HEIGHTMAP!");
        throw std::exception{};
    }

    glm::vec2 heightmapSize = m_heightmap->GetSize();

    m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer("landscapeShader");

    int quadAmount = heightmapSize[0] * heightmapSize[1] * 4;
    int triangleAmount = quadAmount * 2;
    const int SIZE_OF_TRIANGLE_IN_VERTEX = 3;
    const int SIZE_OF_VERTEX_IN_GLFLOAT = 8;
    vertexesNum = triangleAmount * SIZE_OF_TRIANGLE_IN_VERTEX;
    GLfloat* vertexes = (GLfloat*) malloc(
        triangleAmount * SIZE_OF_TRIANGLE_IN_VERTEX  * SIZE_OF_VERTEX_IN_GLFLOAT * sizeof(GLfloat)
    );
    GLfloat* currentVertexesPtr = vertexes;
    
    unsigned indexOrder[] = {0, 1, 2, 1, 3, 2};
    int quadVertexesShift[4][2] = 
    {
        {0 ,0},
        {1, 0},
        {0, 1},
        {1, 1}
    };
    const int QUAD_SIZE_IN_GLFLOAT = 2 * 3 * SIZE_OF_VERTEX_IN_GLFLOAT;

    for (int i = 0; i < heightmapSize[0] * 2; i++)
    {
        for (int j = 0; j < heightmapSize[1] * 2; j++)
        {
            int curPos[] = {i - heightmapSize[0], j - heightmapSize[1]};

            GLfloat currentQuadHeight[4];
            for (int k = 0; k < 4; k++)
            {
                currentQuadHeight[k] = m_heightmap->Get(
                    curPos[0] + quadVertexesShift[k][0],
                    curPos[1] + quadVertexesShift[k][1]
                );
            }

            GLfloat currentQuadNormal[3 * 4];
            for (int k = 0; k < 4; k++)
            {
                ReadNormal(
                    curPos[0] + quadVertexesShift[k][0],
                    curPos[1] + quadVertexesShift[k][1],
                    currentQuadNormal + k * 3
                );
            }

            GLfloat* quad = (GLfloat*) malloc(QUAD_SIZE_IN_GLFLOAT * sizeof(GLfloat));
            for (int k = 0; k < sizeof(indexOrder)  / sizeof(indexOrder[0]); k++)
            {
                const int SHIFT = SIZE_OF_VERTEX_IN_GLFLOAT;

                vec2Int currerntMapTC
                {
                    (curPos[0] + quadVertexesShift[indexOrder[k]][0]),
                    (curPos[1]  + quadVertexesShift[indexOrder[k]][1])
                };
                glm::vec2 currentMapRealPos = PointMapIntToReal(m_heightmap,currerntMapTC);

                glm::vec3 currentPointPos
                {
                    currentMapRealPos[0],
                    currentQuadHeight[indexOrder[k]],
                    currentMapRealPos[1]
                };
                quad[SHIFT * k] = currentPointPos.x;
                quad[SHIFT * k + 1] = currentPointPos.y;
                quad[SHIFT * k + 2] = currentPointPos.z;

                memcpy(
                    quad + SHIFT * k + 3, 
                    currentQuadNormal + indexOrder[k] * 3, 
                    sizeof(GLfloat) * 3
                );

                glm::vec3 currentPointNormal = {
                    currentQuadNormal[indexOrder[k] * 3],
                    currentQuadNormal[indexOrder[k] * 3 + 1],
                    currentQuadNormal[indexOrder[k] * 3 + 2]
                };

                quad[SHIFT * k + 6] = GetHillsSlopeRatio(currentPointPos, currentPointNormal);
                // quad[SHIFT * k + 6] = 1.f - GetPopulationRatioByRealPos(XZ(currentPointPos));
                quad[SHIFT * k + 7] = GetWaterTextureRatio(currentPointPos);
            }
            memcpy(currentVertexesPtr, quad, sizeof(GLfloat) * QUAD_SIZE_IN_GLFLOAT);
            currentVertexesPtr += QUAD_SIZE_IN_GLFLOAT;
            free(quad);
        }
    }

    const int SHIFT = SIZE_OF_VERTEX_IN_GLFLOAT;
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 
                SIZE_OF_VERTEX_IN_GLFLOAT * vertexesNum * sizeof(vertexes[0]), 
                vertexes, 
                GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, SHIFT * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, SHIFT * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, SHIFT * sizeof(GLfloat), (GLvoid*)(6 * sizeof(GLfloat)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, SHIFT * sizeof(GLfloat), (GLvoid*)(7 * sizeof(GLfloat)));
    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    free(vertexes);

    debug("Done!");
}

void Landscape::GenerateHeightmap_PERLIN()
{
    vec2Int heightMapSize = vec2Int(glm::vec2{m_realSize[0], m_realSize[1]} * HeightMap::DENSITY_PER_REAL);
    m_heightmap = new HeightMap(heightMapSize[0], heightMapSize[1]);
    
    noise::module::Perlin hillsPerlin;
    hillsPerlin.SetFrequency(HillsGenerator::FREQUENCY);
    hillsPerlin.SetOctaveCount(HillsGenerator::OCTAVE_COUNT);
    hillsPerlin.SetPersistence(HillsGenerator::PERSISTANCE);
    hillsPerlin.SetLacunarity(HillsGenerator::LACUNARITY);
    
    noise::module::Perlin plainPerlin;
    plainPerlin.SetFrequency(PlainGenerator::FREQUENCY);
    plainPerlin.SetOctaveCount(PlainGenerator::OCTAVE_COUNT);
    plainPerlin.SetPersistence(PlainGenerator::PERSISTANCE);
    plainPerlin.SetLacunarity(PlainGenerator::LACUNARITY);

    vec2Int diamondSquareShiftSize = vec2Int(
        PlainGenerator::DIAMOND_SQUARE_HELP_MAP_SIZE, PlainGenerator::DIAMOND_SQUARE_HELP_MAP_SIZE
    );
    m_savedHeightmap = new HeightMap(diamondSquareShiftSize[0], diamondSquareShiftSize[1]);
    SwapHeightMaps(false);
    GenerateHeightmap_DIAMONDSQUARE(true);
    SwapHeightMaps(false);

    const float BIOME_TO_REAL = GetMapCellRealSize(m_biomemap);

    for (int i = -m_heightmap->m_width; i <= m_heightmap->m_width; i++)
    {
        for (int j = -m_heightmap->m_height; j <= m_heightmap->m_height; j++)
        {
            glm::vec2 curPos = PointMapIntToReal(m_heightmap, vec2Int(i, j));
            glm::vec2 coords = m_heightmap->PointIntToNormalized(vec2Int{i, j});

            std::vector<BiomePoint> surroundingPoints = m_biomemap->GetInterpolated_Normalized(coords);
            BiomePoint bPoint = surroundingPoints[0];
            for (BiomePoint bp : surroundingPoints)
            {
                if (bp.intensity > bPoint.intensity)
                    bPoint = bp;
            }
            for (BiomePoint bp : surroundingPoints)
            {
                if (bPoint.biome != bp.biome)
                    bPoint.intensity -= bp.intensity;
            }
            // bPoint.intensity = clamp(bPoint.intensity, 0.f, bPoint.intensity);
            bPoint.intensity = std::max(0.f, bPoint.intensity);

            float height = 0;
            if (bPoint.biome == Biome::PLAIN || 
                     bPoint.biome == Biome::WATER || 
                     bPoint.biome == Biome::MOUNTAINS ||
                     bPoint.biome == Biome::CITY)
            {
                float mountainRatio = (bPoint.biome == Biome::MOUNTAINS) ?
                    (bPoint.intensity / HillsGenerator::SEMI_MOUNTAIN_AREA_REAL * BIOME_TO_REAL) : 0.f;

                float cityRatio = (bPoint.biome == Biome::CITY) ?
                    clamp((bPoint.intensity / PlainGenerator::SEMI_CITY_AREA_REAL * BIOME_TO_REAL), 0.f, 1.f) : 0.f;
                
                if (mountainRatio < 1.0f)
                {
                    height = (float) plainPerlin.GetValue(curPos[0], 0, curPos[1]);
                    height = height * 0.5 + 0.5;
                    height *= PlainGenerator::MAX_HEIGHT;
                    height += m_savedHeightmap->GetSafe_EdgeMirrored(vec2Int(i, j)) * PlainGenerator::MAX_DISTORTION;
                    
                    float citySupress = PlainGenerator::CITY_HEIGHT_CHANGE_SUPRESSION_FACTOR;
                    height = LERP(height, height * citySupress, cityRatio);
                }

                if (bPoint.biome == Biome::MOUNTAINS)
                {
                    mountainRatio = clamp(mountainRatio, 0.f, 1.f);
                    float mountainHeight = 0.f;
                    mountainHeight = (float) hillsPerlin.GetValue(curPos[0], 0, curPos[1]);
                    mountainHeight = mountainHeight * 0.5 + 0.5;
                    mountainHeight *= HillsGenerator::MAX_HEIGHT;
                    height = LERP(height, mountainHeight, mountainRatio);
                }
            }
            else
            {
                debug("INCORRECT BIOMEMAP POINT!", (int)bPoint.biome);
                throw std::exception{};
            }
            m_heightmap->Set(i,j,height);
        }
    }

    m_heightmap->SaveAsTexture();
}

void Landscape::GenerateHeightmap_DIAMONDSQUARE(bool relative)
{
    //Preparing heightmap
    int iterCount;
    if (!relative)
    {
        if (m_realSize[0] != m_realSize[1])
        {
            debug("LANDSCAPE MUST BE SQUARE");
            throw std::exception{};
        }

        iterCount = (int)log2(m_realSize[0] * HeightMap::DENSITY_PER_REAL) + 1;
    }
    else
    {
        if (m_heightmap->m_width != m_heightmap->m_height)
        {
            debug("RELATIVE LANDSCAPE MUST BE SQUARE");
            throw std::exception{};
        }
        iterCount = (int)log2(m_heightmap->m_width) + 1;
    }
    int landscapeSize[2];
    landscapeSize[0] = landscapeSize[1] = pow(2, iterCount - 1);
    m_heightmap = new HeightMap(landscapeSize[0], landscapeSize[1]);
    m_heightmap->Set(-m_heightmap->m_width,-m_heightmap->m_height,0);
    m_heightmap->Set(m_heightmap->m_width,-m_heightmap->m_height,0);
    m_heightmap->Set(-m_heightmap->m_width,m_heightmap->m_height,0);
    m_heightmap->Set(m_heightmap->m_width,m_heightmap->m_height,0);


    //Generating
    const float spread = PlainGenerator::DIAMOND_SQUARE_SPREAD;
    glm::vec2 heightmapSizeVec = m_heightmap->GetSize();
    int heightmapSize[] = {(int)heightmapSizeVec[0], (int)heightmapSizeVec[1]};
    {
        if (heightmapSize[0] != heightmapSize[1] || heightmapSize[0] == 0)
        {
            debug("WRONG SIZE!");
            throw std::exception{};
        }
        int temp = heightmapSize[0];
        for (int i = 0; i < iterCount - 1; i++)
        {
            if (temp % 2 != 0)
            {
                debug("WRONG SIZE! [2]");
                throw std::exception{};
            }
            temp /= 2;
        }
    }
    
    
    int shiftsSquare[][2] = 
    {
        {0,0},
        {1,0},
        {1,1},
        {0,1}
    };
    int shiftsDiamond[][2] = 
    {
        {-1,0},
        {0,1},
        {1,0},
        {0,-1}
    };

    while (iterCount != 0)
    {
        int curStep = pow(2,iterCount);

        //Square iteration
        float squareDistance;
        if(!relative)
            squareDistance = GetMapCellRealSize(m_heightmap) * sqrt(2) * pow(2, iterCount - 1);
        else
            squareDistance = pow(2, iterCount - 1) / (2 * landscapeSize[0]) * (sqrt(2) / (0.5 + sqrt(2)));

        for (int i = -heightmapSize[0]; i < heightmapSize[0]; i+=curStep)
        {
            for (int j = -heightmapSize[1]; j < heightmapSize[1]; j+=curStep)
            {
                float averageHeight = 0;
                for (int k = 0; k < 4; k++)
                {
                    averageHeight += 0.25 * 
                        m_heightmap->Get(
                            i + shiftsSquare[k][0] * curStep, 
                            j+ shiftsSquare[k][1] * curStep
                        );
                }
                float rndDelta = 2 * (((float)random() / RAND_MAX) - 0.5);
                if (!relative)
                    rndDelta *= squareDistance * spread;
                else
                    rndDelta *= squareDistance / 0.5;
                
                m_heightmap->Set(
                    i + curStep / 2, 
                    j + curStep / 2, 
                    averageHeight + rndDelta
                );
            }
        }

        //Diamond iteration
        int shiftCounter = 0;
        float diamondDistance;
        if(!relative)
            diamondDistance = GetMapCellRealSize(m_heightmap) * pow(2, iterCount - 1);
        else
            diamondDistance = pow(2, iterCount - 1) / (2 * sqrt(2) * landscapeSize[0]) * (sqrt(2) / (0.5 + sqrt(2)));
        
        for (int i = -heightmapSize[0]; i <= heightmapSize[0]; i+=pow(2, iterCount - 1))
        {
            shiftCounter = 1 - shiftCounter;
            for (int j = -heightmapSize[1] + shiftCounter * pow(2, iterCount - 1);
                j <= heightmapSize[1];
                j+=curStep)
            {
                float averageHeight = 0;
                for (int k = 0; k < 4; k++)
                {
                    averageHeight += 0.25 * 
                        m_heightmap->GetSafe_EdgeMirrored(
                            i + shiftsDiamond[k][0] * pow(2, iterCount - 1),
                            j + shiftsDiamond[k][1] * pow(2, iterCount - 1)
                        );
                }
                float rndDelta = 2 * (((float)random() / RAND_MAX) - 0.5);
                if (!relative)
                    rndDelta *= diamondDistance * spread;
                else
                    rndDelta *= diamondDistance / 0.5;

                m_heightmap->Set(i, j, averageHeight + rndDelta);
            }
        }
        iterCount--;
    }
    m_heightmap->SaveAsTexture();
}

void Landscape::Render(RENDER_MODE rmode)
{
    UseShader();

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    // glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // GLint TransformUniformPos;
    // bool firstSprite = true;
                                            
    glm::mat4 viewProj = glm::mat4(1.0f);
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    
    viewProj = cam->GetViewProjectionMatrix();
             
    glm::mat4 transform = glm::mat4(1.0f);
    
    Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
            2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];

    TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
    SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();

    GLint TransformUniformPos, ViewProjMatrixUniformPos, GlobalIlluminationPointUniformPos,
            GlobalIlluminationColorUniformPos, diffuseColorPos, isIsotropicColorPos,
            isDepthOnlyPos, GlobalIlluminationViewProjUniformPos, CameraPositionPos,
            GrassScalePos;
    
    TransformUniformPos = glGetUniformLocation(m_shader->Program, "transform_matrix");
    ViewProjMatrixUniformPos = glGetUniformLocation(m_shader->Program, "view_proj_matrix");
    GlobalIlluminationPointUniformPos = glGetUniformLocation(m_shader->Program, "sun_pos");
    GlobalIlluminationColorUniformPos = glGetUniformLocation(m_shader->Program, "sun_color"); 
    GlobalIlluminationViewProjUniformPos = glGetUniformLocation(m_shader->Program, "sun_view_proj_matrix");
    // diffuseColorPos = glGetUniformLocation(m_shader->Program, "diffuse_color"); 
    // isIsotropicColorPos = glGetUniformLocation(m_shader->Program, "is_isotropic_color"); 
    // isDepthOnlyPos = glGetUniformLocation(m_shader->Program, "is_depth_only"); 
    CameraPositionPos = glGetUniformLocation(m_shader->Program, "camera_position"); 

    glm::mat4 sunProjM = sunComp->projectionMatrix;
    glm::mat4 sunVeiwM = sunTransform->GetViewMatrixFromThisPoint();
    glm::mat4 sunViewProj = sunProjM * sunVeiwM;

    glUniformMatrix4fv(TransformUniformPos, 1, GL_FALSE, glm::value_ptr(transform));
    glUniformMatrix4fv(ViewProjMatrixUniformPos, 1, GL_FALSE, glm::value_ptr(viewProj));
    glUniformMatrix4fv(GlobalIlluminationViewProjUniformPos, 1, GL_FALSE, glm::value_ptr(sunViewProj));
    glUniform3f(GlobalIlluminationPointUniformPos, sunTransform->position.x, sunTransform->position.y, sunTransform->position.z);
    glUniform3f(GlobalIlluminationColorUniformPos, sunComp->color.x, sunComp->color.y, sunComp->color.z);
    glUniform3fv(CameraPositionPos, 1, glm::value_ptr(cam->GetPosition()));
    // glUniform1i(isIsotropicColorPos, !model->isUsingDiffuseTexture);
    // glUniform1i(isDepthOnlyPos, false);
    // glUniform4f(diffuseColorPos, model->diffuseColor.x, model->diffuseColor.y, model->diffuseColor.z, model->diffuseColor.z);
    glUniform1f(glGetUniformLocation(m_shader->Program, "grass_texture_scale"), GRASS_SCALE);
    glUniform1f(glGetUniformLocation(m_shader->Program, "soil_texture_scale"), SOIL_SCALE);
    glUniform1f(glGetUniformLocation(m_shader->Program, "sand_texture_scale"), SAND_SCALE);
    
    //TEXTURES
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, m_grassTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "grass_texture"), 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, m_soilTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "soil_texture"), 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, m_sandTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "sand_texture"), 2);
    // glActiveTexture(GL_TEXTURE1);
    // glBindTexture(GL_TEXTURE_2D, sunComp->shadowTexture);
    // glUniform1i(glGetUniformLocation(m_shader->Program, "shadow_texture"), 1);
    
    //DRAW CALL
    RenderableObjectDataHolder::Render(rmode);
}

void Landscape::RenderMap(MAP map)
{
    switch (map)
    {
        case MAP::NONE:
            break;

        case MAP::HEIGHT:
            Renderer::GetInstance()->RenderTexture(m_heightmap->GetTexture(), 0.0, 1.0, true);
            break;
        
        case MAP::BIOME:
            Renderer::GetInstance()->RenderTexture(m_biomemap->GetTexture(), 0.0, 1.0, true);
            break;

        case MAP::SURFCE:
            Renderer::GetInstance()->RenderTexture(m_surfacemap->GetTexture(), 0.0, 1.0, true);
            break;

        case MAP::ROAD:
            Renderer::GetInstance()->RenderTexture(m_roadScheme->GetTexture(), 0.0, 1.0, true);
            break;

        case MAP::POPULATION:
            Renderer::GetInstance()->RenderTexture(m_populationmap->GetTexture(), 0.0, 1.0, true);
            break;
        
        default:
            debug("WRONG MAP TYPE!");
            throw std::exception{};
    };
}

template<class T>
float Landscape::GetMapCellRealSize(DefaultMap<T>* map)
{
    float realWidthToHeightRatio = m_realSize[0] / m_realSize[1];
    float mapWidthToHeightRatio = map->m_width / map->m_height;
    if (fabsf(realWidthToHeightRatio - mapWidthToHeightRatio) > 0.0001)
    {
        debug("INCORRECT SIZES! [GetMapCellRealSize]");
        throw std::exception{};
    }
    return m_realSize[0] / (2 * map->m_width);
}

void Landscape::SwapHeightMaps(bool recreate)
{
    HeightMap* temp = m_heightmap;
    m_heightmap = m_savedHeightmap;
    m_savedHeightmap = temp;
    if (recreate)
        RecreateMesh();
}

void Landscape::Erosion::CreateDebugDrop(glm::vec3 pos)
{
    Renderer::GetInstance()->CreateDebugSphere(
        0.05, 
        pos,
        glm::vec3(0.1,0.6,0.7),
        "Erosion"
    );
}

void Landscape::Erosion::MoveDebugDropTo(glm::vec3 pos)
{
    Renderer::GetInstance()->DebugSphereMove(pos, "Erosion");
}

void Landscape::Erosion::DebugDropSetColor(glm::vec3 col)
{
    Renderer::GetInstance()->DebugSphereSetColor(col, "Erosion");
}

void Landscape::Erosion::ErosionDroplet::Reset(glm::vec2 pos)
{
    positionReal = pos;
    direction = glm::vec2(0,0);
    velocity = 1;
    sediment = 0;
    water = 1.0;
    stepCounter = 0;
}

Landscape::Erosion::ErosionDroplet::ErosionDroplet()
{
    Reset(glm::vec2{-1.0, -4.0});
}

void Landscape::ApplyErosionToHeightmap()
{
    if (!m_heightmap)
    {
        debug("NO HEIGHTMAP! [ApplyErosionToHeightmap]");
        throw std::exception{};
    }
    m_savedHeightmap = new HeightMap(*m_heightmap);

    // Erosion::CreateDebugDrop(glm::vec3(3, m_heightmap->Get(3, 0), 0));

    const float INERTIA = Erosion::INERTIA;
    // const float MIN_SLOPE = Erosion::MIN_SLOPE;
    const float GRAVITY = Erosion::GRAVITY;
    const float CAPACITY = Erosion::CAPACITY;;
    const int MAX_STEPS = Erosion::MAX_STEPS;
    const float EROSION_RADIUS = Erosion::EROSION_RADIUS;
    const float EVAPORATION = 1.0f - float(pow(Erosion::BASE_EVAPORATION, 1.0 / MAX_STEPS));
    const float STEP_LENGTH = Erosion::STEP_LENGTH;
    const float EROSION = float(pow(Erosion::BASE_EROSION, STEP_LENGTH));
    const float DEPOSITION = float(pow(Erosion::BASE_DEPOSITION, STEP_LENGTH));
    const float NON_EVAPORATE_ANGLE_RADIANS = Erosion::NON_EVAPORATE_ANGLE_DEGREES * M_PI / 180.0;

    Erosion::ErosionDroplet drop;
    HeightMap* HeightMapHolder = new HeightMap(*m_heightmap);

    const float DROPS_COUNT = m_realSize[0] * m_realSize[1] * Erosion::DENSITY;
    for (int d = 0; d < DROPS_COUNT; d++)
    {
        int debugCurrentPart = (int)(((float) d / DROPS_COUNT) * 20);
        int debugPreviousPart = (int)(((float) (d - 1) / DROPS_COUNT) * 20);
        if (debugCurrentPart != debugPreviousPart)
            debug("Eroding in progress:", std::to_string(debugCurrentPart * 5) + "%");

        glm::vec2 startingPos =  glm::vec2{(float) (rand()) / RAND_MAX, (float) (rand()) / RAND_MAX};
        startingPos[0] = LERP(-(m_realSize[0] / 2), (m_realSize[0] / 2), startingPos[0]);
        startingPos[1] = LERP(-(m_realSize[1] / 2), (m_realSize[1] / 2), startingPos[1]);
        if (!IsPointErodable(PointRealToMapNormalized(m_heightmap, startingPos)))
            continue;

        drop.Reset(startingPos);
        while (drop.stepCounter < MAX_STEPS)
        {
            // if (isFollowigDrops)
            // {
            //     float cellPerSecond = 10;
            //     if (!SleepInterval::interval(1000.0 / cellPerSecond, "debugDrop"))
            //         return;
            // }
            glm::vec2 positionTC = m_heightmap->PointNormalizedToIntContinuous(
                PointRealToMapNormalized(m_heightmap,drop.positionReal)
            );
            glm::vec2 positionOld = drop.positionReal;
            glm::vec2 positionOldTC = positionTC; 

            // NOTE: [0] - pos1 [1] - pos2
            std::array<std::array<int, 2>, 2> surroundingRect = HeightMap::GetSurroundingRect(positionTC);
            if (!m_heightmap->IsPointValid(surroundingRect[0][0], surroundingRect[0][1]) ||
                !m_heightmap->IsPointValid(surroundingRect[1][0], surroundingRect[1][1]))
            {
                break;
                // return;
            }
            float heightOld = m_heightmap->GetInterpolated_IntContinuous(m_heightmap->PointNormalizedToIntContinuous(
                PointRealToMapNormalized(m_heightmap, positionOld)
            ));

            // NOTE: (0, 0) -> (0, 1) -> (1, 0) -> (1, 1)
            float rectHeights[4];
            for (int i = 0; i < 4; i++)
            {
                rectHeights[i] = m_heightmap->Get(
                    surroundingRect[(int) (i / 2)][0],
                    surroundingRect[(int) (i % 2)][1]
                );
            }
            
            glm::vec2 interpolateRatio = {
                positionTC.x - surroundingRect[0][0],
                positionTC.y - surroundingRect[0][1]
            };
            glm::vec2 slopeGradient = {
                (rectHeights[2] - rectHeights[0]) * (1 - interpolateRatio[0]) +
                    (rectHeights[3] - rectHeights[1]) * interpolateRatio[0],
                (rectHeights[1] - rectHeights[0]) * (1 - interpolateRatio[1]) +
                    (rectHeights[3] - rectHeights[2]) * interpolateRatio[0]
            };
            slopeGradient = slopeGradient / GetMapCellRealSize(m_heightmap);


            // if (drop.velocity < 0.001)
            // {
            //     drop.velocity = 1.0;
            //     drop.direction = glm::normalize(-slopeGradient);
            // }
            // else
            // {

            // }

            drop.direction = drop.direction * INERTIA - slopeGradient * (1 - INERTIA);
            if (glm::length(drop.direction) < 0.0001)
                drop.direction = glm::vec2((float) rand() / RAND_MAX, (float) rand() / RAND_MAX);
            drop.direction = glm::normalize(drop.direction);

            
            drop.positionReal += drop.direction * GetMapCellRealSize(m_heightmap) * STEP_LENGTH;
            glm::vec2 posNorm = PointRealToMapNormalized(m_heightmap,drop.positionReal);
            if (!m_heightmap->IsPointValid_Normalized(posNorm))
            {
                break;
            }
            positionTC = m_heightmap->PointNormalizedToIntContinuous(posNorm);
            surroundingRect = HeightMap::GetSurroundingRect(positionTC);
            if (!m_heightmap->IsPointValid(surroundingRect[0][0], surroundingRect[0][1]) ||
                !m_heightmap->IsPointValid(surroundingRect[1][0], surroundingRect[1][1]))
            {
                break;
                // return;
            }
            float height = m_heightmap->GetInterpolated_IntContinuous(positionTC);
            float heightDifference = height - heightOld;

            if (heightDifference < 0)
            {
                float slopeAngleRadians = atan2(-heightDifference, GetMapCellRealSize(m_heightmap) * STEP_LENGTH);
                float evaporation = slopeAngleRadians / NON_EVAPORATE_ANGLE_RADIANS;
                evaporation = (evaporation > 1) ? 1 : evaporation;
                evaporation = 1 - sqrt(evaporation);
                evaporation *= EVAPORATION;
                drop.water = drop.water * (1 - evaporation);
            }
            else
            {
                drop.water = drop.water * (1 - EVAPORATION);
            }
            drop.velocity = std::sqrt(std::fabs(drop.velocity * drop.velocity - heightDifference * GRAVITY));
            drop.stepCounter++;

            
            // float capacity = std::max(-heightDifference, MIN_SLOPE * GetMapCellRealSize(m_heightmap));
            // capacity *= drop.velocity * drop.water * CAPACITY;
            float capacity = std::max(-heightDifference, 0.0f) * drop.velocity * drop.water * CAPACITY;
            // debug("CAPACITY", capacity, -heightDifference, drop.velocity, drop.water, drop.sediment);

            if (capacity > drop.sediment && heightDifference < 0)
            {
                // float deltaSediment = std::min((capacity - drop.sediment) * EROSION, -heightDifference);
                float deltaSediment = (capacity - drop.sediment) * EROSION;
                drop.sediment += deltaSediment;
                // debug("Eroded", deltaSediment);
                float weightSum = 0;
                for (int i = (int)ceil(positionOldTC.x - EROSION_RADIUS);
                    i <= (int)floor(positionOldTC.x + EROSION_RADIUS);
                    i++)
                {
                    for (int j = (int)ceil(positionOldTC.y - EROSION_RADIUS);
                        j <= (int)floor(positionOldTC.y + EROSION_RADIUS);
                        j++)
                    {
                        if (!m_heightmap->IsPointValid(i, j))
                            continue;

                        float distanceSqr = pow(positionOldTC.x - i, 2) + pow(positionOldTC.y - j, 2);
                        if (distanceSqr <= pow(EROSION_RADIUS, 2))
                        {
                            weightSum += EROSION_RADIUS - sqrt(distanceSqr);
                        }
                    }
                }
                for (int i = (int)ceil(positionOldTC.x - EROSION_RADIUS);
                    i <= (int)floor(positionOldTC.x + EROSION_RADIUS);
                    i++)
                {
                    for (int j = (int)ceil(positionOldTC.y - EROSION_RADIUS);
                        j <= (int)floor(positionOldTC.y + EROSION_RADIUS);
                        j++)
                    {
                        if (!m_heightmap->IsPointValid(i, j))
                            continue;
                            
                        float distanceSqr = pow(positionOldTC.x - i, 2) + pow(positionOldTC.y - j, 2);
                        if (distanceSqr <= pow(EROSION_RADIUS, 2))
                        {
                            float weightCurrent = EROSION_RADIUS - sqrt(distanceSqr);
                            float heightCurrent = m_heightmap->Get(i, j);
                            float newHeight = heightCurrent - deltaSediment * weightCurrent / weightSum;
                            if (newHeight < m_heightmap->m_minHeight)
                            {
                                newHeight = m_heightmap->m_minHeight;
                            } 
                            // heightmapHolder->Set(i, j, newHeight);s
                            m_heightmap->Set(i, j, newHeight);
                        }
                    }
                }
                // Erosion::DebugDropSetColor(glm::vec3(0.9,0.0,0.0));
            }
            else
            {
                float deltaSediment = (drop.sediment - capacity) * DEPOSITION;
                    // (heightDifference < 0) ? 
                    // (drop.sediment - capacity) * DEPOSITION :
                    // std::min(drop.sediment, heightDifference);
                // debug("Dropped", deltaSediment);
                drop.sediment -= deltaSediment;
                surroundingRect = HeightMap::GetSurroundingRect(positionOldTC);
                if (!m_heightmap->IsPointValid(surroundingRect[0][0], surroundingRect[0][1]) ||
                !m_heightmap->IsPointValid(surroundingRect[1][0], surroundingRect[1][1]))
                {
                    break;
                    // return;
                }

                for (int i = 0; i < 4; i++)
                {
                    int coords[] = {i / 2, i % 2};
                    float currentHeight = m_heightmap->Get(
                        surroundingRect[coords[0]][0], 
                        surroundingRect[coords[1]][1]
                    );

                    float sedimentPart = (1 - std::abs(surroundingRect[coords[0]][0] - positionOldTC.x)) * 
                                            (1 - std::abs(surroundingRect[coords[1]][1] - positionOldTC.y));


                    // heightmapHolder->Set(
                    m_heightmap->Set(
                        surroundingRect[coords[0]][0],
                        surroundingRect[coords[1]][1],
                        currentHeight + sedimentPart * deltaSediment
                    );
                }
                // Erosion::DebugDropSetColor(glm::vec3(0.0,0.9,0.0));
            }
        }
    }

    

    debug("Eroding in progress:100%");
    debug("Erosion complete!");

    debug("Blurring Erosion");
    // *m_heightmap -= *m_savedHeightmap; 
    // // delete m_savedHeightmap;
    // // m_savedHeightmap = new HeightMap(*m_heightmap);
    m_heightmap->BlurHeightmap();
    // // m_savedHeightmap->BlurHeightmap();
    // *m_heightmap += *m_savedHeightmap;
    debug("Blurring complete!");

    RecreateMesh();    
    m_heightmap->SaveAsTexture();
}

void Landscape::Update()
{
    // debug(Component::AnyComponentTyped<CameraComponent>()->GetPosition());
    // glm::vec2 pos = PointRealToMapNormalized(m_heightmap, XZ(Component::AnyComponentTyped<CameraComponent>()->GetPosition()));
    // BiomePoint bp = m_biomemap->Get(m_biomemap->PointNormalizedToInt(pos));
    // debug(bp.intensity, (int)bp.biome);
    // debug(PointRealToMapNormalized(m_heightmap, XZ(Component::AnyComponentTyped<CameraComponent>()->GetPosition())));
}

void Landscape::ApplyShoreImmerse()
{
    if (m_surfacemap == nullptr || m_heightmap == nullptr )
    {
        debug("MAP IS NOT AVAILABLE! [ApplyShoreImmerse]");
        throw std::exception{};
    }


    const float SHORE_BIOMEMAP_CELLS_SIZE = 
        WaterGenerator::SHORE_REAL_SIZE / GetMapCellRealSize(m_surfacemap);
    const float SHALLOW_WATER_BIOMEMAP_CELLS_SIZE = 
        WaterGenerator::SHALLOW_WATER_REAL_SIZE / GetMapCellRealSize(m_surfacemap);
    std::vector<glm::vec2> shorePoints = {{0,0}, {0.2, 0.1}, {0.8, 0.9}, {1,1}};
    glm::vec2 waterShallowCoefs = {0.9, 0.1}; // a * x^2 + b * x
    std::array<glm::vec2, 4> waterShallowBezierPoints = 
    {
        glm::vec2{0.f, 0.f}, 
        glm::vec2{.5f, 0.f},
        glm::vec2{.75f, 1.f},
        glm::vec2{1.f, 1.f}  
    };


    waterLevel = CalculateWaterLevel(); //m_heightmap->m_minHeight;
    CreateWaterPlate();


    //SHARED FUNCTION
    const int APPROACH_SEGMENTS = 500;
    std::function<float(float, std::array<glm::vec2, APPROACH_SEGMENTS>&)> CalculateByPoints = 
        [](float x, std::array<glm::vec2, APPROACH_SEGMENTS>& points)
    {
        int start = 0;
        int end = APPROACH_SEGMENTS - 1;
        while (end - start > 1)
        {
            int middle = (start + end) / 2;
            if (x >= points[middle].x)
                start = middle;
            else
                end = middle;
        }
        float ratio = 0;
        if (x - start > 0.000001)
            ratio = (float)(x - start) / (end - start);
        return LERP(points[start].y, points[end].y, ratio);
    };

    //FOR SHORELINE CALCULATION
    std::array<float, 4> coefs{-1.19209e-07, 2.38889, -4.16668, 2.77778};
/*
    Eigen::MatrixXf A(4, 4);
    Eigen::VectorXf b(4);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (j == 0)
                A(i, j) = 1;
            else
                A(i, j) = pow(shorePoints[i].y, j);
        }
        b(i) = shorePoints[i].x;
    }
    Eigen::VectorXf answer = A.householderQr().solve(b);
    for (int i = 0; i < 4; i++)
    {
        coefs[i] = answer(i);
    }
*/
    std::function<float(float)> polynom = [&coefs](float y)
    {
        float res = 0;
        for (int i = 0; i < 4; i++)
        {
            res += coefs[i] * pow(y, i);
        }
        return res;
    };
    
    // const int APPROACH_SEGMENTS = 500;
    std::array<glm::vec2, APPROACH_SEGMENTS> approachPoints;
    for (int i = 0; i < APPROACH_SEGMENTS; i++)
    {
        approachPoints[i].y = (float) i / APPROACH_SEGMENTS;
        approachPoints[i].x = polynom(approachPoints[i].y);
    }
    
    std::function<float(float)> shoreY = [CalculateByPoints, &approachPoints](float x)
    {
        return CalculateByPoints(x, approachPoints);
    };


    //FOR SHALLOW WATER CALCULATION
    std::function<glm::vec2(float)> waterBezier = [&waterShallowBezierPoints](float t)
    {
        glm::vec2 result;
        glm::vec2 p0 = waterShallowBezierPoints[0];
        glm::vec2 p1 = waterShallowBezierPoints[1];
        glm::vec2 p2 = waterShallowBezierPoints[2];
        glm::vec2 p3 = waterShallowBezierPoints[3];

        result.x = (1-t) * ((1-t) * ((1-t) * p0.x + t * p1.x) + 
            t * ((1-t) * p1.x + t * p2.x)) + 
            t * ((1-t) * ((1-t) * p1.x + t * p2.x) + 
            t * ((1-t) * p2.x + t * p3.x));
        result.y = (1-t) * ((1-t) * ((1-t) * p0.y + t * p1.y) + 
            t * ((1-t) * p1.y + t * p2.y)) + 
            t * ((1-t) * ((1-t) * p1.y + t * p2.y) + 
            t * ((1-t) * p2.y + t * p3.y));

        return result;
    };

    std::array<glm::vec2, APPROACH_SEGMENTS> waterApproachPoints;
    for (int i = 0; i < APPROACH_SEGMENTS; i++)
    {
        float t = (float) i / APPROACH_SEGMENTS;
        waterApproachPoints[i] = waterBezier(t);
    }

    std::function<float(float)> waterY = [CalculateByPoints, &waterApproachPoints](float x)
    {
        return CalculateByPoints(x, waterApproachPoints);
    };

    for (int i = -m_heightmap->m_width; i <= m_heightmap->m_width; i++)
    {
        for (int j = -m_heightmap->m_height; j <= m_heightmap->m_height; j++)
        {
            vec2Int coords = {i, j};
            glm::vec2 vecCoords = m_heightmap->PointIntToNormalized(coords);

            std::vector<SurfacePoint> surroundingPoints = m_surfacemap->GetSurfaceInterpolated_Normalized(vecCoords);
            SurfacePoint surfacePoint = surroundingPoints[0];
            for (SurfacePoint sp : surroundingPoints)
            {
                if (sp.intensity > surfacePoint.intensity)
                    surfacePoint = sp;
            }
                
            if (surfacePoint.surface == Surface::LIQUID)
            {
                noise::module::Perlin perlin;
                perlin.SetFrequency(WaterGenerator::FREQUENCY);
                perlin.SetOctaveCount(WaterGenerator::OCTAVE_COUNT);
                perlin.SetPersistence(WaterGenerator::PERSISTANCE);
                perlin.SetLacunarity(WaterGenerator::LACUNARITY);

                glm::vec2 curPos = PointMapIntToReal(m_heightmap, vec2Int(i, j));
                float height = (float) perlin.GetValue(curPos[0], 0, curPos[1]);
                height /= 1.2;
                height = clamp(height, -1.0f, 1.0f);
                height = height * 0.5 + 0.5;
                height *= WaterGenerator::HEIGHT_DELTA;
                height = waterLevel - WaterGenerator::SEA_DEPTH - height;

                float realIntensity = surfacePoint.intensity;
                for (SurfacePoint sp : surroundingPoints)
                {
                    if (sp.surface == Surface::EARTH)
                        realIntensity -= sp.intensity;
                }
                float x = realIntensity / SHALLOW_WATER_BIOMEMAP_CELLS_SIZE;
                float ratio;
                if (x < 0 || x > 1)
                    ratio = clamp(x, 0.0f, 1.0f);
                else
                    ratio = waterY(x);

                m_heightmap->Set(coords, LERP(waterLevel, height, ratio));
            }
            else if (surfacePoint.surface == Surface::EARTH)
            {
                float realIntensity = surfacePoint.intensity;
                for (SurfacePoint sp : surroundingPoints)
                {
                    if (sp.surface == Surface::LIQUID)
                        realIntensity -= sp.intensity;
                }
                
                float x = realIntensity / SHORE_BIOMEMAP_CELLS_SIZE;
                float ratio;
                if (x < 0 || x > 1)
                    ratio = clamp(x, 0.0f, 1.0f);
                else
                    ratio = shoreY(x);

                float currHeight = m_heightmap->Get(coords);
                m_heightmap->Set(coords, LERP(waterLevel, currHeight, ratio));
            }
        }
    }

    m_heightmap->SaveAsTexture();
}

void Landscape::CreateWaterPlate()
{
    glm::vec3 pos;
    pos.y = waterLevel;
    glm::vec2 p = 
        PointMapIntToReal(m_surfacemap, m_surfacemap->PointNormalizedToInt(glm::vec2{}));
    pos.x = p.x;
    pos.z = p.y;
    glm::vec2 size = {m_realSize[0], m_realSize[1]};
    waterPlate = std::make_unique<WaterPlate>(pos, size);
}

WaterPlate* Landscape::GetWaterPlate()
{
    if (!waterPlate)
    {
        debug("NO WATER PLATE!");
        throw std::exception{};
    }
    return waterPlate.get();
}

bool Landscape::IsRenderingWater()
{
    return static_cast<bool>(waterPlate);
}

float Landscape::GetWaterLevel()
{
    return waterLevel;
}

float Landscape::GetWaterLevelSafeDelta()
{
    return 0.04;
}

bool Landscape::IsPointErodable(glm::vec2 coords)
{
    BiomePoint bp = m_biomemap->Get(m_biomemap->PointNormalizedToInt(coords));
    if (IsBiomeErodable(bp.biome))
        return true;
    else
    {
        SurfacePoint sp = m_surfacemap->GetSurface(m_surfacemap->PointNormalizedToInt(coords));
        if (sp.surface == Surface::EARTH)
        {
            const float SHORE_SURFACEMAP_CELLS_SIZE = 
                WaterGenerator::SHORE_REAL_SIZE / GetMapCellRealSize(m_surfacemap);
            if (sp.intensity / SHORE_SURFACEMAP_CELLS_SIZE < 1.0)
                return true;
        }
    }
    return false;
}

bool Landscape::IsBiomeErodable(Biome b)
{
    // return b <= Biome::GROUND;
    return b == Biome::MOUNTAINS;
}

bool Landscape::IsPointSedimentDropable(glm::vec2 heightmapTC)
{
    vec2Int heightmapTCInt = vec2Int{int(heightmapTC.x), int(heightmapTC.y)};
    glm::vec2 heightmapTCNormalized = m_heightmap->PointIntToNormalized(heightmapTCInt);
    BiomePoint bp = m_biomemap->Get(m_biomemap->PointNormalizedToInt(heightmapTCNormalized));
    if (!(bp.biome <= Biome::GROUND))
        return false;
    return bp.intensity >= WaterGenerator::SHORE_REAL_SIZE / GetMapCellRealSize(m_biomemap);
}

void Landscape::GenerateBiomemap()
{
    vec2Int biomeMapSize = vec2Int(glm::vec2{m_realSize[0], m_realSize[1]} * CGenBiomeMap::DENSITY_PER_REAL);
    m_biomemap = new CGenBiomeMap(biomeMapSize[0], biomeMapSize[1]);

    float m_width = m_biomemap->m_width;
    float m_height = m_biomemap->m_height;

    debug("Generating shoreline points");

    const float COAST_STRAIGHT_SEA_SIZE = 
        Landscape::WaterBiomeGenerator::COAST_STRAIGHT_SEA_SIZE;
    const float COAST_STRAIGHT_MIDPOINT_MAX_ANGLE_DEGREES = 
        Landscape::WaterBiomeGenerator::COAST_STRAIGHT_MIDPOINT_MAX_ANGLE_DEGREES;
    const int COAST_STRAIGHT_MIDPOINT_DIVISION_TIMES = 
        Landscape::WaterBiomeGenerator::COAST_STRAIGHT_MIDPOINT_DIVISION_TIMES;

    class Shoreline
    {
        public:
            std::vector<glm::vec2> points;
            std::vector<vec2Int> coords;
    };

    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            m_biomemap->Set(i, j, Biome::EMPTY);
        }
    }  

    Shoreline shoreline;
    shoreline.points.push_back(glm::vec2 {0.0, 1 - COAST_STRAIGHT_SEA_SIZE});
    shoreline.points.push_back(glm::vec2 {1.0, 1 - COAST_STRAIGHT_SEA_SIZE});

    for (int division = 0; division < COAST_STRAIGHT_MIDPOINT_DIVISION_TIMES; division++)
    {
        const int CHANGEABLE_DIMENSION = 1;
        auto iter = shoreline.points.begin();
        while ((iter + 1) != shoreline.points.end())
        {
            glm::vec2 start = *iter;
            glm::vec2 end = *(iter + 1);

            float maxAngle = COAST_STRAIGHT_MIDPOINT_MAX_ANGLE_DEGREES;
            float deviationAngle = (((float) (rand()) / RAND_MAX) * 2 - 1) * maxAngle;

            glm::vec2 forward = glm::normalize(end - start);
            glm::vec2 leftCross = glm::vec2(-forward.y, forward.x);
            glm::vec2 toMidpoint = cosf(glm::radians(deviationAngle)) * forward +
                                   sinf(glm::radians(deviationAngle)) * leftCross;

            float toMidpointLen = glm::length(end - start) / 2 / cosf(glm::radians(deviationAngle));
            glm::vec2 midpoint = toMidpoint * toMidpointLen + start;
            
            //Check if new point is not too close to border
            const float DELTA = std::max(
                COAST_STRAIGHT_SEA_SIZE * 0.2f,
                0.0000001f + 1.0f / (2 * std::min(m_width, m_height) + 1)
            );

            if (midpoint[CHANGEABLE_DIMENSION] < DELTA ||
                midpoint[CHANGEABLE_DIMENSION] > 1.0 - DELTA)
            {
                continue;
            }
            else
            {
                iter = shoreline.points.insert(iter + 1, midpoint);
                iter++;
            }
        }   
    }

    debug("Done!");
    debug("Drawning shoreline");

    float stepLength = 1.0f / (std::max(m_width, m_height) * 3 + 1); 
    for (int i = 0; i < shoreline.points.size() - 1; i++)
    {
        glm::vec2 start = shoreline.points[i];
        glm::vec2 end = shoreline.points[i+1];
        glm::vec2 curPoint = start;
        while (glm::length(curPoint - start) - stepLength < glm::length(end - start))
        {
            glm::vec2 realPoint = curPoint;
            if (glm::length(curPoint - start) > glm::length(end - start))
                realPoint = end;
            vec2Int curCoords = m_biomemap->PointNormalizedToInt(realPoint);
            if (m_biomemap->Get(curCoords) == Biome::EMPTY)
            {
                m_biomemap->Set(curCoords, Biome::SHORELINE);
                shoreline.coords.push_back(curCoords);
            }

            curPoint += stepLength * glm::normalize(end - start);
        }
    }

    debug("Done!");
    debug("Filling ground and water");

    //Setting starting points for sea and ground
    vec2Int seaBreeder = {0, m_height};
    m_biomemap->Set(seaBreeder,Biome::WATER);
    vec2Int groundBreeder = {0, -m_height};
    m_biomemap->Set(groundBreeder, Biome::GROUND);

    //Creating functions for breadth-first search

    //FIND NEIGHBOURS FUNCTIONS
        std::function biomeFindNeighbs_FindEmpty
        {
            [this]
            (vec2Int coords) 
            {
                return m_biomemap->biomeFindNeighbs_FindBiome(coords, Biome::EMPTY);
            }
        };

        std::function biomeFindNeighbs_FindWater
        {
            [this]
            (vec2Int coords) 
            {
                return m_biomemap->biomeFindNeighbs_FindBiome(coords, Biome::WATER);
            }
        };

        std::function biomeFindNeighbs_FindGround
        {
            [this]
            (vec2Int coords) 
            {
                return m_biomemap->biomeFindNeighbs_FindBiome(coords, Biome::GROUND);
            }
        };

    //APPLY FUNCTIONS
        std::function biomeApplyToReachedPoints_SetWater
        {
            [this](vec2Int coords, float distance)
            {
                m_biomemap->Set(coords, Biome::WATER);
            }
        };

        //SHORELINE EXCLUDED
        std::function biomeApplyToReachedPoints_SetWater_DistanceWise
        {
            [this](vec2Int coords, float distance)
            {
                if (m_biomemap->Get(coords) == Biome::WATER)
                {
                    m_biomemap->Set(coords, BiomePoint{Biome::WATER, distance});
                }
            }
        };

        //SHORELINE INCLUDED
        std::function biomeApplyToReachedPoints_SetGround_DistanceWise
        {
            [this](vec2Int coords, float distance)
            {
                Biome b = m_biomemap->Get(coords);
                if (b == Biome::GROUND || b == Biome::SHORELINE)
                {
                    m_biomemap->Set(coords, BiomePoint{Biome::GROUND, distance});
                }
            }
        };

        std::function biomeApplyToReachedPoints_SetGround
        {
            [this](vec2Int coords, float distance)
            {
                m_biomemap->Set(coords, Biome::GROUND);
            }
        };

    //Fill water
    biomeBreadthSearch(
        std::vector<vec2Int> {seaBreeder},
        biomeFindNeighbs_FindEmpty,
        biomeApplyToReachedPoints_SetWater,
        m_biomemap->biomeFindDistance_FromOrigin
    );
    //Fill ground
    biomeBreadthSearch(
        std::vector<vec2Int> {groundBreeder},
        biomeFindNeighbs_FindEmpty,
        biomeApplyToReachedPoints_SetGround,
        m_biomemap->biomeFindDistance_FromOrigin
    );
    //Fill shoreline with ground
    auto it = shoreline.coords.begin();
    while (it != shoreline.coords.end())
    {
        m_biomemap->Set(*it, Biome::GROUND);
        it++;
    }

    debug("Done!");


    debug("Generating mountain and plain Voronoi areas");

    // const float DENSITY = MountainBiomeGenerator::POINTS_DENSITY;
    const float AMOUNT = MountainBiomeGenerator::POINTS_AMOUNT;
    const float PLAIN_RADIUS = MountainBiomeGenerator::CENTER_PLAIN_RADIUS;
    const float MIN_CHANCE = MountainBiomeGenerator::MOUNTAIN_MIN_PROBABILITY;
    const float MAX_CHANCE = MountainBiomeGenerator::MOUNTAIN_MAX_PROBABILITY;

    struct VoronoiPoint
    {
        BiomePoint biome;
        glm::vec2 coordsNormalized;
    };
    
    std::vector<VoronoiPoint> points;
    // int pointsAmount = int(m_realSize[0] * m_realSize[1] * MountainBiomeGenerator::POINTS_DENSITY);
    int pointsAmount = MountainBiomeGenerator::POINTS_AMOUNT;
    {
        int i = 0;
        while (i < pointsAmount)
        {
            glm::vec2 coords = {rndNormalized(), rndNormalized()};

            if (m_biomemap->Get(m_biomemap->PointNormalizedToInt(coords)).biome != Biome::GROUND)
            {
                continue;
            }

            float distToCenter = glm::length(coords - glm::vec2(0.5));
            float mountainProbability = 0.f;
            if (distToCenter > PLAIN_RADIUS)
            {
                mountainProbability = distToCenter - PLAIN_RADIUS;
                mountainProbability = mountainProbability / (sqrt(2) * 0.5 - PLAIN_RADIUS);
                mountainProbability = LERP(MIN_CHANCE, MAX_CHANCE, mountainProbability);
            }  
            points.emplace_back(); 
            points[i].biome = (rndNormalized() < mountainProbability) ? Biome::MOUNTAINS : Biome::PLAIN;
            points[i].coordsNormalized = coords;

            i++;
        }
    }

    noise::module::Perlin mountainPerlin;
    mountainPerlin.SetFrequency(MountainBiomeGenerator::PERLIN_FREQUENCY_REAL);
    mountainPerlin.SetOctaveCount(MountainBiomeGenerator::PERLIN_OCTAVE_COUNT);
    mountainPerlin.SetPersistence(MountainBiomeGenerator::PERLIN_PERSISTANCE);
    mountainPerlin.SetLacunarity(MountainBiomeGenerator::PERLIN_LACUNARITY);
    float BIOMEMAP_RATIO = MountainBiomeGenerator::PERLIN_AFFECT_RATIO;

    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            vec2Int curCoords{i,j};
            BiomePoint curPoint = m_biomemap->Get(curCoords);
            
            if (curPoint.biome != Biome::GROUND)
                continue;

            glm::vec2 curCoordsReal = PointMapIntToReal(m_biomemap, curCoords);
            float perlin = mountainPerlin.GetValue(curCoordsReal.x, 0, curCoordsReal.y);
            perlin /= 1.45;
            perlin = perlin * 0.5f + 0.5f;
            perlin = clamp(perlin, 0.f, 1.f);
            perlin *= BIOMEMAP_RATIO; 
            float modifier;

            float minDistance = glm::length(2.f * glm::vec2(m_width, m_height));
            Biome biome;
            for (int p = 0; p < points.size(); p++)
            {
                float len = 
                    (curCoords - m_biomemap->PointNormalizedToInt(points[p].coordsNormalized)).length();

                modifier = 0.f;
                if (points[p].biome == Biome::PLAIN)
                    modifier = 1 + perlin;
                else if (points[p].biome == Biome::MOUNTAINS)
                    modifier = 1 - perlin;
                else
                    modifier = 1.f;
                len *= modifier;

                if (minDistance > len)
                {
                    minDistance = len;
                    biome = points[p].biome;
                }
            }

            m_biomemap->Set(curCoords, BiomePoint{biome, curPoint.intensity});
        }
    }

    //Calculating mountains biome real intensity
    std::vector<vec2Int> startPoints;
    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            vec2Int curCoords{i,j};
            BiomePoint bp = m_biomemap->Get(curCoords);
            if (bp.biome == Biome::PLAIN || bp.biome == Biome::WATER)
                startPoints.push_back(curCoords);
        }
    }

    std::function biomeFindNeighbs_FindMountains
    {
        [this]
        (vec2Int coords) 
        {
            return m_biomemap->biomeFindNeighbs_FindBiome(coords, Biome::MOUNTAINS);
        }
    };

    std::function biomeApplyToReachedPoints_SetMountains_DistanceWise
    {
        [this](vec2Int coords, float distance)
        {
            Biome b = m_biomemap->Get(coords);
            if (b == Biome::MOUNTAINS)
            {
                m_biomemap->Set(coords, BiomePoint{Biome::MOUNTAINS, distance});
            }
        }
    };

    biomeBreadthSearch(
        startPoints,
        biomeFindNeighbs_FindMountains,
        biomeApplyToReachedPoints_SetMountains_DistanceWise,
        m_biomemap->biomeFindDistance_FromOrigin
    );

    debug("Done!");


    debug("Generating water map");
    m_surfacemap = SurfaceMap::ScanBiomemap(*m_biomemap);
    debug("Done");

    debug("Calibrating biome map");
    m_biomemap->CalibrateBiomemap();
    debug("Done");

    // m_biomemap->CalibrateBiomemap();
    // BiomeMap test(*m_biomemap);
    // test.CalibrateBiomemap();
    // for (int i = -test.m_width; i <= test.m_width; i++)
    // {
    //     for (int j = -test.m_height; j <= test.m_height; j++)
    //     {
    //         vec2Int coords{i, j};
    //         BiomePoint bp1 = test.Get(coords);
    //         BiomePoint bp2 = m_biomemap->Get(coords);
    //         // m_biomemap->Set(coords, BiomePoint{Biome::PLAIN, 2.f});
    //         if (bp1.biome != bp2.biome)
    //         {
    //             debug("FFFF");
    //             throw std::exception{};
    //         }
    //         if (fabs(bp1.intensity - bp2.intensity) > 0.001)
    //         {
    //             if (coords.length() < 20)
    //                 continue;

    //             if (fabs(bp1.intensity - bp2.intensity) < 5.1f)
    //                 m_biomemap->Set(coords, BiomePoint{Biome::MOUNTAINS, 2.f});
    //             else
    //                 m_biomemap->Set(coords, BiomePoint{Biome::WATER, 2.f});
    //         }
    //     }
    // }
    // m_biomemap->SaveAsTexture();
}

MAP& operator++(MAP& map)
{
    return map = static_cast<MAP>(static_cast<int>(map) + 1);
}

void Landscape::GenerateRoad()
{
    debug("Generating road scheme!");
    m_roadScheme = new RoadScheme(m_heightmap->m_width / m_heightmap->m_height);

    // STARTING POINTS
    float len = RoadSchemeGenerator::MIN_GENERATION_SEGMENT_LENGTH; //real length
    glm::vec2 posFrom = PointMapIntToReal(m_populationmap, m_populationmap->shoreCenterPos);
    glm::vec2 posTo = PointMapIntToReal(m_populationmap, m_populationmap->realCenterPos);
    posTo = glm::normalize(posTo - posFrom) * len + posFrom;

    posFrom = PointRealToMapNormalized(m_heightmap, posFrom);
    posTo = PointRealToMapNormalized(m_heightmap, posTo);

    RoadNode* from = m_roadScheme->graph.AddNode(posFrom);
    RoadNode* to = m_roadScheme->graph.AddNode(posTo);
    m_roadScheme->graph.AddSection(from, to, 0.f);  


    // //Testing angles
    // float testAngleDegrees = 90.f;
    // {
    //     float angle = glm::radians(testAngleDegrees);
    //     float _len = glm::length(from->pos - to->pos);
    //     glm:;vec2 forward = glm::normalize(posFrom - posTo);
    //     glm::vec2 left = LeftNormal(forward);
    //     glm::vec2 newPos = posTo + (cosf(angle) * forward + sinf(angle) * left) * _len;
    //     RoadNode* newNode = m_roadScheme->graph.AddNode(newPos);
    //     m_roadScheme->graph.AddSection(newNode, to, 0.f);  
    // }
    
    
    // float l = 7.f;
    //     float lNorm = glm::length(PointRealToMapNormalized(m_biomemap, glm::vec2{0,l}) - 
    //     PointRealToMapNormalized(m_biomemap, glm::vec2{0,0}));
    // glm::vec2 centerPos = m_populationmap->PointIntToNormalized(m_populationmap->realCenterPos);
    // auto center = m_roadScheme->graph.AddNode(centerPos);
    // auto bot = m_roadScheme->graph.AddNode(centerPos - glm::vec2{0, lNorm});
    // auto nl = m_roadScheme->graph.AddNode(centerPos - glm::vec2{lNorm / 2, lNorm * 2});
    // auto nr = m_roadScheme->graph.AddNode(centerPos - glm::vec2{-lNorm / 2, lNorm * 2});
    // m_roadScheme->graph.AddSection(center, bot, 0.f);
    // m_roadScheme->graph.AddSection(bot, nl, 0.f);
    // m_roadScheme->graph.AddSection(bot, nr, 0.f);
    // debug("POINTS");
    // debug(center->pos);
    // debug(bot->pos);
    // debug(nl->pos);
    // debug(nr->pos);

    // //SOME crossroads
    // const int counter = 10;
    // std::array<std::array<RoadNode*, (size_t)counter>, (size_t)counter> nodes;
    // for (int i = 0; i < counter; i++)
    // {
    //     for (int j = 0; j < counter; j++)
    //     {
    //         glm::vec2 pos{(1.0f / (counter + 1)) * (i + 1), (1.0f / (counter + 1)) * (j + 1)};
    //         nodes[i][j] = m_roadScheme->graph.AddNode(pos);
    //     }
    // }
    // for (int i = 0; i < counter; i++)
    // {
    //     for (int j = 0; j < counter; j++)
    //     {
    //         std::vector<vec2Int> friends;
    //         if (i + 1 < counter)
    //             friends.push_back(vec2Int{i+1, j});
    //         if (j + 1 < counter)
    //             friends.push_back(vec2Int{i, j+1});
    //         for (vec2Int p: friends)
    //             m_roadScheme->graph.AddSection(nodes[i][j], nodes[p.x][p.y], 0.f);
    //     }
    // }
    

    float ratio = 0.4;
    // //Y
    // auto root = m_roadScheme->graph.AddNode(glm::vec2{0.5,0.2} * ratio);
    // auto n2 = m_roadScheme->graph.AddNode(glm::vec2{0.5,0.5} * ratio);
    // auto n3 = m_roadScheme->graph.AddNode(glm::vec2{0.77,0.66} * ratio);
    // auto n4 = m_roadScheme->graph.AddNode(glm::vec2{0.3,0.7} * ratio);
    // m_roadScheme->graph.AddSection(root, n2, 0.f);
    // m_roadScheme->graph.AddSection(n3, n2, 0.f);
    // m_roadScheme->graph.AddSection(n2, n4, 0.f);
    // // // CROSS
    // auto center = m_roadScheme->graph.AddNode(glm::vec2{0.5,0.1} * ratio);
    // // auto c1 = m_roadScheme->graph.AddNode(glm::vec2{0.5,0.01} * ratio);
    // auto c2 = m_roadScheme->graph.AddNode(glm::vec2{0.4,0.1} * ratio);
    // // auto c3 = m_roadScheme->graph.AddNode(glm::vec2{0.6,0.1} * ratio);
    // m_roadScheme->graph.AddSection(root, center, 0.f);
    // // m_roadScheme->graph.AddSection(c1, center, 0.f);
    // m_roadScheme->graph.AddSection(c2, center, 0.f);
    // // m_roadScheme->graph.AddSection(c3, center, 0.f);
    // //э cross
    // auto n5 = m_roadScheme->graph.AddNode(glm::vec2{0.9,0.9} * ratio);
    // // auto n6 = m_roadScheme->graph.AddNode(glm::vec2{0.7,0.9} * ratio);
    // // auto n7 = m_roadScheme->graph.AddNode(glm::vec2{0.5,0.9} * ratio);
    // // auto n8 = m_roadScheme->graph.AddNode(glm::vec2{0.9,0.5} * ratio);
    // // auto n9 = m_roadScheme->graph.AddNode(glm::vec2{0.7,0.5} * ratio);
    // m_roadScheme->graph.AddSection(n3, n5, 0.f);
    // // m_roadScheme->graph.AddSection(n3, n6, 0.f);
    // // m_roadScheme->graph.AddSection(n3, n7, 0.f);
    // // m_roadScheme->graph.AddSection(n3, n8, 0.f);
    // // m_roadScheme->graph.AddSection(n3, n9, 0.f);

    // //DEBUG TEST
    // auto center = m_roadScheme->graph.AddNode(glm::vec2{0.427512, 0.489259});
    // auto c1 = m_roadScheme->graph.AddNode(glm::vec2{0.51304, 0.468397});
    // auto c2 = m_roadScheme->graph.AddNode(glm::vec2{0.360223, 0.42509});
    // auto c3 = m_roadScheme->graph.AddNode(glm::vec2{.376695, 0.501654});
    // m_roadScheme->graph.AddSection(center, c1, 0.f);
    // m_roadScheme->graph.AddSection(center, c2, 0.f);
    // m_roadScheme->graph.AddSection(center, c3, 0.f);

    // // circle
    // auto r1 = m_roadScheme->graph.AddNode(glm::vec2{0.3,0.3} * ratio);
    // auto r2 = m_roadScheme->graph.AddNode(glm::vec2{0.1,0.5} * ratio);
    // auto r3 = m_roadScheme->graph.AddNode(glm::vec2{0.3,0.7} * ratio);
    // auto r4 = m_roadScheme->graph.AddNode(glm::vec2{0.5,0.5} * ratio);
    // m_roadScheme->graph.AddSection(r1, r2, 45.f);
    // m_roadScheme->graph.AddSection(r2, r3, 45.f);
    // m_roadScheme->graph.AddSection(r3, r4, 45.f);
    // m_roadScheme->graph.AddSection(r4, r1, 45.f);

    // // square
    // auto r1 = m_roadScheme->graph.AddNode(glm::vec2{0.3,0.3} * ratio);
    // auto r2 = m_roadScheme->graph.AddNode(glm::vec2{0.4,0.5} * ratio);
    // auto r3 = m_roadScheme->graph.AddNode(glm::vec2{0.6,0.5} * ratio);
    // auto r4 = m_roadScheme->graph.AddNode(glm::vec2{0.6,0.3} * ratio);
    // auto r5 = m_roadScheme->graph.AddNode(glm::vec2{0.7,0.7} * ratio);
    // m_roadScheme->graph.AddSection(r1, r2, 0.f);
    // m_roadScheme->graph.AddSection(r2, r3, 0.f);
    // m_roadScheme->graph.AddSection(r3, r4, 0.f);
    // m_roadScheme->graph.AddSection(r4, r1, 0.f);
    // m_roadScheme->graph.AddSection(r3, r5, 0.f);


    //Road scheme generation logic

    //For all iterations
    PopulationMap populationDistances(*m_populationmap);
    std::vector<RoadNode*> roadNodes;
    std::vector<RoadSection*> roadSections;

    std::function UpdateGraph
    {
        [&populationDistances, &roadNodes, &roadSections, this]()
        {
            roadNodes = m_roadScheme->graph.GetNodes();
            roadSections = m_roadScheme->graph.GetSections();
            std::sort(
                roadNodes.begin(),
                roadNodes.end(),
                [] (RoadNode* a, RoadNode* b) 
                { 
                    return (a->pos.x + 10000*a->pos.y) > (b->pos.x + 10000*b->pos.y);
                }
            );
            std::sort(
                roadSections.begin(),
                roadSections.end(),
                [] (RoadSection* a, RoadSection* b) 
                { 
                    return (a->start->pos.x + 10*a->start->pos.y + 100*a->end->pos.x + 1000*a->end->pos.y ) > 
                        (b->start->pos.x + 10*b->start->pos.y + 100*b->end->pos.x + 1000*b->end->pos.y );
                }
            );
            for (RoadNode* nodePtr: roadNodes)
            {
                nodePtr->data = 0.f;
            }

            for (int i = -populationDistances.m_width; i <= populationDistances.m_width; i++)
            {
                for (int j = -populationDistances.m_height; j <= populationDistances.m_height; j++)
                { 
                    vec2Int coords{i, j};
                    PopulationPoint pp = m_populationmap->Get(coords);
                    if (pp.population != Population::POPULATED)
                        continue;
                    
                    glm::vec2 coordsNorm = m_populationmap->PointIntToNormalized(coords);
                    RoadSection* nearestSecPtr = nullptr;
                    float minDist = 1000;
                    for (RoadSection* sectionPtr: roadSections)
                    { 
                        RoadSection& section = *sectionPtr;
                        float curDist = DistanceToSection(coordsNorm, section.start->pos, section.end->pos);

                        if (minDist > curDist)
                        {
                            nearestSecPtr = sectionPtr;
                            minDist = curDist;
                        }
                    }

                    glm::vec2 sectionPoint = NearestPointOfSection(coordsNorm, nearestSecPtr->start->pos, nearestSecPtr->end->pos);
                    float len = glm::length(nearestSecPtr->end->pos - nearestSecPtr->start->pos);
                    float population = m_populationmap->Get(coords).intensity;
                    float populationWeighed = population * minDist;

                    nearestSecPtr->start->data += glm::length(nearestSecPtr->end->pos - sectionPoint) / len * populationWeighed;
                    nearestSecPtr->end->data += glm::length(nearestSecPtr->start->pos - sectionPoint) / len * populationWeighed;
                    populationDistances.Set(coords, PopulationPoint(Population::POPULATED, minDist));
                }   
            }
        }
    };

    //For each iteration
    float final = RoadSchemeGenerator::SECTION_FINAL_MAX_POPULATION_RATIO;
    final *= PopulationGenerator::POPULATION_VOLUME_REAL;

    enum class EXPANSION_TYPE
    {
        GROWTH,
        LINKING
    } expansionType;
    expansionType = EXPANSION_TYPE::GROWTH;

    enum class LINKING_PHASE
    {
        NONE,
        MERGE_CROSSROAD,
        MERGE_MIDROAD,
        CHECKING_FRINGE,
        EXTEND_DEADEND

    } linkingPhase;
    linkingPhase = LINKING_PHASE::NONE;

    std::vector<RoadNode*> permitedEndNodes;
    std::vector<RoadNode*> deniedDeadendNodes;
    std::vector<RoadNode*> deniedExtensions;

    while(true)
    {
        // debug("1");
        UpdateGraph();
        RoadSection* breedingSection = nullptr;
        RoadNode* breedingNode = nullptr;
        RoadNode* breedingEndNode = nullptr;
        glm::vec2 breedingEndPos;

        if (expansionType == EXPANSION_TYPE::GROWTH)
        {
            // debug("2.1");
            float maxPopulation = 0.f;
            for (RoadSection* sectionPtr: roadSections)
            {
                float curPopulation = sectionPtr->start->data + sectionPtr->end->data;
                if (curPopulation > maxPopulation)
                {
                    if (sectionPtr->blocked && sectionPtr->start->blocked && sectionPtr->end->blocked)
                        continue;

                    // debug("2_1_1", sectionPtr->blocked, sectionPtr->start->blocked, sectionPtr->end->blocked);
                    maxPopulation = curPopulation;
                    breedingSection = sectionPtr;
                }
            }
            // debug("2.2");
            float realSectionPopulation = 0.f;
            for (int i = -populationDistances.m_width; i <= populationDistances.m_width; i++)
            {
                for (int j = -populationDistances.m_height; j <= populationDistances.m_height; j++)
                { 
                    vec2Int coords{i, j};
                    PopulationPoint pp = m_populationmap->Get(coords);
                    if (pp.population != Population::POPULATED)
                        continue;
                    glm::vec2 coordsNorm = m_populationmap->PointIntToNormalized(coords);
                    RoadSection& section = *breedingSection;
                    float curDist = DistanceToSection(coordsNorm, section.start->pos, section.end->pos);
                    if (Approx(curDist, populationDistances.Get(coords).intensity))
                        realSectionPopulation += pp.intensity;
                }
            }
            if (realSectionPopulation <= final)
            {
                expansionType = EXPANSION_TYPE::LINKING;
                continue;
            }
            // debug("2.49");
        }
        else
        {
            // debug("2.5");
            linkingPhase = LINKING_PHASE::NONE;
            breedingNode = nullptr;

            for (RoadNode* n: roadNodes)
            {
                if (!n->IsEndNode())
                    continue;

                if (std::find(permitedEndNodes.begin(), permitedEndNodes.end(), n) != permitedEndNodes.end())
                {
                    continue;
                }

                breedingNode = n;
                break;
            }
            if (breedingNode == nullptr)
                break;

            // debug("2_5_2", breedingNode->pos.x, breedingNode->pos.y);

            if (linkingPhase == LINKING_PHASE::NONE)
            {
                //try merge to crossroad
                RoadNode* nearest = nullptr;
                for (RoadNode* other: roadNodes)
                {
                    if (other == breedingNode)
                        continue;

                    bool linked = false;
                    for (RoadSection* otherSec: other->edges)
                    {
                        if (otherSec->start == breedingNode || otherSec->end == breedingNode)
                        {
                            linked = true;
                            break;
                        }
                    }
                    if (linked)
                        continue;

                    if (other->edges.size() > 3)
                        continue;

                    if (nearest == nullptr)
                        nearest = other;
                    else
                    {
                        if (glm::length(nearest->pos - breedingNode->pos) > glm::length(breedingNode->pos - other->pos))
                            nearest = other;
                    }
                }

                if (nearest == nullptr)
                {
                    debug("NEAREST IS NULL!");
                    throw std::exception{};
                }
                glm::vec2 breedingNodeReal = PointMapNormalizedToReal(m_heightmap, breedingNode->pos);
                glm::vec2 nearestReal = PointMapNormalizedToReal(m_heightmap, nearest->pos);

                if (glm::length(breedingNodeReal - nearestReal) < RoadSchemeGenerator::MERGE_TO_CROSSROAD_DISTANCE_LINKING)
                {
                    RoadSection* change = breedingNode->edges[0];
                    RoadNode* newStart = (change->start == breedingNode) ? change->end : change->start;

                    if (IsCrossroadJoinable(nearest, newStart->pos) && IsCrossroadJoinable(newStart, nearest->pos))
                    {
                        linkingPhase = LINKING_PHASE::MERGE_CROSSROAD;
                        newStart->blocked = false;
                        m_roadScheme->graph.RemoveSection(change->start, change->end);
                        delete breedingNode;
                        breedingEndNode = nearest;
                        breedingNode = newStart;
                        breedingSection = newStart->edges[0];
                        UpdateGraph();
                    }
                }
            }
            // debug("2.75");
            if (linkingPhase == LINKING_PHASE::NONE)
            {
                //try merge to midroad
                RoadNode* friendNode = (breedingNode->edges[0]->start == breedingNode) 
                        ? breedingNode->edges[0]->end
                        : breedingNode->edges[0]->start;
                RoadSection* current = breedingNode->edges[0];
                glm::vec2 nearestPoint{100};
                for (RoadSection* other: roadSections)
                {
                    if (other->start == breedingNode || other->end == breedingNode)
                        continue;

                    bool linked = false;
                    for (RoadSection* friendSection: friendNode->edges)
                    {
                        if (friendSection == other)
                        {
                            linked = true;
                            break;
                        }
                    }
                    if (linked)
                        continue;


                    float abs = RoadSchemeGenerator::MIN_DISTANCE_BETWEEN_CROSSROADS;
                    float sectionLenReal = glm::length(PointMapNormalizedToReal(m_heightmap, other->start->pos) -
                        PointMapNormalizedToReal(m_heightmap, other->end->pos));
                    abs = abs * 2.f / sectionLenReal;
                    float prop = RoadSchemeGenerator::HIGHWAY_MIN_SEGMENT_SUBDIVISION_PART;
                    float _border = std::max(prop, abs);

                    if (_border >= 0.5)
                        continue;

                    glm::vec2 p = NearestPointOfSection(breedingNode->pos, other->start->pos, other->end->pos);
                    
                    float part = glm::length(p - other->start->pos) / glm::length(other->end->pos -  other->start->pos);
                    part = clamp(part, _border, 1 - _border);
                    
                    p = LERP(other->start->pos, other->end->pos, part);

                    if (glm::length(nearestPoint) > 10)
                        nearestPoint = p;
                    else
                    {
                        if (glm::length(breedingNode->pos - p) < glm::length(breedingNode->pos - nearestPoint))
                            nearestPoint = p;
                    }
                }

                if (glm::length(nearestPoint) < 9.99)
                {
                    glm::vec2 breedingNodeReal = PointMapNormalizedToReal(m_heightmap, breedingNode->pos);
                    glm::vec2 nearestReal = PointMapNormalizedToReal(m_heightmap, nearestPoint);
                    if (glm::length(breedingNodeReal - nearestReal) < RoadSchemeGenerator::MERGE_TO_MIDROAD_DISTANCE_LINKING)
                    {
                        linkingPhase = LINKING_PHASE::MERGE_MIDROAD;
                        if (friendNode->edges.size() < 2)
                        {
                            debug("TOO FEW EDGES");
                            throw std::exception{};
                        }
                        breedingSection = friendNode->edges[0] == current 
                            ? friendNode->edges[1]
                            : friendNode->edges[0];
                        // debug("2_8");
                        // float debug_m = 10.f;
                        // for (RoadSection* sectionPtr: roadSections)
                        // {
                        //     RoadSection& section = *sectionPtr;

                        //     // std::pair<glm::vec2, bool> intersection = IntersectSegments(
                        //     //         section.start->pos, section.end->pos, breedingNode->pos, breedingEndPos
                        //     //     );
                        //     float dd = DistanceToSection(nearestPoint, section.start->pos, section.end->pos);
                        //     debug_m = std::min(debug_m, dd);
                        // }
                        // debug("2_8_1", debug_m);
                        // debug_m = 10.f;
                        m_roadScheme->graph.RemoveNode(breedingNode);
                        breedingNode = friendNode;
                        breedingEndPos = nearestPoint;
                        friendNode->blocked = false;
                        UpdateGraph();
                        
                        // for (RoadSection* sectionPtr: roadSections)
                        // {
                        //     RoadSection& section = *sectionPtr;

                        //     // std::pair<glm::vec2, bool> intersection = IntersectSegments(
                        //     //         section.start->pos, section.end->pos, breedingNode->pos, breedingEndPos
                        //     //     );
                        //     float dd = DistanceToSection(breedingEndPos, section.start->pos, section.end->pos);
                        //     debug_m = std::min(debug_m, dd);
                        // }
                        // debug("2_8_2", debug_m);
                    }
                }
            }
            // debug("2.85");
            if (linkingPhase == LINKING_PHASE::NONE)
            {
                // debug("2_95");
                //try extend 
                if (std::find(deniedExtensions.begin(), deniedExtensions.end(), breedingNode) == deniedExtensions.end())
                {
                    // debug("2_99");
                    linkingPhase = LINKING_PHASE::EXTEND_DEADEND;
                    breedingSection = breedingNode->edges.at(0);
                }
            }
            // debug("3");
            if (linkingPhase == LINKING_PHASE::NONE)
            {
                //check for correct deadend
                if (std::find(deniedDeadendNodes.begin(), deniedDeadendNodes.end(), breedingNode) == deniedDeadendNodes.end())
                {
                    float radius = RoadSchemeGenerator::DEADEND_SECTION_BIOME_RADIUS_REAL;
                    radius /= GetMapCellRealSize(m_biomemap);
                    BiomePoint bp = m_biomemap->Get(m_biomemap->PointNormalizedToInt(breedingNode->pos));
                    if (bp.biome != Biome::CITY || bp.intensity < radius)
                    {
                        linkingPhase = LINKING_PHASE::CHECKING_FRINGE;
                        breedingSection = breedingNode->edges[0];
                    }
                }
            }
            // debug("4");
            if (linkingPhase == LINKING_PHASE::NONE)
            {
                // permitedEndNodes.push_back(breedingNode);
                if (std::find(deniedDeadendNodes.begin(), deniedDeadendNodes.end(), breedingNode) != deniedDeadendNodes.end())
                {
                    std::remove(deniedDeadendNodes.begin(), deniedDeadendNodes.end(), breedingNode);
                }
                if (std::find(deniedExtensions.begin(), deniedExtensions.end(), breedingNode) != deniedExtensions.end())
                {
                    std::remove(deniedExtensions.begin(), deniedExtensions.end(), breedingNode);
                }
                // debug("4_5");
                m_roadScheme->graph.RemoveNode(breedingNode);
                continue;
            }
        }
        // debug("5", (int)linkingPhase);
        // if (expansionType == EXPANSION_TYPE::LINKING)
        //     debug(breedingNode->pos);

        std::array<RoadNode*, 2> points{breedingSection->start, breedingSection->end};
        float proportion = 0.f;
        if (expansionType == EXPANSION_TYPE::GROWTH)
        {
            std::array<float, 2> populations{points[0]->data, points[1]->data};
            float segemntLen = glm::length(PointMapNormalizedToReal(m_heightmap, points[0]->pos) - 
                PointMapNormalizedToReal(m_heightmap, points[1]->pos));
            float proportionBorder = RoadSchemeGenerator::HIGHWAY_MIN_SEGMENT_SUBDIVISION_PART;
            float absoluteBorder = RoadSchemeGenerator::MIN_DISTANCE_BETWEEN_CROSSROADS;
            proportion = populations[1] / (populations[0] + populations[1]);

            float newProportion = proportion;
            if (proportion < proportionBorder || proportion > 1 - proportionBorder)
                newProportion = roundf(proportion);
            else if (proportion * segemntLen < absoluteBorder || (1 - proportion) * segemntLen < absoluteBorder)
                newProportion = roundf(proportion);

            // debug("5_1", newProportion);
            if (breedingSection->blocked || breedingSection->start->blocked || breedingSection->end->blocked)
            {
                if (Approx(newProportion, 0.f))
                {
                    if (breedingSection->start->blocked)
                    {
                        if (breedingSection->blocked)
                        {
                            newProportion = 1.f;
                        }
                        else
                        {
                            newProportion = std::max(proportionBorder, absoluteBorder / segemntLen);

                            if (proportion < proportionBorder || 
                                proportion > 1 - proportionBorder ||
                                proportion * segemntLen < absoluteBorder || 
                                (1 - proportion) * segemntLen < absoluteBorder)
                            {
                                breedingSection->blocked = true;
                                continue;
                            }
                        }
                    }
                }
                else if (Approx(newProportion, 1.f))
                {
                    if (breedingSection->end->blocked)
                    {
                        if (breedingSection->blocked)
                        {
                            newProportion = 0.f;
                        }
                        else
                        {
                            newProportion = std::min(1.f - proportionBorder, 1.f - absoluteBorder / segemntLen);

                            if (proportion < proportionBorder || 
                                proportion > 1 - proportionBorder ||
                                proportion * segemntLen < absoluteBorder || 
                                (1 - proportion) * segemntLen < absoluteBorder)
                            {
                                breedingSection->blocked = true;
                                continue;
                            }
                        }
                    }
                }
                else
                {
                    if (breedingSection->blocked)
                        newProportion = roundf(proportion);
                }
            }
            proportion = newProportion;
        }
        else if (expansionType == EXPANSION_TYPE::LINKING && 
            (linkingPhase == LINKING_PHASE::MERGE_CROSSROAD ||
            linkingPhase == LINKING_PHASE::MERGE_MIDROAD || 
            linkingPhase == LINKING_PHASE::CHECKING_FRINGE || 
            linkingPhase == LINKING_PHASE::EXTEND_DEADEND))
        {
            proportion = (breedingSection->start == breedingNode) ? 0.f : 1.f;
        }
        // debug("6");
        enum class START_TYPE
        {
            CROSSROAD,
            MIDDLEROAD
        };

        glm::vec2 newSegmentStartPos;
        START_TYPE startType;
        glm::vec2 newSegmentMapDirection;
        if (expansionType == EXPANSION_TYPE::GROWTH)
        {
            if (Approx(proportion, 0.f) || Approx(proportion, 1.f))
            {
                int id = int(roundf(proportion));
                newSegmentMapDirection = breedingSection->StraightDir(points[1 - id]);
                newSegmentStartPos = points[id]->pos;
                startType = START_TYPE::CROSSROAD;
            }
            else
            {
                newSegmentStartPos = LERP(points[0]->pos, points[1]->pos, proportion);
                startType = START_TYPE::MIDDLEROAD;
                newSegmentMapDirection = LeftNormal(breedingSection->StraightDir());
            }
        }
        else if (expansionType == EXPANSION_TYPE::LINKING && 
            (linkingPhase == LINKING_PHASE::MERGE_CROSSROAD ||
            linkingPhase == LINKING_PHASE::MERGE_MIDROAD || 
            linkingPhase == LINKING_PHASE::CHECKING_FRINGE || 
            linkingPhase == LINKING_PHASE::EXTEND_DEADEND))
        {
            // debug("7");
            newSegmentStartPos = breedingNode->pos;
            startType = START_TYPE::CROSSROAD;
            if (linkingPhase == LINKING_PHASE::MERGE_CROSSROAD)
                newSegmentMapDirection = breedingEndNode->pos - breedingNode->pos;
            else if (linkingPhase == LINKING_PHASE::MERGE_MIDROAD)
                newSegmentMapDirection = breedingEndPos - breedingNode->pos;
            else
                newSegmentMapDirection = (breedingSection->start == breedingNode)
                    ?   breedingSection->start->pos - breedingSection->end->pos
                    :   breedingSection->end->pos - breedingSection->start->pos;
        }
        // debug("8");
        glm::vec2 newSegmentStartPosReal = PointMapNormalizedToReal(m_heightmap, newSegmentStartPos);
        std::vector<std::pair<glm::vec2, float>> possibleRayEndings;
        bool extendingDeadend = expansionType == EXPANSION_TYPE::LINKING && linkingPhase == LINKING_PHASE::EXTEND_DEADEND;
        if 
        (
            expansionType == EXPANSION_TYPE::GROWTH || 
            (expansionType == EXPANSION_TYPE::LINKING && 
                (linkingPhase == LINKING_PHASE::CHECKING_FRINGE || linkingPhase == LINKING_PHASE::EXTEND_DEADEND)
            )
        )
        {
            float newSegmentLengthReal;
            if (extendingDeadend)
                newSegmentLengthReal = RoadSchemeGenerator::MAX_GENERATION_SEGMENT_LENGTH;
            else
            {
                newSegmentLengthReal = RoadSchemeGenerator::MAX_GENERATION_SEGMENT_LENGTH - RoadSchemeGenerator::MIN_GENERATION_SEGMENT_LENGTH;
                newSegmentLengthReal = rnd() * newSegmentLengthReal + RoadSchemeGenerator::MIN_GENERATION_SEGMENT_LENGTH;
            }
            int raysCount = RoadSchemeGenerator::HIGHWAY_RAYS_PER_SEGMENT;
            float deltaAngleRadians = 2 * M_PI / (raysCount - 1);
            float newSectionEfficiency = -1;
            // debug("9", extendingDeadend);
            for (int i = 0; i < raysCount; i++)
            {
                glm::vec2 forward = newSegmentMapDirection;
                glm::vec2 left = LeftNormal(newSegmentMapDirection);
                float dispersion = (extendingDeadend)
                    ? glm::radians(RoadSchemeGenerator::EXTEND_DEADEND_ANGLE_DISPERSION_DEGREES)
                    : glm::radians(359.f / 2.f - RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES());
                float angleRadians = deltaAngleRadians * i;
                angleRadians = LERP(-dispersion, dispersion, angleRadians / (2 * M_PI));
                glm::vec2 curDirection = cosf(angleRadians) * forward + sinf(angleRadians) * left;
                
                glm::vec2 curDirectionReal = PointMapNormalizedToReal(m_heightmap, glm::vec2{0.5,0.5} + curDirection * 0.5f) -
                    PointMapNormalizedToReal(m_heightmap, glm::vec2{0.5,0.5});
                curDirectionReal = glm::normalize(curDirectionReal);
                float currentEfficiency = 0.f;
                std::vector<std::pair<float, float>> steps;
                float reachedDistance = RoadSchemeGenerator::RAY_GENERATION_STEP_REAL;
                while(reachedDistance < newSegmentLengthReal)
                {
                    float pointEfficieny = 0.f;
                    glm::vec2 pointPosNorm = PointRealToMapNormalized(m_heightmap, 
                            newSegmentStartPosReal + curDirectionReal * reachedDistance
                        );

                    if (!m_populationmap->IsPointValid_Normalized(pointPosNorm))
                        break;

                    vec2Int pointPosInt = m_populationmap->PointNormalizedToInt(pointPosNorm);

                    if (m_populationmap->Get(pointPosInt).population != Population::POPULATED)
                    {
                        break;
                    }

                    pointEfficieny = m_populationmap->Get(pointPosInt).intensity;
                    pointEfficieny *= populationDistances.Get(pointPosInt).intensity;
                    pointEfficieny /= glm::length(pointPosNorm - newSegmentStartPos);
                    currentEfficiency += pointEfficieny;

                    steps.emplace_back(reachedDistance, currentEfficiency);
                    reachedDistance += RoadSchemeGenerator::RAY_GENERATION_STEP_REAL;
                }
                reachedDistance -= RoadSchemeGenerator::RAY_GENERATION_STEP_REAL;
                // debug("9_1", steps.size());
                //AFTER REWORK THIS WILL CHANGE
                if (steps.size() == 0)
                    continue;

                int bestStepInd = -1;
                float bestStepEfficiency = -1.f;
                for (int j = 0; j < steps.size(); j++)
                {
                    if (extendingDeadend)
                    {
                        if (steps.at(j).first < RoadSchemeGenerator::MIN_DISTANCE_BETWEEN_CROSSROADS)
                        {
                            continue;
                        }
                    }
                    else
                    {
                        if (steps.at(j).first < RoadSchemeGenerator::MIN_GENERATION_SEGMENT_LENGTH)
                        {
                            continue;
                        }
                    }
                    float stepEfficiency = steps.at(j).second / (j + 1);

                    float ratio = RoadSchemeGenerator::RAY_GENERATION_LONG_PROMOTION_RATIO;
                    stepEfficiency *= powf(steps.at(j).first / steps.at(steps.size() - 1).first, ratio);

                    if (stepEfficiency > bestStepEfficiency)
                    {
                        bestStepInd = j;
                        bestStepEfficiency = stepEfficiency;
                    }
                }
                
                if (bestStepInd  < 0)
                    continue;

                if (extendingDeadend)
                    bestStepInd = steps.size() - 1;

                glm::vec2 currentRayEnding = PointRealToMapNormalized(m_heightmap, 
                        newSegmentStartPosReal + curDirectionReal * steps.at(bestStepInd).first
                    );

                if (m_heightmap->IsPointValid_Normalized(currentRayEnding))
                    possibleRayEndings.push_back(std::pair<glm::vec2, float>{currentRayEnding, bestStepEfficiency});
                else
                    continue;
            }
            
            std::sort(possibleRayEndings.begin(), possibleRayEndings.end(),
                [] (std::pair<glm::vec2, float> a, std::pair<glm::vec2, float> b) 
                { 
                    return a.second > b.second;
                }
            );
            // debug("10");
        }
        else
        {
            if (linkingPhase == LINKING_PHASE::MERGE_CROSSROAD)
                possibleRayEndings.emplace_back(breedingEndNode->pos, 1.f);
            else if (linkingPhase == LINKING_PHASE::MERGE_MIDROAD)
                possibleRayEndings.emplace_back(breedingEndPos, 1.f);
            else
            {
                debug("WRONG LINKING PHASE!");
                throw std::exception{};
            }
        }

        enum class RESULT
        {
            PROCESSING,
            FAILED,
            REPEAT,
            SUCCESS
        };
        enum class END_TYPE
        {
            EMPTY = 0,
            JOIN_CROSS = 1,
            CREATE_CROSS = 2,
            UNCLASSIFIED_CROSS,
        };
        
        struct EndInfo
        {
            RESULT result = RESULT::PROCESSING;
            END_TYPE endType = END_TYPE::EMPTY;
            RoadNode* endCross = nullptr;
            RoadSection* endSection = nullptr;
        };
        std::vector<std::pair<EndInfo, glm::vec2>> possibleExtensions;

        // debug("11", possibleRayEndings.size());
        EndInfo endInfo;
        while(endInfo.result != RESULT::SUCCESS && possibleRayEndings.size() > 0)
        {
            // debug("11_05", (int)endInfo.result);
            if (endInfo.result == RESULT::REPEAT)
                endInfo.result = RESULT::PROCESSING;

            glm::vec2 curEnd;
            if (endInfo.result == RESULT::PROCESSING)
                curEnd = possibleRayEndings.at(0).first;
            else
            {
                possibleRayEndings.erase(possibleRayEndings.begin());
                endInfo.result = RESULT::PROCESSING;
                endInfo.endType = END_TYPE::EMPTY;
                endInfo.endCross = nullptr;
                endInfo.endSection = nullptr;
                continue;
            }
            // debug("11_1", curEnd);
            float curLength = glm::length(PointMapNormalizedToReal(m_heightmap, curEnd) - newSegmentStartPosReal);
            if (curLength <= RoadSchemeGenerator::MIN_DISTANCE_BETWEEN_CROSSROADS)
            {
                endInfo.result = RESULT::FAILED;
                continue;
            }
            // debug("11_2", (int)endInfo.endType);
            if (endInfo.endType == END_TYPE::EMPTY)
            {
                float mergeDist = RoadSchemeGenerator::MERGE_TO_CROSSROAD_DISTANCE_GROWTH;
                for (RoadNode* n: roadNodes)
                {
                    glm::vec2 d = PointMapNormalizedToReal(m_heightmap, n->pos);
                    d -= PointMapNormalizedToReal(m_heightmap, curEnd);
                    if (glm::length(d) < mergeDist)
                    {
                        if (IsCrossroadJoinable(n, newSegmentStartPos))
                        {
                            RoadNode* _newSegmentStart = nullptr;
                            if (startType == START_TYPE::CROSSROAD)
                            {
                                int id = int(roundf(proportion));
                                _newSegmentStart = points[id];
                            }
                            if (startType != START_TYPE::CROSSROAD || IsCrossroadJoinable(_newSegmentStart, n->pos))
                            {
                                possibleRayEndings.at(0).first = n->pos;
                                endInfo.result = RESULT::REPEAT;
                                endInfo.endType = END_TYPE::JOIN_CROSS;
                                endInfo.endSection = nullptr;
                                endInfo.endCross = n;
                                break;
                            }
                        }
                    }
                }
                if (endInfo.result == RESULT::REPEAT)
                    continue;
            }
            // debug("11_3", (int)endInfo.endType);
            for (RoadSection* sectionPtr: roadSections)
            {
                RoadSection& section = *sectionPtr;

                std::pair<glm::vec2, bool> intersection = IntersectSegments(
                        section.start->pos, section.end->pos, newSegmentStartPos, curEnd
                    );

                // debug("11_31", section.start->pos);
                // debug("11_32", section.end->pos);
                // debug("11_33", newSegmentStartPos);
                // debug("11_34", curEnd);
                // debug("11_35", intersection.first);
                // debug("11_36", intersection.second);
                if (intersection.second == true)
                {
                    if (!Approx(intersection.first, newSegmentStartPos))
                    {
                        bool shouldRepeat = false;

                        if (endInfo.endType == END_TYPE::EMPTY)
                            shouldRepeat = true;
                        else 
                        {
                            if (endInfo.endType == END_TYPE::UNCLASSIFIED_CROSS || endInfo.endType == END_TYPE::CREATE_CROSS)
                                shouldRepeat = !(endInfo.endSection == sectionPtr);
                            else
                            {
                                shouldRepeat = !(section.end == endInfo.endCross || section.start == endInfo.endCross);
                            }

                            if (glm::length(newSegmentStartPos - intersection.first) >= glm::length(newSegmentStartPos - possibleRayEndings.at(0).first))
                                shouldRepeat = false;
                        }

                        if (shouldRepeat)
                        {
                            possibleRayEndings.at(0).first = intersection.first;
                            endInfo.result = RESULT::REPEAT;
                            endInfo.endType = END_TYPE::UNCLASSIFIED_CROSS;
                            endInfo.endSection = sectionPtr;
                            endInfo.endCross = nullptr;
                            break;
                        }
                    }
                }
            }
            if (endInfo.result == RESULT::REPEAT)
                continue;
            // debug("11_4", (int)endInfo.endType);
            if (endInfo.endType == END_TYPE::UNCLASSIFIED_CROSS)
            {
                float mergeDist = RoadSchemeGenerator::MERGE_TO_CROSSROAD_DISTANCE_GROWTH;
                for (RoadNode* n: endInfo.endSection->GetEdges())
                {
                    glm::vec2 d = PointMapNormalizedToReal(m_heightmap, n->pos);
                    d -= PointMapNormalizedToReal(m_heightmap, curEnd);
                    if (glm::length(d) < mergeDist)
                    {
                        if (IsCrossroadJoinable(n, newSegmentStartPos))
                        {
                            RoadNode* _newSegmentStart = nullptr;
                            if (startType == START_TYPE::CROSSROAD)
                            {
                                int id = int(roundf(proportion));
                                _newSegmentStart = points[id];
                            }
                            if (startType != START_TYPE::CROSSROAD || IsCrossroadJoinable(_newSegmentStart, n->pos))
                            {
                                possibleRayEndings.at(0).first = n->pos;
                                endInfo.result = RESULT::REPEAT;
                                endInfo.endType = END_TYPE::JOIN_CROSS;
                                endInfo.endSection = nullptr;
                                endInfo.endCross = n;
                                break;
                            }
                        }
                        else if (glm::length(d) < RoadSchemeGenerator::INTERSECTING_ROADS_DISTANCE_REAL)
                        {
                            endInfo.result = RESULT::FAILED;
                            break;
                        }
                    }
                }
                if (endInfo.result == RESULT::REPEAT)
                    continue;
            }
            // debug("11_5", (int)endInfo.endType);
            if (endInfo.endType == END_TYPE::UNCLASSIFIED_CROSS)
            {
                if (endInfo.endSection == nullptr)
                {
                    debug("END SECTION NEEDED BUT IT IS NULL");
                    throw std::exception{};
                }

                auto edges = endInfo.endSection->GetEdges();

                float sectionLenReal = glm::length(PointMapNormalizedToReal(m_heightmap, edges[0]->pos) -
                    PointMapNormalizedToReal(m_heightmap, edges[1]->pos));
                // debug("11_55", sectionLenReal);
                // debug("11_56", newSegmentStartPos.x, newSegmentStartPos.y);
                // debug("11_57", curEnd.x, curEnd.y);
                // debug("11_58", edges[0]->pos.x, edges[0]->pos.y);
                // debug("11_59", edges[1]->pos.x, edges[1]->pos.y);
                float minProportion = RoadSchemeGenerator::MIN_DISTANCE_BETWEEN_CROSSROADS / sectionLenReal;
                minProportion = std::max(minProportion, RoadSchemeGenerator::HIGHWAY_MIN_SEGMENT_SUBDIVISION_PART);

                if (minProportion > 0.5f)
                {
                    endInfo.result = RESULT::FAILED;
                    continue;
                }

                float endProportion = 
                    glm::length(edges[0]->pos - curEnd) / glm::length(edges[0]->pos - edges[1]->pos);
                endProportion = (endProportion < 0.5f) ? 
                    std::max(minProportion, endProportion) :
                    std::min(1 - minProportion, endProportion);
                glm::vec2 newEndPos = LERP(edges[0]->pos, edges[1]->pos, endProportion);

                float minAngleRadians = M_PI + 0.1f; //min angle to other roads in this crossroad
                for (RoadNode* node: edges)
                {
                    glm::vec2 dirNew = glm::normalize(newSegmentStartPos - newEndPos);
                    glm::vec2 dirEdge = glm::normalize(node->pos - newEndPos);
                    float angle = glm::dot(dirNew, dirEdge);
                    angle = SafeAcosf(angle);
                    minAngleRadians = std::min(minAngleRadians, angle);
                }
                float borderRadians = RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES();
                borderRadians = glm::radians(borderRadians);
                if (borderRadians > minAngleRadians)
                {
                    endInfo.result = RESULT::FAILED;
                    continue;
                }

                if (startType == START_TYPE::CROSSROAD)
                {
                    RoadNode* _newSegmentStart = nullptr;
                    int id = int(roundf(proportion));
                    _newSegmentStart = points[id];

                    if (!IsCrossroadJoinable(_newSegmentStart, newEndPos))
                    {
                        endInfo.result = RESULT::FAILED;
                        continue;
                    }
                }

                // debug("11_59_3", glm::degrees(minAngleRadians));
                possibleRayEndings.at(0).first = newEndPos;
                endInfo.result = RESULT::REPEAT;
                endInfo.endType = END_TYPE::CREATE_CROSS;
                endInfo.endSection = endInfo.endSection;
                endInfo.endCross = nullptr;
                continue;
            }
            // debug("11_6", (int)endInfo.endType);
            if (extendingDeadend)
            {
                // debug("11_8");
                if (endInfo.result == RESULT::PROCESSING)
                {
                    endInfo.result = RESULT::FAILED;
                    RoadNode* friendNode = (breedingNode->edges[0]->start == breedingNode)
                            ? breedingNode->edges[0]->end
                            : breedingNode->edges[0]->start;
                    
                    {
                        float existLen = breedingSection->GetLength();
                        float extensionLen = glm::length(breedingNode->pos - possibleRayEndings.at(0).first);
                        float border = RoadSchemeGenerator::EXTEND_DEADEND_MAX_LENGTH_RATIO;
                        if (extensionLen / existLen > border)
                        {
                            continue;
                        }
                    } 

                    if (endInfo.endType == END_TYPE::CREATE_CROSS)
                    {
                        if (endInfo.endSection->start != friendNode && endInfo.endSection->end != friendNode)
                        {
                            possibleExtensions.emplace_back(endInfo, possibleRayEndings.at(0).first);
                        }
                    }
                    else if (endInfo.endType == END_TYPE::JOIN_CROSS)
                    {
                        bool adjacent = false;
                        for (RoadSection* s: friendNode->edges)
                        {
                            if (s->start == endInfo.endCross || s->end == endInfo.endCross)
                            {
                                adjacent = true;
                                break;
                            }
                        }
                        if (!adjacent)
                            possibleExtensions.emplace_back(endInfo, possibleRayEndings.at(0).first);
                    }
                    // debug("11_9", (int)endInfo.endType, possibleExtensions.size());
                }
            }

            if (endInfo.endType == END_TYPE::EMPTY && endInfo.result == RESULT::PROCESSING)
            {
                if (startType == START_TYPE::CROSSROAD)
                {
                    RoadNode* _newSegmentStart = nullptr;
                    int id = int(roundf(proportion));
                    _newSegmentStart = points[id];

                    if (!IsCrossroadJoinable(_newSegmentStart, possibleRayEndings.at(0).first))
                    {
                        endInfo.result = RESULT::FAILED;
                        continue;
                    }
                }
                else if (startType == START_TYPE::MIDDLEROAD)
                {
                    RoadSection* _section = RoadSection::SectionByEdges(points[0], points[1]);
                    if (!IsSectionJoinable(_section, proportion, possibleRayEndings.at(0).first))
                    {
                        endInfo.result = RESULT::FAILED;
                        continue;
                    }
                }
            }

            if (endInfo.result == RESULT::PROCESSING)
                endInfo.result = RESULT::SUCCESS;
        }
        // debug("12");
        if (expansionType == EXPANSION_TYPE::LINKING && linkingPhase == LINKING_PHASE::CHECKING_FRINGE)
        {
            // debug("12_1");
            int s = possibleRayEndings.size();
            float part = 1 - RoadSchemeGenerator::DEADEND_SECTION_BLOCKED_RAYS_MIN_PART;
            if (s <= part * RoadSchemeGenerator::HIGHWAY_RAYS_PER_SEGMENT)
            {
                permitedEndNodes.push_back(breedingNode);
            }
            else
            {
                deniedDeadendNodes.push_back(breedingNode);
            }
        }
        else
        {
            // debug("12_2", possibleRayEndings.size());
            if (possibleRayEndings.size() == 0 && !extendingDeadend)
            {
                if (startType == START_TYPE::MIDDLEROAD)
                    breedingSection->blocked = true;
                else
                {
                    int id = int(roundf(proportion));
                    points[id]->blocked = true;
                }
            }
            else
            {
                if (extendingDeadend && possibleExtensions.size() == 0)
                {
                    //  debug("NO");
                    deniedExtensions.push_back(breedingNode);
                }
                else
                {
                    // debug("13", (int)endInfo.endType, (int)startType);
                    if (extendingDeadend)
                    {
                        // debug("14", possibleExtensions.size());
                        std::sort(
                            possibleExtensions.begin(), 
                            possibleExtensions.end(),
                            [&breedingNode] (const std::pair<EndInfo, glm::vec2>& a, const std::pair<EndInfo, glm::vec2>& b) 
                            { 
                                glm::vec2 aPos, bPos;
                                if (a.first.endType == END_TYPE::JOIN_CROSS)
                                    aPos = a.first.endCross->pos;
                                else if (a.first.endType == END_TYPE::CREATE_CROSS)
                                    aPos = a.second;
                                else
                                {
                                    debug("WRONG END TYPE FOR EXTENSION");
                                    throw std::exception{};
                                }
                                if (b.first.endType == END_TYPE::JOIN_CROSS)
                                    bPos = b.first.endCross->pos;
                                else if (b.first.endType == END_TYPE::CREATE_CROSS)
                                    bPos = b.second;
                                else
                                {
                                    debug("WRONG END TYPE FOR EXTENSION");
                                    throw std::exception{};
                                }
                                return glm::length(aPos - breedingNode->pos) < glm::length(bPos - breedingNode->pos);
                            }
                        );
                        endInfo = possibleExtensions.at(0).first;
                        // debug("14_1", (int)endInfo.endType);
                    }
                    
                    RoadNode* newSegmentStart = nullptr;
                    if (startType == START_TYPE::CROSSROAD)
                    {
                        int id = int(roundf(proportion));
                        newSegmentStart = points[id];
                    }
                    else
                    {
                        if (extendingDeadend)
                        {
                            debug("SHOULD NOT BE STARTING MIDROAD!");
                            throw std::exception{};
                        }
                        newSegmentStart = m_roadScheme->graph.DivideSection(points[0], points[1], proportion);
                        if (proportion < 0.00001 || proportion > 0.9999)
                        {
                            debug("ER 2");
                            throw std::exception{};
                        }
                    }

                    RoadNode* newSegmentEnd = nullptr;
                    if (endInfo.endType == END_TYPE::JOIN_CROSS)
                    {
                        if (endInfo.endCross == nullptr)
                        {
                            debug("END CROSS NEEDED BUT IT IS NULL");
                            throw std::exception{};
                        }
                        newSegmentEnd = endInfo.endCross;
                        if (extendingDeadend)
                        {
                            // debug("14_2_1", newSegmentEnd->pos.x, newSegmentEnd->pos.y);
                        }
                    }
                    else if (endInfo.endType == END_TYPE::CREATE_CROSS)
                    {
                        if (endInfo.endSection == nullptr)
                        {
                            debug("END SECTION NEEDED BUT IT IS NULL");
                            throw std::exception{};
                        }
                        glm::vec2 endPos = (extendingDeadend) 
                            ? possibleExtensions.at(0).second 
                            : possibleRayEndings.at(0).first;
                        float endProportion = glm::length(endInfo.endSection->start->pos - endPos);
                        endProportion  /= glm::length(endInfo.endSection->start->pos - endInfo.endSection->end->pos);
                        if (endProportion < 0.00001 || endProportion > 0.9999)
                        {
                            debug("ER", extendingDeadend);
                            throw std::exception{};
                        }
                        newSegmentEnd = m_roadScheme->graph.DivideSection(
                                endInfo.endSection->start, 
                                endInfo.endSection->end, 
                                endProportion
                            );
                        if (extendingDeadend)
                        {
                            // debug("14_2_2", endProportion);
                        }
                    }
                    else
                    {
                        if (extendingDeadend)
                        {
                            debug("SHOULD NOT BE EMPTY END!");
                            throw std::exception{};
                        }
                        newSegmentEnd = m_roadScheme->graph.AddNode(possibleRayEndings.at(0).first);
                    }
                    // debug("15");
                    // debug("15_1", newSegmentStart->edges.size());
                    // debug("15_2", newSegmentEnd->edges.size());

                    //check for some frequent mistakes
                    if (expansionType == EXPANSION_TYPE::LINKING)
                    {
                        if (linkingPhase == LINKING_PHASE::MERGE_CROSSROAD || linkingPhase == LINKING_PHASE::MERGE_MIDROAD)
                        {
                            if (endInfo.endType == END_TYPE::EMPTY)
                            {
                                debug("NOT CORRECT END TYPE!");
                                throw std::exception{};
                            }
                        }
                    }

                    m_roadScheme->graph.AddSection(newSegmentStart, newSegmentEnd, 0.f);
                    // debug("15_3", newSegmentStart->pos);
                    // debug("15_4", newSegmentStart->edges.size());
                    // debug("15_5", newSegmentEnd->pos);
                    // debug("15_6", newSegmentEnd->edges.size());

                    // debug("16_1", (int)endInfo.endType);
                    // debug("16_2", (int)expansionType);
                    // debug("16_3", (int)linkingPhase);

                    if (newSegmentEnd->edges.size() > 4)
                    {
                        debug("TOO MUCH ROADS");
                        throw std::exception{};
                    }

                    for (RoadNode* rn: m_roadScheme->graph.GetNodes())
                    {
                        float minAngleDegrees = 361.f;
                        for (RoadSection* rs: rn->edges)
                        {
                            for (RoadSection* otherRs: rn->edges)
                            {
                                if (otherRs == rs)
                                    continue;
                                minAngleDegrees = std::min(
                                    minAngleDegrees, 
                                    glm::degrees(SafeAcosf(glm::dot(otherRs->StraightDir(rn), rs->StraightDir(rn))))
                                );
                            }
                        }
                        
                        float border = RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES();
                        if (minAngleDegrees < border)
                        {
                            // debug("ROADS ARE TO CLOSE IN CROSSROAD! angle is ", glm::degrees(alphaRadians), "degrees");
                            // debug("CROSSROAD", crossroad->pos);
                            // for (RoadSection* s: crossroad->edges)
                            // {
                            //     RoadNode* n = (s->start == crossroad)
                            //         ? s->end
                            //         : s->start;
                            //     debug("ROAD", n->pos);
                            // }
                            debug("CLOSE ROADS!", minAngleDegrees, border, rn->pos);
                            throw std::exception{};
                        }
                    }
                }
            }
        }
    }

    //CHECKING LENGTHS
    for (RoadSection* r: m_roadScheme->graph.GetSections())
    {
        if (r->GetLength() < 0.00001)
        {
            debug("WRONG LEN!");
            throw std::exception{};
        }
    }

    RoundOffTurns();
    
    m_roadScheme->SaveAsTexture();

    debug("Done!");
    debug("Generating road mesh!");

    //roadDirection is facing right from foundation
    std::function CreateSliceFromCrossroad
    {
        [this]
        (RoadSection& section, bool isEndNode)
        {
            RoadNode* crossroad = (isEndNode) ? section.end : section.start;
            RoadSlice result = CreateSliceFromCrossroad_Uncorrected(section, isEndNode);
            auto planarInds = RoadSlice::GetPlanarBorderIndexes();
            std::array<glm::vec2, 2> planarRealPoints{
                XZ(result.border[planarInds[0]]), 
                XZ(result.border[planarInds[1]])
            };

            // Correcting slice side edges
            auto turns = GetLeftRightTurns(crossroad, section);
            for (int i = 0; i < 2; i++)
            {
                bool correctingLeft = (i == 0);
                RoadSection* neighb;
                if (correctingLeft) 
                    neighb = (turns[0].size() > 0) ? turns[0][0] : turns[1][turns[1].size() - 1];
                else
                    neighb = (turns[1].size() > 0) ? turns[1][0] : turns[0][turns[0].size() - 1];

                if ((correctingLeft && turns[0].size() != 0) || (!correctingLeft && turns[1].size() != 0))
                {
                    float borderAngle = RoadModelGenerator::CROSSROAD_CORRECTION_MIN_ANGLE_DEGREES;
                    float angle = SafeAcosf(glm::dot(neighb->StraightDir(crossroad), section.StraightDir(crossroad)));
                    angle = glm::degrees(angle);
                    if (angle < borderAngle)
                        continue;
                }
                
                glm::vec2 correctingDirection = correctingLeft ?
                    LeftNormal(section.StraightDir(crossroad)) : 
                    -LeftNormal(section.StraightDir(crossroad));
                glm::vec2 planarNormalizedPoint = PointRealToMapNormalized(m_heightmap, planarRealPoints[0]);
                glm::vec2 correctingPlanarPoint = (glm::dot(planarNormalizedPoint - crossroad->pos, correctingDirection) > 0) ?
                    planarRealPoints[0] : planarRealPoints[1];
                bool curFlag = Approx(glm::length(correctingPlanarPoint - XZ(result.border[planarInds[0]])), 0.f);
                int curShift = (curFlag) ? -1 : 1;
                int curI = (curFlag) ? planarInds[0] : planarInds[1];
                curI += curShift;
                RoadSlice neighbSlice = CreateSliceFromCrossroad_Uncorrected(*neighb, (*neighb->end) == (*crossroad));
                bool neighbFlag = Approx(glm::length(correctingPlanarPoint - XZ(neighbSlice.border[planarInds[0]])), 0.f);
                int neighbShift = (neighbFlag) ? -1 : 1;
                int neighbI = (neighbFlag) ? planarInds[0] : planarInds[1];
                neighbI += neighbShift;
                    
                while (curI >= 0 && curI < result.border.size())
                {
                    result.border[curI] = (neighbSlice.border[neighbI] + result.border[curI]) * 0.5f;
                    curI += curShift;
                    neighbI += neighbShift;
                }
            }
            
            //result
            return result;
        }
    };

    int vertexesAmount = 0;
    int availableAmount = 1;
    const int VERTEX_SIZE_IN_FLOAT = RoadSlice::VERTEX_SIZE_IN_FLOAT;
    GLfloat* generalBuffer = new GLfloat[availableAmount * VERTEX_SIZE_IN_FLOAT];
    const int FLOAT_SIZE_IN_BYTE = sizeof(generalBuffer[0]);
    // debug("B_1");
    const float CONNECTION_STRAIGHT_LEN = Landscape::RoadModelGenerator::SEGMENT_STRAIGHT_REAL_LENGTH;
    const float CONNECTION_ROUND_ANGLE = Landscape::RoadModelGenerator::SEGMENT_ROUND_MAX_ANGLE_DEGREES;
    for (auto sectionPtr : m_roadScheme->graph.GetSections())
    {   
        RoadSection& section = *sectionPtr;
        std::vector<RoadSlice> slices;
        if (!section.IsRoundTurn() && !section.IsStraight())
        {
            debug(section.startDir, section.endDir);
            debug("WRONG ROAD TURN!");
            throw std::exception{};
        }
        
        if (section.IsStraight())
        {
            
            glm::vec2 dir = section.StraightDir();
            RoadSlice startSlice = (section.start->IsCrossroad()) ? 
                CreateSliceFromCrossroad(section, false) : 
                CreateSliceFromSchemePoint(section.start->pos, dir, 1.f);
            RoadSlice endSlice = (section.end->IsCrossroad()) ? 
                CreateSliceFromCrossroad(section, true) : 
                CreateSliceFromSchemePoint(section.end->pos, dir, 1.f);
            float mapSegmentLen = glm::length(PointMapNormalizedToReal(m_heightmap, section.end->pos) - 
                PointMapNormalizedToReal(m_heightmap, section.start->pos));
            float step = CONNECTION_STRAIGHT_LEN / mapSegmentLen;
            
            float current = 0;
            while (current <= 1.f)
            {
                bool skipFlag = false;

                RoadSlice slice;
                if (current == 0.f) 
                    slice = startSlice;
                else if (current == 1.f) 
                    slice = endSlice;
                else
                {
                    glm::vec2 pos = LERP(section.start->pos, section.end->pos, current);
                    slice = CreateSliceFromSchemePoint(pos, dir, 1.f);

                    //Check for not crossing crossroad slice
                    glm::vec2 realStraightDir = section.StraightDir();
                    realStraightDir = PointMapNormalizedToReal(m_heightmap, section.end->pos) - 
                        PointMapNormalizedToReal(m_heightmap, section.start->pos);
                    realStraightDir = glm::normalize(realStraightDir);
                    std::vector<float> cosSigns;

                    glm::vec2 startDir1 = 
                        XZ(slice.GetBoundRealPoints()[0] - startSlice.GetBoundRealPoints()[0]);
                    glm::vec2 startDir2 = 
                        XZ(slice.GetBoundRealPoints()[1] - startSlice.GetBoundRealPoints()[0]);
                    glm::vec2 startSliceNormal = LeftNormal(
                        XZ(startSlice.GetBoundRealPoints()[1] - startSlice.GetBoundRealPoints()[0]));
                    glm::vec2 endDir1 = 
                        XZ(slice.GetBoundRealPoints()[0] - endSlice.GetBoundRealPoints()[0]);
                    glm::vec2 endDir2 = 
                        XZ(slice.GetBoundRealPoints()[1] - endSlice.GetBoundRealPoints()[0]);
                    glm::vec2 endSliceNormal = LeftNormal(
                        XZ(endSlice.GetBoundRealPoints()[1] - endSlice.GetBoundRealPoints()[0]));
                    cosSigns.push_back(glm::dot(realStraightDir, startSliceNormal));
                    cosSigns.push_back(glm::dot(startDir1, startSliceNormal));
                    cosSigns.push_back(glm::dot(startDir2, startSliceNormal));
                    cosSigns.push_back(glm::dot(realStraightDir, endSliceNormal));
                    cosSigns.push_back(-glm::dot(endDir1, endSliceNormal));
                    cosSigns.push_back(-glm::dot(endDir2, endSliceNormal));

                    float sign = cosSigns[0];
                    for (int i = 0; i < cosSigns.size(); i++)
                    {
                        if (sign * cosSigns[i] < 0)
                        {
                            skipFlag = true;
                        }
                    }
                }
                
                if (!skipFlag)
                {
                    slices.push_back(slice);
                }
                
                if (current == 1.0f)
                {
                    break;
                }
                current = clamp(current + step, 0.f, 1.f);
            }
        }
        else if (section.IsRoundTurn())
        {
            if (section.start->IsCrossroad() || section.end->IsCrossroad())
            {
                debug("NO CROSSROADS FOR CIRCLE TURNS!");
                throw std::exception{};
            }
            glm::vec2 middle = (section.start->pos + section.end->pos) * 0.5f;
            float alphaAngleRadians = SafeAcosf(glm::dot(section.startDir, section.StraightDir()));
            float middleLen = glm::length(middle - section.start->pos);
            float normalLen = middleLen * glm::tan(M_PI * 0.5f - alphaAngleRadians);
            float radius = sqrtf(middleLen * middleLen + normalLen * normalLen);
            glm::vec2 middleToCenter = LeftNormal(section.StraightDir());
            if (glm::dot(middleToCenter, section.startDir) > 0.0f)
                middleToCenter *= -1.f;
            glm::vec2 center = middle + middleToCenter * normalLen;
            
            float fullSectionTurnDegrees = glm::degrees(2 * alphaAngleRadians);
            float fullSectionCurveLength = alphaAngleRadians * radius;
            float angleStep = CONNECTION_ROUND_ANGLE / fullSectionTurnDegrees;
            float lengthStep = CONNECTION_STRAIGHT_LEN / fullSectionCurveLength;
            float step = std::min(angleStep, lengthStep);
            float current = 0;
            while (current <= 1.f)
            {
                glm::vec2 leftToStraight = LeftNormal(section.StraightDir());
                float startAngle = SafeAcosf(glm::dot(leftToStraight, section.startDir)); // To straight left
                float endAngle = SafeAcosf(glm::dot(leftToStraight, section.endDir)); //To straight left
                float curAngle = LERP(startAngle, endAngle, current);
                glm::vec2 dir = cosf(curAngle) * leftToStraight + sinf(curAngle) * section.StraightDir();
                glm::vec2 centerToPos = LeftNormal(dir);
                if (glm::dot(centerToPos, middleToCenter) > 0.0f)
                    centerToPos *= -1.f;
                glm::vec2 pos = center + centerToPos * radius;
                RoadSlice slice = CreateSliceFromSchemePoint(pos, dir, 1.f);
                slices.push_back(slice);
                
                if (current == 1.0f)
                {
                    break;
                }
                current = clamp(current + step, 0.f, 1.f);
            }
        }

        for (int i = 0; i < slices.size(); i++)
        {
            RoadSlice* prev = (i > 0) ? &(slices[i-1]) : nullptr;
            RoadSlice* next = (i < slices.size() - 1) ? &(slices[i+1]) : nullptr;
            
            if (prev != nullptr && next != nullptr)
                slices[i].CalculateNormals(prev, next);
            else
            {
                glm::vec3 point;
                if (prev == nullptr && section.start->IsCrossroad())
                {
                    point = X0Z(PointMapNormalizedToReal(m_heightmap, section.start->pos));
                }
                else if (next == nullptr && section.end->IsCrossroad())
                {
                    point = X0Z(PointMapNormalizedToReal(m_heightmap, section.end->pos));
                }
                else
                {
                    auto planarInds = RoadSlice::GetPlanarBorderIndexes();
                    point = (slices[i].border[planarInds[0]] + slices[i].border[planarInds[1]]) * 0.5f;
                }
                glm::vec2 curHeightmapPoint = 
                    m_heightmap->PointNormalizedToIntContinuous(PointRealToMapNormalized(m_heightmap, XZ(point)));
                
                glm::vec3 normal = GetNormal(int(floor(curHeightmapPoint.x)), int(floor(curHeightmapPoint.y)));
                glm::vec3 realFoundation;
                if (section.IsStraight())
                    realFoundation = X0Z(LeftNormal(section.StraightDir()));
                else if (section.IsRoundTurn())
                    realFoundation = X0Z(LeftNormal((prev == nullptr) ? section.startDir : section.endDir)); 
                else
                {
                    debug("WRONG TURN, BROTHER");
                    throw std::exception{};
                }

                realFoundation *= -1;

                slices[i].CalculateNormals(normal, realFoundation);
                for (auto ind : RoadSlice::GetPlanarBorderIndexes())
                {
                    slices[i].normals[ind] = normal;
                }
            }
        }
        
        for (int i = 1; i < slices.size(); i++)
        {
            glm::vec2 roadDir = (section.IsStraight()) ? section.StraightDir() : glm::vec2{0};
            int extraVertexes = RoadSlice::GetConnectionVertexesCount();
            // int extraVertexes = RoadSlice::GetSliceVertexesCount();
            UpdateBufferCapacity(vertexesAmount, extraVertexes, availableAmount, generalBuffer);
            RoadSlice::AddRoadConnectionToBuffer(
                slices[i], slices[i-1], generalBuffer + vertexesAmount * VERTEX_SIZE_IN_FLOAT, roadDir);
            // slice.AddFacesToBuffer(generalBuffer + vertexesAmount * VERTEX_SIZE_IN_FLOAT);
            vertexesAmount += extraVertexes;
        }
    }
    // debug("B_2");
    for (auto nodePtr : m_roadScheme->graph.GetNodes())
    {   
        RoadNode& node = *nodePtr;
        if (!node.IsCrossroad())
            continue;

        glm::vec3 nodeRealPos = X0Z(PointMapNormalizedToReal(m_heightmap,node.pos));
        std::vector<std::array<glm::vec3, 2>> sectionsEdgeCenter;
        auto planarInds = RoadSlice::GetPlanarBorderIndexes();

        for (RoadSection* sectionPtr : node.edges)
        {
            std::array<glm::vec3, 2> cur;

            RoadSection& section = *sectionPtr;
            RoadSlice slice = CreateSliceFromCrossroad(section, (*section.end) == (node));
            std::array<glm::vec3, 2> planarPoints{slice.border[planarInds[0]], slice.border[planarInds[1]]};
            glm::vec3 sectionCenter = (planarPoints[0] + planarPoints[1]) * 0.5f;
            cur[1] = sectionCenter;
            
            //orientating clockwise
            glm::vec2 dirOne = glm::normalize(XZ(planarPoints[0]) - XZ(nodeRealPos));
            glm::vec2 dirCenter = glm::normalize(XZ(sectionCenter) - XZ(nodeRealPos));
            if (glm::dot(LeftNormal(dirCenter), dirOne) > 0.f)
                cur[0] = planarPoints[0];
            else
                cur[0] = planarPoints[1];

            sectionsEdgeCenter.push_back(cur);
        }

        //Sorting into sequence
        for (int i = 0; i < sectionsEdgeCenter.size() - 2; i++)
        {
            auto cur = sectionsEdgeCenter[i];
            glm::vec2 curLeft = LeftNormal(XZ(glm::normalize(cur[1] - nodeRealPos)));
            int nearestLeftInd = i + 1;
            float minAngleRadians = 2 * M_PI;
            for (int j = nearestLeftInd; j < sectionsEdgeCenter.size(); j++)
            {
                float curAngle = SafeAcosf(glm::dot(
                            glm::normalize(XZ(cur[1] - nodeRealPos)), 
                            glm::normalize(XZ(sectionsEdgeCenter[j][1] - nodeRealPos))
                        )
                    );
                if (glm::dot(curLeft, glm::normalize(XZ(sectionsEdgeCenter[j][1] - nodeRealPos))) < 0)
                    curAngle = 2 * M_PI - curAngle;
                
                if (curAngle < minAngleRadians)
                {
                    minAngleRadians = curAngle;
                    nearestLeftInd = j;
                }
            }
            auto temp = sectionsEdgeCenter[i + 1];
            sectionsEdgeCenter[i + 1] = sectionsEdgeCenter[nearestLeftInd];
            sectionsEdgeCenter[nearestLeftInd] = temp;
        }

        for (auto edgeCenter : sectionsEdgeCenter)
        {
            nodeRealPos.y += edgeCenter[1].y / sectionsEdgeCenter.size();
        }

        {
            int size = sectionsEdgeCenter.size();
            bool centerFlag = false;
            int i = 0;
            while (i < size)
            {
                int extraVertexes = RoadSlice::TRIANGLE_SIZE_IN_VERTEX;
                UpdateBufferCapacity(vertexesAmount, extraVertexes, availableAmount, generalBuffer);
                std::array<glm::vec3, 3> triangle;
                std::array<glm::vec3, 3> normals;
                if (centerFlag)
                {
                    triangle = std::array<glm::vec3, 3>
                    {
                        sectionsEdgeCenter[i][1], 
                        sectionsEdgeCenter[(i + 1) % size][1], 
                        nodeRealPos
                    };
                }
                else
                {
                    triangle = std::array<glm::vec3, 3>
                    {
                        sectionsEdgeCenter[i][1], 
                        sectionsEdgeCenter[(i + 1) % size][1], 
                        sectionsEdgeCenter[i][0]
                    };
                }

                GLfloat* ptr = generalBuffer  + vertexesAmount * VERTEX_SIZE_IN_FLOAT;
                for (int j = 0; j < triangle.size(); j++)
                {
                    memcpy(ptr, glm::value_ptr(triangle[j]), sizeof(triangle[0]));
                    ptr += (sizeof(triangle[0]) / sizeof(triangle[0].x));

                    glm::vec2 curHeightmapPoint = m_heightmap->PointNormalizedToIntContinuous(
                        PointRealToMapNormalized(m_heightmap, XZ(nodeRealPos)));
                    glm::vec3 normal = GetNormal(int(floor(curHeightmapPoint.x)), int(floor(curHeightmapPoint.y)));
                    memcpy(ptr, glm::value_ptr(normal), sizeof(normal));
                    ptr += (sizeof(normal) / sizeof(normal.x));

                    glm::vec2 asphaltTexCoords = XZ(triangle[j]);
                    memcpy(ptr, glm::value_ptr(asphaltTexCoords), sizeof(asphaltTexCoords));
                    ptr += (sizeof(asphaltTexCoords) / sizeof(asphaltTexCoords.x));


                    float markingTexCoords = 0;
                    glm::vec4 bezierAnchorsCoords{0};
                    glm::vec4 bezierHandlesCoords{0};

                    if (!centerFlag)
                    {
                        markingTexCoords = 1.3;
                        if (j == 2)
                        {
                            bezierAnchorsCoords = glm::vec4 {XZ(triangle[0]), XZ(triangle[1])};
                            bezierHandlesCoords = glm::vec4 {XZ(nodeRealPos), XZ(triangle[2])};
                        }
                    }

                    memcpy(ptr, (void*)(&markingTexCoords), sizeof(markingTexCoords));
                    ptr += (sizeof(markingTexCoords) / sizeof(markingTexCoords));
                    
                    memcpy(ptr, glm::value_ptr(bezierAnchorsCoords), sizeof(bezierAnchorsCoords));
                    ptr += (sizeof(bezierAnchorsCoords) / sizeof(bezierAnchorsCoords.x));

                    memcpy(ptr, glm::value_ptr(bezierHandlesCoords), sizeof(bezierHandlesCoords));
                    ptr += (sizeof(bezierHandlesCoords) / sizeof(bezierHandlesCoords.x));
                }
                vertexesAmount += extraVertexes;   

                i += (centerFlag) ? 1 : 0;
                centerFlag = !centerFlag;
            }
        }
    }
    // debug("B_3");
    m_road = new RoadModel();
    m_road->vertexesNum = vertexesAmount;
    glBindVertexArray(m_road->VAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_road->VBO);
    glBufferData(
        GL_ARRAY_BUFFER, 
        m_road->vertexesNum * FLOAT_SIZE_IN_BYTE * VERTEX_SIZE_IN_FLOAT, 
        generalBuffer, 
        GL_STATIC_DRAW);

    glVertexAttribPointer(
        0, 
        RoadSlice::POSITION_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)0
    );
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(
        1, 
        RoadSlice::NORMAL_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)(RoadSlice::POSITION_SIZE_IN_FLOAT * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(1);
    
    glVertexAttribPointer(
        2, 
        RoadSlice::ASPHALT_TEX_COORDS_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)((RoadSlice::POSITION_SIZE_IN_FLOAT + RoadSlice::NORMAL_SIZE_IN_FLOAT) * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(2);
    
    glVertexAttribPointer(
        3, 
        RoadSlice::MARKING_TEX_COORDS_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)((RoadSlice::POSITION_SIZE_IN_FLOAT + 
            RoadSlice::ASPHALT_TEX_COORDS_SIZE_IN_FLOAT +
            RoadSlice::NORMAL_SIZE_IN_FLOAT) * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(3);
    
    glVertexAttribPointer(
        4, 
        RoadSlice::BEZIER_ANCHORS_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)((RoadSlice::POSITION_SIZE_IN_FLOAT + 
            RoadSlice::NORMAL_SIZE_IN_FLOAT +
            RoadSlice::ASPHALT_TEX_COORDS_SIZE_IN_FLOAT +
            RoadSlice::MARKING_TEX_COORDS_SIZE_IN_FLOAT) * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(4);
    
    glVertexAttribPointer(
        5, 
        RoadSlice::BEZIER_HANDLES_SIZE_IN_FLOAT, 
        GL_FLOAT, 
        GL_FALSE, 
        VERTEX_SIZE_IN_FLOAT * sizeof(GLfloat), 
        (GLvoid*)((RoadSlice::POSITION_SIZE_IN_FLOAT + 
            RoadSlice::NORMAL_SIZE_IN_FLOAT + 
            RoadSlice::ASPHALT_TEX_COORDS_SIZE_IN_FLOAT +
            RoadSlice::MARKING_TEX_COORDS_SIZE_IN_FLOAT +
            RoadSlice::BEZIER_ANCHORS_SIZE_IN_FLOAT) * sizeof(GLfloat))
    );
    glEnableVertexAttribArray(5);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    delete [] generalBuffer;

    debug("Done!");
}

void Landscape::RenderRoad(RENDER_MODE mode)
{
    m_road->Render(mode);
}

void Landscape::RenderBuildings(RENDER_MODE mode)
{
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (Building* b: m_buildings)
        b->Render(mode);
    // glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

bool Landscape::PopulationGenerator::IsBiomePopulatable(Biome b)
{
    return b == Biome::PLAIN;
}

void Landscape::GeneratePopulationMap()
{
    debug("Generating population map");

    if (m_biomemap == nullptr)
    {
        debug("NO BIOME MAP!");
        throw std::exception{};
    }
    m_populationmap = new PopulationMap(*m_biomemap);

    float LOWLAND_R = PopulationGenerator::POPULATION_RATIO_LOWLAND;
    float CENTER_R = PopulationGenerator::POPULATION_RATIO_DISTANCE_TO_CENTER;
    float SHORE_R = PopulationGenerator::POPULATION_RATIO_DISTANCE_TO_SHORE;
    float BIOME_R = PopulationGenerator::POPULATION_RATIO_BIOME_POSSIBILITY;
    float PERLIN_R = PopulationGenerator::POPULATION_RATIO_PERLIN;
    float RATIO_SUM = PopulationGenerator::POPULATION_RATIO_POSITIVE_SUM;

    float maxSurfIntensity = m_surfacemap->GetMaxIntensity(Surface::EARTH);

    noise::module::Perlin populationPerlin;
    populationPerlin.SetFrequency(PopulationGenerator::POPULATION_PERLIN_FREQUENCY_REAL);
    populationPerlin.SetOctaveCount(PopulationGenerator::POPULATION_PERLIN_OCTAVE_COUNT);
    populationPerlin.SetPersistence(PopulationGenerator::POPULATION_PERLIN_PERSISTANCE);
    populationPerlin.SetLacunarity(PopulationGenerator::POPULATION_PERLIN_LACUNARITY);

    for (int i = -m_populationmap->m_width; i <= m_populationmap->m_width; i++)
    {
        for (int j = -m_populationmap->m_height; j <= m_populationmap->m_height; j++)
        {
            PopulationPoint result{};

            vec2Int pos{i, j};
            glm::vec2 curRealCoords = PointMapIntToReal(m_biomemap, pos);            

            if (!PopulationGenerator::IsBiomePopulatable(m_biomemap->Get(pos)))
            {
                result.population = Population::UNUSABLE;
                m_populationmap->Set(pos, result);
                continue;
            }

            result.population = Population::POPULATED;
            result.intensity = 0;

            float waterQuality;
            SurfacePoint sp = m_surfacemap->GetSurface(pos);
            waterQuality = (sp.surface == Surface::EARTH) ? 1 - (sp.intensity / maxSurfIntensity) : 0.f;

            float lowlandQuality;
            const int LOWLAND_ITERATIONS = 2;
            const float STEP_REAL_SIZE = 1;
            float currentHeight = -10000.f;
            float minHeight = 100000;
            float maxHeight = -100000;
            bool closeToUnfitLand = false;
            for (int k = -LOWLAND_ITERATIONS; k <= LOWLAND_ITERATIONS; k++)
            {
                for (int l = -LOWLAND_ITERATIONS; l <= LOWLAND_ITERATIONS; l++)
                {
                    glm::vec2 temp = curRealCoords + glm::vec2{k,l} * STEP_REAL_SIZE;

                    glm::vec2 check = PointRealToMapNormalized(m_biomemap, temp, true);
                    if(!m_biomemap->IsPointValid_IntContinuous(m_biomemap->PointNormalizedToIntContinuous(check, true)))
                    {
                        closeToUnfitLand = true;
                        continue;
                    }
                    if (!PopulationGenerator::IsBiomePopulatable(m_biomemap->Get(m_biomemap->PointNormalizedToInt(check))))
                    {
                        closeToUnfitLand = true;
                        continue;
                    }
                    
                    temp = PointRealToMapNormalized(m_heightmap, temp);
                    temp = m_heightmap->PointNormalizedToIntContinuous(temp, true);
                    if(!m_heightmap->IsPointValid(vec2Int(temp)))
                    {
                        closeToUnfitLand = true;
                        continue;
                    }

                    float height = m_heightmap->Get(vec2Int(temp));
                    minHeight = std::min(height, minHeight);
                    maxHeight = std::max(height, maxHeight);
                    if (k == 0 && l == 0)
                        currentHeight = height;
                }
            } 
            if (currentHeight < -1000)
            {
                debug("CURRENT POINT IS NOT AVAILABLE!", i, j);
                throw std::exception{};
            }
            lowlandQuality = 1 - (currentHeight - minHeight) / (maxHeight - minHeight);

            float biomeQuality = 0.f;
            const int BIOMES_ITERATIONS = 3;
            for (int k = -BIOMES_ITERATIONS; k <= BIOMES_ITERATIONS; k++)
            {
                for (int l = -BIOMES_ITERATIONS; l <= BIOMES_ITERATIONS; l++)
                {
                    glm::vec2 posReal = curRealCoords + glm::vec2{k,l} * STEP_REAL_SIZE;
                    glm::vec2 check = PointRealToMapNormalized(m_biomemap, posReal, true);
                    if(!m_biomemap->IsPointValid_IntContinuous(m_biomemap->PointNormalizedToIntContinuous(check, true)))
                    {
                        biomeQuality = 1.f;
                        continue;
                    }
                    if (!PopulationGenerator::IsBiomePopulatable(m_biomemap->Get(m_biomemap->PointNormalizedToInt(check))))
                    {
                        biomeQuality = 1.f;
                        continue;
                    }
                    
                    glm::vec2 heightPos = PointRealToMapNormalized(m_heightmap, posReal);
                    heightPos = m_heightmap->PointNormalizedToIntContinuous(heightPos, true);
                    if(!m_heightmap->IsPointValid(vec2Int(heightPos)))
                    {
                        biomeQuality = 1.f;
                        continue;
                    }

                    glm::vec2 waterPos = PointRealToMapNormalized(m_surfacemap, posReal);
                    vec2Int waterPosInt = m_surfacemap->PointNormalizedToInt(waterPos);
                    SurfacePoint sp = m_surfacemap->GetSurface(waterPosInt);
                    float shoreRatio = 
                        sp.intensity * GetMapCellRealSize(m_surfacemap) / 
                            (WaterGenerator::SHORE_REAL_SIZE + RoadModelGenerator::ROAD_FULL_REAL_WIDTH / 2);
                    if(shoreRatio < 1.01f)
                    {
                        biomeQuality = 1.f;
                        continue;
                    }
                }
            } 

            float perlinQuality;
            perlinQuality = (float) populationPerlin.GetValue(i, 0, j);
            perlinQuality = perlinQuality * 0.5 + 0.5;

            

            result.intensity = 
                biomeQuality * BIOME_R + 
                lowlandQuality * LOWLAND_R + 
                waterQuality * SHORE_R + 
                perlinQuality * PERLIN_R;

            result.population = Population::EMPTY;

            // result.population = Population::POPULATED;
            // result.intensity = perlinQuality;

            m_populationmap->Set(pos, result);
        } 
    }

    glm::vec2 cityCenterPosNormalized{0.5f, 0.5};
    vec2Int cityCenterPos = m_biomemap->PointNormalizedToInt(cityCenterPosNormalized);
    m_populationmap->realCenterPos = cityCenterPos;
    vec2Int newCenter;
    float minWaterDistance = 10000;
    float minCenterDistance = 10000;
    for (int i = -m_populationmap->m_width; i <= m_populationmap->m_width; i++)
    {
        for (int j = -m_populationmap->m_height; j <= m_populationmap->m_height; j++)
        {
            vec2Int pos{i, j};
            PopulationPoint pp = m_populationmap->Get(pos);
            if (pp.population != Population::EMPTY || pp.intensity < 0)
                continue;

            float waterDist = m_surfacemap->GetSurface(
                    m_surfacemap->PointNormalizedToInt(m_populationmap->PointIntToNormalized(pos))
                ).intensity;
            float centerDist = (pos - cityCenterPos).length();
            if (waterDist < minWaterDistance)
            {
                minWaterDistance = waterDist;
                minCenterDistance = centerDist;
                newCenter = pos;
                continue;
            }
            if (Approx(waterDist, minWaterDistance) && centerDist < minCenterDistance)
            {
                minWaterDistance = waterDist;
                minCenterDistance = centerDist;
                newCenter = pos;
            }
        }
    }
    cityCenterPos = newCenter;
    m_populationmap->shoreCenterPos = cityCenterPos;

    if (m_populationmap->Get(cityCenterPos).population != Population::EMPTY)
    {
        debug("CITY CENTER MUST BE AVAILABLE!");
        throw std::exception{};
    }

    m_populationmap->Set(cityCenterPos, PopulationPoint(Population::POPULATED, 1.f));
    const float maxPopulation = PopulationGenerator::POPULATION_VOLUME_REAL / GetMapCellRealSize(m_biomemap);
    float populationLeft = maxPopulation;

    std::function populationFindNeighbs_FindEmpty
    {
        [this]
        (vec2Int coords) 
        {
            return m_populationmap->populationFindNeighbs_FindPopulation(coords, Population::EMPTY);
        }
    };
    std::function populationApplyToReachedPoints_SetPopulated
    {
        [this, &populationLeft, &maxPopulation]
        (vec2Int coords, float distance)
        {
            float population = PopulationGenerator::GetPointPopulationByLeftPart(populationLeft / maxPopulation);
            populationLeft -= population;
            if (populationLeft <= 0.f)
            {
                populationLeft = 0.f;
                m_populationmap->Set(coords, PopulationPoint(Population::EMPTY, 0));
            }
            else
            {
                m_populationmap->Set(coords, PopulationPoint(Population::POPULATED, population));
            }
        }
    };
    std::function populationFindDistance
    {
        [this, &populationLeft, &maxPopulation, CENTER_R, RATIO_SUM]
        (vec2Int current, vec2Int previous, float prevDist, vec2Int origin)
        {
            float resultWeight = 0.f;

            float populationLeftPart = populationLeft / maxPopulation;
            float centerQuality = PopulationGenerator::GetPointCenterDistancePartByLeftPart(populationLeftPart);
            float currentQuality = m_populationmap->Get(current).intensity;
            currentQuality += centerQuality * CENTER_R;
            if (currentQuality <= 0)
                resultWeight = 100000000.f;

            currentQuality = currentQuality / RATIO_SUM;
            currentQuality = 1.f - currentQuality;

            return prevDist + currentQuality;
        }
    };

    biomeBreadthSearch(
        std::vector<vec2Int> {cityCenterPos},
        populationFindNeighbs_FindEmpty,
        populationApplyToReachedPoints_SetPopulated,
        populationFindDistance
    );

    m_populationmap->SaveAsTexture();

    debug("Done!");


    debug("Creating CITY biome");
    
    if (m_populationmap->m_width != m_biomemap->m_width || m_populationmap->m_height != m_biomemap->m_height)
    {
        debug("POPULATION AND BIOME MAP DO NOT MATCH");
        throw std::exception{};
    }

    for (int i = -m_biomemap->m_width; i <= m_biomemap->m_width; i++)
    {
        for (int j = -m_biomemap->m_height; j <= m_biomemap->m_height; j++)
        {
            vec2Int coords{i, j};
            PopulationPoint pp = m_populationmap->Get(coords);
            if (pp.population == Population::POPULATED)
            {
                m_biomemap->Set(coords, BiomePoint{Biome::CITY, 0.f});
            }
        }
    }
    
    m_biomemap->CalibrateBiomemap();

    debug("Done!");
}

float Landscape::PopulationGenerator::GetPointPopulationByLeftPart(float populationLeft)
{
    if (populationLeft <= 0.1)
        return 0.1;
    if (populationLeft >= 0.9)
        return 0.9f;
    return populationLeft;
}

float Landscape::PopulationGenerator::GetPointCenterDistancePartByLeftPart(float populationLeft)
{
    return populationLeft;
}

float Landscape::CalculateWaterLevel()
{
    float level = 100000.f;
    for (int i = -m_heightmap->m_width; i <= m_heightmap->m_width; i++)
    {
        for (int j = -m_heightmap->m_height; j <= m_heightmap->m_height; j++)
        {
            vec2Int heightmapCoords{i, j};
            vec2Int biomemapCoords = m_biomemap->PointNormalizedToInt(
                    m_heightmap->PointIntToNormalized(heightmapCoords)
                );
            if (m_biomemap->Get(biomemapCoords).biome != Biome::PLAIN)
                continue;
            level = std::min(level, m_heightmap->Get(heightmapCoords) - WaterGenerator::GROUND_HEIGHT_ABOVE_WATER);
        }
    }
    return level;
}

float Landscape::GetPopulationRatioByRealPos(glm::vec2 p)
{
    vec2Int coords = m_populationmap->PointNormalizedToInt(PointRealToMapNormalized(m_heightmap, p));
    PopulationPoint pp = m_populationmap->Get(coords);
    if (pp.population != Population::POPULATED)
        return 0.f;
    return pp.intensity;

}

bool Landscape::IsCrossroadJoinable(RoadNode* cross)
{
    return cross->edges.size() < RoadSchemeGenerator::MAX_CROSSROAD_ROADS;
}

bool Landscape::IsCrossroadJoinable(RoadNode* cross, glm::vec2 posNew)
{
    if (!IsCrossroadJoinable(cross))
        return false;
    
    float minAngleRadians = M_PI + 0.1f; //min angle to other roads in this crossroad
    for (RoadSection* s: cross->edges)
    {
        glm::vec2 dir = glm::normalize(posNew - cross->pos);
        float angle = glm::dot(s->StraightDir(cross), dir);
        angle = SafeAcosf(angle);
        minAngleRadians = std::min(minAngleRadians, angle);
    }
    float borderRadians = RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES();
    borderRadians = glm::radians(borderRadians);
    return borderRadians < minAngleRadians;
}

bool Landscape::IsSectionJoinable(RoadSection* section, float ratio, glm::vec2 posNew)
{
    glm::vec2 midRoadPoint = LERP(section->start->pos, section->end->pos, ratio);
    float minAngleRadians = M_PI + 0.1f; //min angle to other roads in this crossroad
    for (RoadNode* n: section->GetEdges())
    {
        glm::vec2 dir = glm::normalize(posNew - midRoadPoint);
        float angle = glm::dot(glm::normalize(n->pos - midRoadPoint), dir);
        angle = SafeAcosf(angle);
        minAngleRadians = std::min(minAngleRadians, angle);
    }
    float borderRadians = RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES();
    borderRadians = glm::radians(borderRadians);
    return borderRadians < minAngleRadians;
}

void Landscape::RoundOffTurns()
{
    for (RoadSection* section: m_roadScheme->graph.GetSections())
    {
        if (section->start->IsTurn() || section->end->IsTurn())
        {
            if (section->start->IsEndNode() || section->end->IsEndNode())
                continue;

            RoadSection* other = nullptr;
            if (section->start->IsTurn())
            {
                other = (section->start->edges[0] == section)
                    ? section->start->edges[1]
                    : section->start->edges[0];
            }
            else
            {
                other = (section->end->edges[0] == section)
                    ? section->end->edges[1]
                    : section->end->edges[0];
            }

            if (RoadSection::AreParrallel(section, other))
                continue;

            RoadNode* newnode = m_roadScheme->graph.DivideSection(section->start, section->end, 0.5f);
        }
    }

    if (m_realSize[0] != m_realSize[1])
    {
        debug("NORMALIZED COORDS CAN NOT BE USED IN THIS CASE");
        throw std::exception{};
    }

    const float minRadiusNormalized = RoadModelGenerator::GET_TURN_MIN_RADIUS_RATIO()
        *  RoadSlice::GetDefaultPlanarSlicewidth() * 0.5f / m_realSize[0];
    const float preferRadiusNormalized = RoadModelGenerator::TURN_PREFER_RADIUS_RATIO
        *  RoadSlice::GetDefaultPlanarSlicewidth() * 0.5f / m_realSize[0];

    
    // const float TURN_RADIUS_TO_ROAD_WIDTH = 1.5f;
    // float ratio = std::max(
    //     TURN_RADIUS_TO_ROAD_WIDTH, 
    //     RoadSlice::GetDefaultFullSlicewidth() / RoadSlice::GetDefaultPlanarSlicewidth()
    // );
    // ratio = TURN_RADIUS_TO_ROAD_WIDTH;
    
    // float radiusNormalized = ratio * RoadSlice::GetDefaultPlanarSlicewidth() * 0.5f / m_realSize[0];
    bool flag = true;
    while (flag)
    {
        flag = false;
        for (RoadNode* t: m_roadScheme->graph.GetNodes())
        {
            //REWORK TWO ROUND TURNS
            if (!t->IsTurn())
                continue;

            if (t->edges[0]->IsStraight() && t->edges[1]->IsStraight())
            {
                if (Approx(fabsf(glm::dot(t->edges[0]->StraightDir(), t->edges[1]->StraightDir())), 1.f))
                    continue;
            }
            else if (t->edges[0]->IsRoundTurn() && t->edges[1]->IsRoundTurn())
            {
                continue;
            }
            else
            {
                RoadSection* round = (t->edges[0]->IsRoundTurn())
                    ? t->edges[0]
                    : t->edges[1];
                RoadSection* straight = (t->edges[0]->IsRoundTurn())
                    ? t->edges[1]
                    : t->edges[0];
                glm::vec2 dir = (round->start == t)
                    ? round->startDir
                    : round->endDir;
                if (Approx(fabsf(glm::dot(glm::normalize(dir), straight->StraightDir())), 1.f))
                    continue;
                else
                {
                    debug("STRAIGHT AND ROUND MUST ALLIGN!");
                    throw std::exception{};
                }
            }

            RoadNode* turn = t;
            std::array<RoadNode*, 2> roads;
            std::array<glm::vec2, 2> dirs;
            for (int i = 0; i < 2; i++)
            {
                RoadNode* other = (turn->edges[i]->start == t) ? turn->edges[i]->end : turn->edges[i]->start;
                roads[i] = other;
                dirs[i] = glm::normalize(roads[i]->pos - turn->pos);
            }
            float angleRadians = SafeAcosf(glm::dot(dirs[0], dirs[1]));

            float currentRoadLenNormalized = std::min(
                    glm::length(roads[0]->pos - turn->pos), 
                    glm::length(roads[1]->pos - turn->pos)
                );
            float curRadiusNormalized = glm::tan(angleRadians / 2) * currentRoadLenNormalized;

            bool generateNewNodes = true;
            if (Approx(curRadiusNormalized, minRadiusNormalized))
            {
                curRadiusNormalized = minRadiusNormalized;
                generateNewNodes = false;
            }
            if (curRadiusNormalized < minRadiusNormalized)
            {
                debug("ROAD IS TOO SHORT FOR THIS TURN!");
                throw std::exception{};
            }
            else
            {
                curRadiusNormalized = clamp(curRadiusNormalized, minRadiusNormalized, preferRadiusNormalized);
            }
            
            if (generateNewNodes)
            {
                std::array<glm::vec2, 2> newNodesPoses;
                std::array<RoadNode*, 2> newNodes;
                float dist = curRadiusNormalized / tanf(angleRadians / 2);
                for (int i = 0; i < 2; i++)
                {
                    newNodesPoses[i] = turn->pos + dirs[i] * dist;
                    newNodes[i] = m_roadScheme->graph.AddNode(newNodesPoses[i]);
                }

                for (int i = 0; i < 2; i++)
                {
                    m_roadScheme->graph.AddSection(roads[i], newNodes[i], 0.f);
                }
                RoadSection* s = m_roadScheme->graph.AddSection(newNodes[0], newNodes[1], -dirs[0], dirs[1]);
                m_roadScheme->graph.RemoveNode(t);
            }
            else
            {
                RoadSection* s = m_roadScheme->graph.AddSection(roads[0], roads[1], -dirs[0], dirs[1]);
                m_roadScheme->graph.RemoveNode(t);
            }

            flag = true;
            break;
        }
    }
}

float Landscape::RoadModelGenerator::GET_TURN_MIN_RADIUS_RATIO()
{
    return RoadSlice::GetDefaultFullSlicewidth() / RoadSlice::GetDefaultPlanarSlicewidth();
}

float Landscape::RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES()
{
    float minDist = RoadSchemeGenerator::MIN_DISTANCE_BETWEEN_CROSSROADS / 2;
    float minRadius =  GET_TURN_MIN_RADIUS_RATIO() * RoadSlice::GetDefaultPlanarSlicewidth() * 0.5f;
    float angle = glm::degrees(glm::atan(minRadius / minDist) * 2.f);
    return angle;
}

std::array<std::vector<RoadSection*>, 2> Landscape::GetLeftRightTurns(RoadNode* crossroad, RoadSection& section)
{
    if (!crossroad->IsCrossroad())
    {
        debug("IT IS NOT CROSSROAD!");
        throw std::exception{};
    }
    glm::vec2 leftDir = LeftNormal(section.StraightDir(crossroad));
    glm::vec2 rightDir = -leftDir;
    std::vector<RoadSection*> leftTurns, rightTurns;
    for (RoadSection* neighbSection : crossroad->edges)
    {   
        if (*neighbSection == section)
            continue;
        if (!neighbSection->IsStraight())
        {
            debug("NEIGHBS SHOULD BE STRAIGHT!");
            throw std::exception{};
        }

        if (glm::dot(neighbSection->StraightDir(crossroad), leftDir) > 0)
            leftTurns.push_back(neighbSection);
        else
            rightTurns.push_back(neighbSection);
    }
    std::sort(leftTurns.begin(), leftTurns.end(),
        [&section, &crossroad] (RoadSection* a, RoadSection* b) 
        { 
            float cosA = glm::dot(a->StraightDir(crossroad), section.StraightDir(crossroad));
            float cosB = glm::dot(b->StraightDir(crossroad), section.StraightDir(crossroad));
            return cosA > cosB;
        }
    );
    std::sort(rightTurns.begin(), rightTurns.end(),
        [&section, &crossroad] (RoadSection* a, RoadSection* b) 
        { 
            float cosA = glm::dot(a->StraightDir(crossroad), section.StraightDir(crossroad));
            float cosB = glm::dot(b->StraightDir(crossroad), section.StraightDir(crossroad));
            return cosA > cosB;
        }
    );

    std::array<std::vector<RoadSection*>, 2> turns;
    turns[0] = leftTurns;
    turns[1] = rightTurns;

    return turns;
}

//roadDirection is facing right from foundation
RoadSlice Landscape::CreateSliceFromSchemePoint (glm::vec2 roadPosition, glm::vec2 roadDirection, float scale)
{
    float angleRadians  = atan2f(roadDirection.y, roadDirection.x) + M_PI * 0.5;
    glm::vec2 crossDir{cosf(angleRadians), sinf(angleRadians)};
    float sliceLen = RoadSlice::GetDefaultFullSlicewidth() * scale;
    glm::vec2 centerRealXZ = PointMapNormalizedToReal(m_heightmap, roadPosition);
    
    std::array<glm::vec3, 2> edgeRealPositions;
    for (int i = 0; i < 2; i++)
    {
        glm::vec2 edgeRealXZ = centerRealXZ;
        edgeRealXZ += crossDir * sliceLen * 0.5f * (float)(i * 2 - 1);
        glm::vec2 currentCoordsNormalized = PointRealToMapNormalized(m_heightmap, edgeRealXZ);
        if (!m_heightmap->IsPointValid_Normalized(currentCoordsNormalized))
        {
            debug("ROAD IS GOING OUT OF MAP!");
            throw std::exception{};
        }
        edgeRealPositions[i] = glm::vec3{
            edgeRealXZ.x,
            m_heightmap->GetInterpolated_IntContinuous(
                m_heightmap->PointNormalizedToIntContinuous(
                    currentCoordsNormalized
                )
            ) + RoadModelGenerator::ROAD_HEIGHT_ABOVE_HEIGHTMAP,
            edgeRealXZ.y
        };
    }
    return RoadSlice((edgeRealPositions[0] + edgeRealPositions[1]) * 0.5f,
        edgeRealPositions[1] - edgeRealPositions[0], scale);
}

//DO NOT USE IT! IT IS ONLY FOR INNER PURPOSES!
std::array<glm::vec2, 2> Landscape::GetPlanarPointsForCrossroadSection(RoadNode* crossroad, RoadSection& section)
{
    if (!crossroad->IsCrossroad())
    {
        debug("IT IS NOT CROSSROAD!");
        throw std::exception{};
    }
    glm::vec2 leftDir = LeftNormal(section.StraightDir(crossroad));
    glm::vec2 rightDir = -leftDir;
    auto turns = GetLeftRightTurns(crossroad, section);
    std::vector<RoadSection*> leftTurns = turns[0];
    std::vector<RoadSection*> rightTurns = turns[1];
    

    RoadSection* leftNeighb = (leftTurns.size() > 0) ? 
        leftTurns[0] : rightTurns[rightTurns.size() - 1];
    RoadSection* rightNeighb = (rightTurns.size() > 0) ? 
        rightTurns[0] : leftTurns[leftTurns.size() - 1];
    
    std::array<glm::vec2, 2> planarRealPoints;
    std::array<RoadSection*, 2> LRNeighbs = {leftNeighb, rightNeighb};
    int distantNeighbInd = -1; // -1 for none, 0 for left, 1  for right
    distantNeighbInd = (leftTurns.size() == 0) ? 0 : distantNeighbInd;
    distantNeighbInd = (rightTurns.size() == 0) ? 1 : distantNeighbInd;

    for (int i = 0; i < 2; i++)
    {
        float normalLen = RoadSlice::GetDefaultPlanarSlicewidth() / 2;
        RoadSection* curNeighb = LRNeighbs[i];
        float alphaRadians = glm::dot(
            curNeighb->StraightDir(crossroad), section.StraightDir(crossroad));
        alphaRadians = SafeAcosf(alphaRadians);
        
        float border = RoadModelGenerator::GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES();
        if (glm::degrees(alphaRadians) < border)
        {
            debug("ROADS ARE TO CLOSE IN CROSSROAD! angle is ", glm::degrees(alphaRadians), "degrees");
            debug("CROSSROAD", crossroad->pos);
            for (RoadSection* s: crossroad->edges)
            {
                RoadNode* n = (s->start == crossroad)
                    ? s->end
                    : s->start;
                debug("ROAD", n->pos);
            }
            throw std::exception{};
        }

        alphaRadians /= 2;

        if (distantNeighbInd == i)
        {
            alphaRadians = M_PI_2 - alphaRadians;
            float subSectionLen = normalLen * tanf(alphaRadians);
            planarRealPoints[i] = subSectionLen * (-section.StraightDir(crossroad)) + 
                PointMapNormalizedToReal(m_heightmap, crossroad->pos); 
        }
        else
        {
            float subSectionLen = normalLen / tanf(alphaRadians);
            planarRealPoints[i] = subSectionLen * section.StraightDir(crossroad) + 
                PointMapNormalizedToReal(m_heightmap, crossroad->pos); 
        }
        
        planarRealPoints[i] += ((i == 0) ? leftDir : rightDir) * normalLen; 
    }

    //Checking if straightDir is not matching right of planarRealPoints
    glm::vec2 l = LeftNormal(
        PointRealToMapNormalized(m_heightmap,planarRealPoints[1]) - 
        PointRealToMapNormalized(m_heightmap,planarRealPoints[0]));
    if (glm::dot(l, section.StraightDir()) > 0)
    {
        glm::vec2 tmp = planarRealPoints[0];
        planarRealPoints[0] = planarRealPoints[1];
        planarRealPoints[1] = tmp;
    }
    
    return planarRealPoints;
}

//roadDirection is facing right from foundation
RoadSlice Landscape::CreateSliceFromCrossroad_Uncorrected(RoadSection& section, bool isEndNode)
{
    RoadNode* crossroad = (isEndNode) ? section.end : section.start;
    std::array<glm::vec2, 2> planarRealPoints = GetPlanarPointsForCrossroadSection(crossroad, section);
    
    //Generating Slice
    glm::vec2 slicePos = PointRealToMapNormalized(
        m_heightmap,(planarRealPoints[0] + planarRealPoints[1]) * 0.5f);
    glm::vec2 foundation = glm::normalize(
        PointRealToMapNormalized(m_heightmap,planarRealPoints[1]) - 
        PointRealToMapNormalized(m_heightmap,planarRealPoints[0]));
    glm::vec2 direction = -LeftNormal(foundation);
    
    float scale = glm::length(planarRealPoints[1] - planarRealPoints[0]) / RoadSlice::GetDefaultPlanarSlicewidth();
    RoadSlice result = CreateSliceFromSchemePoint(slicePos, direction, scale);
    auto planarInds = RoadSlice::GetPlanarBorderIndexes();
    glm::vec3 sectionCenter = (result.border[planarInds[0]] + result.border[planarInds[1]]) * 0.5f;
    //Correcting edges cross pane
    for (int i = 0; i < planarInds.size(); i++)
    {
        result.border[planarInds[i]] = X0Z(planarRealPoints[i]);
        result.border[planarInds[i]].y = 
            m_heightmap->GetInterpolated_IntContinuous(
                m_heightmap->PointNormalizedToIntContinuous(
                    PointRealToMapNormalized(m_heightmap, planarRealPoints[i])
                )
            ) + RoadModelGenerator::ROAD_HEIGHT_ABOVE_HEIGHTMAP;
    }
    
    //result
    return result;
}

std::vector<Basement> Landscape::GenerateBasements()
{
    debug("Generating basements!");

    std::vector<Basement> basements;
    const float preferRadiusNormalized = RoadModelGenerator::TURN_PREFER_RADIUS_RATIO
        *  RoadSlice::GetDefaultPlanarSlicewidth() * 0.5f / m_realSize[0];
    const float distRoadToBasementNorm = BuildingsGenerator::DIST_BETWEEN_ROAD_AND_BASEMENT_REAL / m_realSize[0];
    // debug("1");
    std::function<RoadNode* (std::vector<RoadNode*>& ,int)> GetNode
    {
        [] (std::vector<RoadNode*>& path, int index)
        {
            if (index == -1)
                return path[path.size() - 1];
            else if (index >= 0 && index < path.size())
                return path[index];
            else if (index == path.size())
                return path[0];
            else
            {
                debug ("WRONG INDEX!", index, path.size());
                throw std::exception{};
            }
        }
    };

    auto roadSections = m_roadScheme->graph.GetSections();
    std::sort(
        roadSections.begin(),
        roadSections.end(),
        [] (RoadSection* a, RoadSection* b) 
        { 
            return (a->start->pos.x + 10*a->start->pos.y + 100*a->end->pos.x + 1000*a->end->pos.y ) > 
                (b->start->pos.x + 10*b->start->pos.y + 100*b->end->pos.x + 1000*b->end->pos.y );
        }
    );
    std::unordered_set<RoadSection*> reached{};
    // debug("2");
    std::vector<Polygon> roadBlockBasements; 
    const float SHIFT = preferRadiusNormalized + distRoadToBasementNorm - 0.00001;
    for (RoadSection* section: roadSections)
    {
        Polygon cur;
        if (section->IsStraight())
        {
            float roadWidthNorm = RoadSlice::GetDefaultPlanarSlicewidth() * 0.5 / m_realSize[0];
            glm::vec2 left = LeftNormal(section->StraightDir());
            for (int i = 0; i < 4; i++)
            {
                bool isLeft = i < 2;
                bool isEnd = (i % 3 != 0);

                if ((isEnd && section->end->IsCrossroad()) || (!isEnd && section->start->IsCrossroad()))
                {
                    RoadNode* curNode = (isEnd) ? section->end : section->start;
                    auto planarPoints = GetPlanarPointsForCrossroadSection(curNode, *section);
                    for (int j = 0; j < 2; j++)
                    {
                        glm::vec2 curPlanarPoint = PointRealToMapNormalized(m_heightmap, planarPoints[j]);
                        bool planarPointIsOnLeft = glm::dot(curPlanarPoint - curNode->pos, left) > 0;
                        if ((isLeft && planarPointIsOnLeft) || (!planarPointIsOnLeft && !isLeft))
                        {
                            glm::vec2 point = PointRealToMapNormalized(m_heightmap, planarPoints[j]);
                            point = curNode->pos + 
                                glm::normalize(point - curNode->pos) * 
                                    glm::length(point - curNode->pos) / roadWidthNorm *
                                    (preferRadiusNormalized + distRoadToBasementNorm - 0.00001f);
                            cur.points.push_back(PointMapNormalizedToReal(m_heightmap, point));
                            break;
                        }
                    }
                }
                else
                {
                    glm::vec2 point = isEnd ? section->end->pos : section->start->pos;
                    point += LeftNormal(section->StraightDir()) * (isLeft ? 1.f : -1.f) * SHIFT;
                    cur.points.push_back(PointMapNormalizedToReal(m_heightmap, point));
                }
            }
        }
        else if (section->IsRoundTurn())
        {
            glm::vec2 center = section->RoundTurnCenter();
            float angleRadians = SafeAcosf(glm::dot(glm::normalize(section->startDir), glm::normalize(-section->endDir)));
            angleRadians /= 2;
            float x = glm::length(section->start->pos - section->end->pos) * 0.5f;
            float l = x / sinf(angleRadians) / cosf(angleRadians);
            l *= (2 * preferRadiusNormalized + distRoadToBasementNorm - 0.00001) / (preferRadiusNormalized);
            glm::vec2 hordeCenter = (section->start->pos + section->end->pos) * 0.5f;
            glm::vec2 antiCenter = center + glm::normalize(hordeCenter - center) * l;

            for (int i = 0; i < 4; i++)
            {
                glm::vec2 p;
                glm::vec2 dir;
                float _len = 2 * preferRadiusNormalized + distRoadToBasementNorm - 0.00001;
                switch (i)
                {
                    case 0:
                        p = center;
                        break;
                    case 1:
                        dir = glm::normalize(section->start->pos - center);
                        p = center + dir * _len;
                        break;
                    case 2:
                        p = antiCenter;
                        break;
                    case 3:
                        dir = glm::normalize(section->end->pos - center);
                        p = center + dir * _len;
                        break;
                };
                cur.points.push_back(PointMapNormalizedToReal(m_heightmap, p));
            }
        }
        else
        {
            debug("STRANGE TURN");
            throw std::exception{};
        }
        roadBlockBasements.push_back(cur);
    }

    std::vector<Basement> innerBasements;
    while (true)
    {
        RoadSection* startSec = nullptr;
        for (RoadSection* rs: roadSections)
        {
            if (reached.find(rs) != reached.end())
                continue;
            startSec = rs;
        } 

        if (startSec == nullptr)
            break;

        for (int side = 0; side < 2; side++)
        {
            bool leftSide = (side == 0);

            bool result = false;
            std::vector<RoadNode*> path{startSec->start, startSec->end};        
            while(true)
            {
                RoadNode* curEnd = path[path.size() - 1];
                RoadSection* cur = RoadSection::SectionByEdges(path[path.size() - 2], path[path.size() - 1]);

                if (curEnd->IsEndNode())
                {
                    result = false;
                    reached.emplace(cur);
                    break;
                }
                else if (curEnd->IsTurn() || curEnd->IsCrossroad())
                {
                    RoadSection* next;
                    if (curEnd->IsTurn())
                    {
                        next = (curEnd->edges[0] == cur)
                            ? curEnd->edges[1]
                            : curEnd->edges[0];
                    }
                    else
                    {
                        auto turns = GetLeftRightTurns(curEnd, *cur);
                        if (leftSide)
                        {
                            next = (turns[1].size() != 0)
                                ? turns[1][0]
                                : turns[0][turns[0].size() - 1];
                        }
                        else
                        {
                            next = (turns[0].size() != 0)
                                ? turns[0][0]
                                : turns[1][turns[1].size() - 1];
                        }
                    }

                    if (reached.find(next) != reached.end())
                    {
                        result = false;
                        break;
                    }

                    if (next->Other(curEnd) == path[0])
                    {
                        result = true;
                        break;
                    }
                    else
                    {
                        path.push_back(next->Other(curEnd));
                    }
                }
                else
                {
                    debug("WHAT ARE YOU???");
                    throw std::exception{};
                }
            }

            if (result == true)
            {
                Polygon protoBasement;
                
                for (int i = 0; i < path.size(); i++)
                {
                    RoadNode* n = path[i];
                    glm::vec3 posReal;

                    if (!n->IsCrossroad() && !n->IsTurn())
                    {
                        debug("WHAT THE FUCK IS IT?");
                        throw std::exception{};
                    }

                    std::array<RoadNode*, 2> edges{GetNode(path, i - 1), GetNode(path, i + 1)};
                    std::array<RoadSection*, 2> sections{
                        RoadSection::SectionByEdges(n, edges[0]), 
                        RoadSection::SectionByEdges(n, edges[1])
                    };
                    // debug("2_5");
                    if (RoadSection::AreParrallel(sections[0], sections[1]))
                    {
                        // debug("2_6");
                        RoadNode* next = GetNode(path, i + 1);
                        RoadSection* sec = RoadSection::SectionByEdges(n, next);
                        glm::vec2 dir = sec->StraightDir(n);
                        dir = LeftNormal(dir);
                        dir = (leftSide) ? dir : -dir;
                        dir *= preferRadiusNormalized + distRoadToBasementNorm;
                        posReal = X0Z(PointMapNormalizedToReal(m_heightmap, n->pos + dir));
                        // debug("2_7");
                    }
                    else
                    {
                        // debug("2_8");
                        if (n->IsCrossroad())
                        {
                            auto planarPoints1 = GetPlanarPointsForCrossroadSection(n, *sections[0]);
                            auto planarPoints2 = GetPlanarPointsForCrossroadSection(n, *sections[1]);

                            bool success = false;
                            glm::vec2 point{-1, -1};
                            for (int l = 0; l < 2; l++)
                            {
                                for (int p = 0; p < 2; p++)
                                {
                                    if (Approx(planarPoints1[l], planarPoints2[p], 0.01))
                                    {
                                        point = PointRealToMapNormalized(m_heightmap, planarPoints1[l]);
                                        glm::vec2 dir = glm::normalize(point - n->pos);
                                        float dist = glm::length(point - n->pos) * RoadModelGenerator::TURN_PREFER_RADIUS_RATIO;
                                        dist += distRoadToBasementNorm * (glm::length(point - n->pos) / preferRadiusNormalized);
                                        point = n->pos + dir * dist;
                                        success = true;
                                    }
                                }
                            }

                            if (!success)
                            {
                                debug("CORRECT POINT NOT FOUND !");
                                throw std::exception{};
                            }

                            posReal = X0Z(PointMapNormalizedToReal(m_heightmap, point));
                        }
                        // debug("2_9");
                        if (n->IsTurn())
                        {
                            RoadNode* next = GetNode(path, i + 1);
                            RoadSection* sec = RoadSection::SectionByEdges(n, next);

                            if (sec->IsRoundTurn())
                            {
                                glm::vec2 center = sec->RoundTurnCenter();
                                bool turnToLeft = 
                                    glm::dot(LeftNormal(sec->StraightDir(n)), glm::normalize(center - n->pos)) > 0;
                                
                                bool simpleTurn = (turnToLeft && leftSide) || (!turnToLeft && !leftSide);
                                if (simpleTurn)
                                {
                                    glm::vec2 hordeCenter = (sec->start->pos + sec->end->pos) * 0.5f;
                                    glm::vec2 dir = glm::normalize(center - hordeCenter);
                                    posReal = 
                                        X0Z(PointMapNormalizedToReal(m_heightmap, center + dir * distRoadToBasementNorm));
                                }
                                else
                                {
                                    glm::vec2 shift = n->pos - center;
                                    posReal = 
                                        X0Z(PointMapNormalizedToReal(m_heightmap, center + shift * 2.f));

                                    float angleRadians = SafeAcosf(glm::dot(glm::normalize(sec->startDir), glm::normalize(-sec->endDir)));
                                    angleRadians /= 2;
                                    float x = glm::length(sec->start->pos - sec->end->pos) * 0.5f;
                                    float l = x / sinf(angleRadians) / cosf(angleRadians);
                                    l *= (2 * preferRadiusNormalized + distRoadToBasementNorm) / (preferRadiusNormalized);
                                    glm::vec2 hordeCenter = (sec->start->pos + sec->end->pos) * 0.5f;
                                    glm::vec2 antiCenter = center + glm::normalize(hordeCenter - center) * l;

                                    protoBasement.points.push_back(XZ(posReal));
                                    posReal =
                                        X0Z(PointMapNormalizedToReal(m_heightmap, antiCenter));
                                }
                            }
                            else
                            {
                                RoadNode* prev = GetNode(path, i - 1);
                                RoadSection* secPrev = RoadSection::SectionByEdges(n, prev);

                                if (secPrev->IsRoundTurn())
                                {
                                    glm::vec2 center = secPrev->RoundTurnCenter();
                                    bool turnToLeft = 
                                        glm::dot(LeftNormal(secPrev->StraightDir(prev)), glm::normalize(center - prev->pos)) > 0;
                                    
                                    bool simpleTurn = (turnToLeft && leftSide) || (!turnToLeft && !leftSide);
                                    if (simpleTurn)
                                    {
                                        glm::vec2 hordeCenter = (secPrev->start->pos + secPrev->end->pos) * 0.5f;
                                        glm::vec2 dir = glm::normalize(center - hordeCenter);
                                        posReal = 
                                            X0Z(PointMapNormalizedToReal(m_heightmap, center + dir * distRoadToBasementNorm));
                                    }
                                    else
                                    {
                                        glm::vec2 shift = n->pos - center;
                                        posReal = 
                                            X0Z(PointMapNormalizedToReal(m_heightmap, center + shift * 2.f));
                                    }
                                }
                                else
                                {
                                    debug("SHOULD BE ROUND TURN!");
                                    throw std::exception{};
                                }
                            }
                        }
                        // debug("2_10");
                    }

                    protoBasement.points.push_back(XZ(posReal));
                }
                // debug("2_11");
                Polygon clearProtoBasement;
                for (int o = 0; o < protoBasement.points.size(); o++)
                {
                    if (Approx(protoBasement.points[o], protoBasement.GetPoint(o + 1), 0.001))
                        continue;
                    clearProtoBasement.points.push_back(protoBasement.points[o]);
                }

                // m_buildings.push_back(new Building(clearProtoBasement, 6));
                // debug("2_12");
                std::unique_ptr<Rect> inscribed = clearProtoBasement.Insribe();
                if (inscribed)
                {
                    bool success = true;
                    for (Basement blocked: innerBasements)
                    {
                        if (blocked.CheckRectIntersect(*inscribed) || blocked.IsPointInside(inscribed->center))
                        {
                            success = false;
                            break;
                        }
                    }
                    
                    if (success)
                    {
                        Basement basement = *inscribed;
                        basement.height = CalculateBasementHeight(basement);
                        innerBasements.emplace_back(basement);
                    }

                    // Polygon temp;
                    // temp.points = basement.GetPoints();
                    // m_buildings.push_back(new Building(temp, 10));
                }  
                // debug("2_13");
            }
        }
        reached.emplace(startSec);
    }
    
    // debug("4");
    std::function rand_gen
    {
        [] () 
        {
            // return a uniformly distributed random value
            return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
        }
    };
    std::function normalRandom
    {
        [&rand_gen] ()
        {
            // return a normally distributed random value
            double v1=rand_gen();
            double v2=rand_gen();
            return cosf(2*3.14*v2)*sqrt(-2.*log(v1));
        }
    };
    std::function rndNormalized
    {
        [&normalRandom] (double sigma, double Mi)
        {
            return normalRandom() * sigma + Mi;
        }
    };

    const int MAX_ATTEMPTS = BuildingsGenerator::OUT_MAX_ATTEMPTS; 
    const float MIN_AREA_REAL = BuildingsGenerator::OUT_MIN_AREA_REAL;
    const float MAX_AREA_REAL = BuildingsGenerator::OUT_MAX_AREA_REAL;
    const float MIN_SIDES_RATIO = BuildingsGenerator::OUT_MIN_SIDES_RATIO;
    const float MAX_SIDES_RATIO = BuildingsGenerator::OUT_MAX_SIDES_RATIO;
    const float NORMAL_DISP = BuildingsGenerator::OUT_NORMAL_DISP;
    const float OUT_BUILDINGS_EXTRA_SHIFT_REAL = BuildingsGenerator::OUT_OUT_BUILDINGS_EXTRA_SHIFT_REAL;


    int counter = 0;
    float extraShift = OUT_BUILDINGS_EXTRA_SHIFT_REAL / m_realSize[0];
    roadSections = m_roadScheme->graph.GetSections();
    std::vector<Basement> outBasements;
    while (counter < MAX_ATTEMPTS)
    {
        RoadSection* section = roadSections.at(rnd() * roadSections.size());
        std::unique_ptr<Rect> result(nullptr);

        if (section->IsStraight())
        {
            float centerPart = (float)rndNormalized(sqrt(NORMAL_DISP), 0);
            centerPart = centerPart * 0.5f + 0.5f;
            if (centerPart < 0 || centerPart > 1.f)
                continue;
            
            float sidesRatio = LERP(MIN_SIDES_RATIO, MAX_SIDES_RATIO, rnd());
            float area = LERP(MIN_AREA_REAL, MAX_AREA_REAL, rnd());
            float b = sqrtf((area / sidesRatio) - area);
            glm::vec2 size{area / b, b};
            if (size.x < size.y)
            {
                float tmp = size.x;
                size.x = size.y;
                size.y = tmp;
            }

            glm::vec2 center = LERP(section->start->pos, section->end->pos, centerPart);
            glm::vec2 normDir = LeftNormal(section->StraightDir());
            if (rnd() > 0.5)
                normDir *= -1.f;
            center += normDir * (extraShift + preferRadiusNormalized + distRoadToBasementNorm + size.y * 0.5f / m_realSize[0]);
            if (!m_heightmap->IsPointValid_Normalized(center))
            {
                counter++;
                continue;
            }

            center = PointMapNormalizedToReal(m_heightmap, center);

            glm::vec2 roadDir = section->StraightDir();
            if (glm::dot(roadDir, glm::vec2{0.f, 1.f}) < 0.f)
                roadDir *= -1.f;
            float angleRad = SafeAcosf(glm::dot(roadDir, glm::vec2{1.f, 0.f}));

            result = std::make_unique<Rect>(center, size, angleRad);

            bool bad = false;
            for (glm::vec2 _point: result->GetPoints())
            {
                if (!m_heightmap->IsPointValid_Normalized(PointRealToMapNormalized(m_heightmap, _point)))
                {
                    bad = true;
                    break;
                }
            }
            if (bad)
            {
                counter++;
                continue;
            }
        }
        else if (section->IsRoundTurn())
        {
            float centerPart = (float)rndNormalized(sqrt(NORMAL_DISP), 0);
            centerPart = centerPart * 0.5f + 0.5f;
            if (centerPart < 0 || centerPart > 1.f)
                continue;
            
            float sidesRatio = LERP(MIN_SIDES_RATIO, MAX_SIDES_RATIO, rnd());
            float area = LERP(MIN_AREA_REAL, MAX_AREA_REAL, rnd());
            float b = sqrtf((area / sidesRatio) - area);
            glm::vec2 size{area / b, b};
            if (size.x < size.y)
            {
                float tmp = size.x;
                size.x = size.y;
                size.y = tmp;
            }

            float l;
            {
                float angleRadians = SafeAcosf(glm::dot(glm::normalize(section->startDir), glm::normalize(-section->endDir)));
                angleRadians /= 2;
                float x = glm::length(section->start->pos - section->end->pos) * 0.5f;
                l = x / sinf(angleRadians) / cosf(angleRadians);
                l *= (2 * preferRadiusNormalized + distRoadToBasementNorm) / (preferRadiusNormalized);
            }
            glm::vec2 turnCenter = section->RoundTurnCenter();
            glm::vec2 cetnerDir = glm::normalize(LERP(
                glm::normalize(section->start->pos - turnCenter), 
                glm::normalize(section->end->pos - turnCenter),
                centerPart));
            
            glm::vec2 center = turnCenter + (extraShift + l + size.y * 0.5f / m_realSize[0]) * cetnerDir;
            if (!m_heightmap->IsPointValid_Normalized(center))
            {
                counter++;
                continue;
            }
            center = PointMapNormalizedToReal(m_heightmap, center);

            float angleRad = SafeAcosf(glm::dot(cetnerDir, glm::vec2{0.f, 1.f}));
            if (glm::dot(cetnerDir, glm::vec2{1.f, 0.f}) > 0.f)
                angleRad *= -1.f;

            result = std::make_unique<Rect>(center, size, angleRad);

            bool bad = false;
            for (glm::vec2 _point: result->GetPoints())
            {
                if (!m_heightmap->IsPointValid_Normalized(PointRealToMapNormalized(m_heightmap, _point)))
                {
                    bad = true;
                    break;
                }
            }
            if (bad)
            {
                counter++;
                continue;
            }
        }

        if (result)
        {
            bool success = true;

            for (int i = 0; i < 5; i++)
            {
                glm::vec2 p = (i == 4) ? result->center : result->GetPoint(i);
                p = PointRealToMapNormalized(m_heightmap, p);
                PopulationPoint pp = m_populationmap->Get(m_populationmap->PointNormalizedToInt(p));
                if (pp.population != Population::POPULATED)
                {
                    success = false;
                    break;
                }
            }
            

            if (success)
            {
                for (Polygon blocked: roadBlockBasements)
                {
                    if (blocked.CheckRectIntersect(*result) || blocked.IsPointInside(result->center))
                    {
                        success = false;
                        break;
                    }
                }
            }

            if (success)
            {
                for (Basement blocked: innerBasements)
                {
                    if (blocked.CheckRectIntersect(*result) || blocked.IsPointInside(result->center))
                    {
                        success = false;
                        break;
                    }
                }
            }

            if (success)
            {
            //   debug("3");      
            for (Basement blocked: outBasements)
                {
                    if (blocked.CheckRectIntersect(*result) || blocked.IsPointInside(result->center))
                    {
                        success = false;
                        break;
                    }
                }
            }
            
            if (success)
            {
                Basement b = *result;
                b.height = CalculateBasementHeight(b);
                outBasements.push_back(b);
            }
            
        }
        counter++;
    }

    // debug("5");
    for (Basement b: outBasements)
    {
        basements.push_back(b);
    }

    for (Basement b: innerBasements)
    {
        basements.push_back(b);
    }

    debug("Done!");
    
    return basements;
}

void Landscape::GenerateBuildings()
{
    std::vector<Basement> basements = GenerateBasements();
    
    debug("Generating buildings!");

    std::sort(
        basements.begin(), 
        basements.end(),
        [this] (const Basement& a, const Basement& b) 
        { 
            vec2Int aPos = m_populationmap->PointNormalizedToInt(PointRealToMapNormalized(m_biomemap, a.center));
            vec2Int bPos = m_populationmap->PointNormalizedToInt(PointRealToMapNormalized(m_biomemap, b.center));
            return m_populationmap->Get(aPos).intensity > m_populationmap->Get(bPos).intensity;
        }
    );

    std::vector<int> corruptedInds;
    for (size_t i = 0; i < basements.size(); i++)
    {
        int intersectedCount = 0;
        for (int j = 0; j < basements.size(); j++)
        {
            if (i == j)
                continue;

            if (basements[i].CheckRectIntersect(basements[j]))
            {
                intersectedCount++;
            }
        }  
        if (intersectedCount > 1)
        {
            corruptedInds.push_back(i);
            break;
        }
    }

    for (int i = corruptedInds.size() - 1; i >= 0; i--)
    {
        basements.erase(basements.begin() + corruptedInds[i]);
    }
    
    std::vector<Basement*> generatedBasemetns;
    Building::InitBuildingTypes();
    int scyscraperCount = int(BuildingsGenerator::SCYSCRAPERS_AMOUNT_PART * basements.size());
    scyscraperCount = std::min(scyscraperCount, BuildingsGenerator::SCYSCRAPERS_MAX_AMOUNT);
    for (Basement b: basements)
    {
        bool good = true;
        for (int i = 0; i < 4; i++)
        {
            glm::vec2 p = b.GetPoint(i);
            p = PointRealToMapNormalized(m_heightmap, p);
            if (!m_heightmap->IsPointValid_Normalized(p))
            {
                good = false;
            }
        }
        if (!good)
            continue;
        
        if (scyscraperCount > 0)
        {
            BUILDING_TYPE curType = BUILDING_TYPE::SCYSCRAPER;
            float minSize = Building::buildingTypes.at(curType).GetRealSize().x;
            if (b.size.x < minSize || b.size.y < minSize)
                continue;

            scyscraperCount--;
            float height = LERP(
                    BuildingsGenerator::SCYSCRAPER_MIN_HEIGHT_REAL, BuildingsGenerator::SCYSCRAPER_MAX_HEIGHT_REAL, rnd()
                );
            m_buildings.push_back(new Building(b, curType, height));
            generatedBasemetns.push_back(&b);
        }
        else
        {
            BUILDING_TYPE curType = BUILDING_TYPE::HIGH_RISE;
            float minSize = Building::buildingTypes.at(curType).GetRealSize().x;
            if (b.size.x < minSize || b.size.y < minSize)
                continue;

            float height = LERP(
                    BuildingsGenerator::PANEL_MIN_HEIGHT_REAL, BuildingsGenerator::PANEL_MAX_HEIGHT_REAL, rnd()
                );
            m_buildings.push_back(new Building(b, curType, height));
            generatedBasemetns.push_back(&b);
        }
    }

    AddBuildingsToTexture();

    debug("Done!");
}

float Landscape::CalculateBasementHeight(Basement& b)
{
    const bool DEBUG_BASEMENT = false;
    float basementHeight = (DEBUG_BASEMENT) ? -10000 : 10000;
    for (glm::vec2 point: b.GetPoints())
    {
        glm::vec2 p = PointRealToMapNormalized(m_heightmap, point);
        if (!m_heightmap->IsPointValid_Normalized(point))
            return 0.f;
        float h = m_heightmap->GetInterpolated_IntContinuous(
            m_heightmap->PointNormalizedToIntContinuous(PointRealToMapNormalized(m_heightmap, point))
        );
        basementHeight = (DEBUG_BASEMENT) 
            ? std::max(basementHeight, h)
            : std::min(basementHeight, h);
    }

    return basementHeight;
}

void Landscape::AddBuildingsToTexture()
{
    float widthToHeightRatio = m_roadScheme->widthToHeightRatio;
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
            SetColor(vec2Int{i, j}, color);
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
    for (RoadSection* section : m_roadScheme->graph.GetSections())
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

    for (unsigned j = 0; j < TEX_HEIGHT; j++) 
    {
        for (unsigned i = 0; i < TEX_WIDTH; i++) 
        {
            glm::vec2 coords = glm::vec2{(float)i / TEX_HEIGHT, (float)j / TEX_WIDTH};
            coords = PointMapNormalizedToReal(m_heightmap, coords);
            for (Building* building: m_buildings)
            {
                if (building->basement.IsPointInside(coords))
                {
                    unsigned color = (((94u << 16) | (94u << 8) | (75)) << 8) | 255;
                    SetColor(vec2Int{i, j}, color);
                    break;
                }
            }
        }
    }

    //Saving
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &m_roadScheme->GetTextureRef());
    glBindTexture(GL_TEXTURE_2D, m_roadScheme->GetTexture());
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