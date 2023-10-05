#pragma once
#include "cities_generator/global.h"
#include "../ShaderDataHolders/RenderableObjectDataHolder.h"
#include <libnoise/noise.h>
#include "HeightMap.h"
#include "BiomeMap.h"
#include <memory>
#include "WaterPlate.h"
#include "SurfaceMap.h"
#include "RoadScheme.h"
#include "../ShaderDataHolders/Archetypes/BasicGLTFModel.h"
#include "RoadModel.h"
#include "PopulationMap.h"
#include "Building.h"

enum MAP
{
    NONE = 0,
    HEIGHT = 1,
    BIOME = 2,
    SURFCE = 3,
    ROAD = 4,
    POPULATION = 5,
    ENUM_SIZE = 6
};
MAP& operator++(MAP& map);


class Landscape : public RenderableObjectDataHolder
{

    public:
        Landscape();
        void Init();
        void Update();
        virtual void Render(RENDER_MODE mode);
        void RenderMap(MAP map);
        ~Landscape();
        void SwapHeightMaps(bool recreate = true);
        WaterPlate* GetWaterPlate(); 
        bool IsRenderingWater();
        float GetWaterLevel(); 
        float GetWaterLevelSafeDelta();
        void RenderRoad(RENDER_MODE mode); 
        void RenderBuildings(RENDER_MODE mode); 

    private:         
        void GenerateHeightmap_PERLIN();
        void GenerateHeightmap_DIAMONDSQUARE(bool = false);
        void ApplyErosionToHeightmap();
        void RecreateMesh();
        void ReadNormal(int x, int y, GLfloat output[]);
        glm::vec3 GetNormal(int x, int y);
        float GetHillsGrassRatio(glm::vec3 pos, glm::vec3 normal);
        float GetHillsSlopeRatio(glm::vec3 pos, glm::vec3 normal);

        struct HillsGenerator
        {
            HillsGenerator();
            float GetGrassRatioByCos(float flatAngleCos);

            private:
                float m_fromDegrees, m_toDegrees, m_minRockRatio, m_maxRockRatio;
                float m_bordersCos[2];
                static constexpr float FREQUENCY = 0.03;
                static constexpr int OCTAVE_COUNT = 7;
                static constexpr float PERSISTANCE = 0.5;
                static constexpr float LACUNARITY = 2;
                static constexpr float MAX_HEIGHT = 30;
                static constexpr float SEMI_MOUNTAIN_AREA_REAL = 10;

                static constexpr float MAX_MOUNTAIN_GRASS_TEXTURE_ALPHA = 0.4f;

            friend void Landscape::GenerateHeightmap_PERLIN();
            friend float Landscape::GetHillsGrassRatio(glm::vec3 pos, glm::vec3 normal);
        };

        class PlainGenerator
        {
            static constexpr float FREQUENCY = 0.03;
            static constexpr int OCTAVE_COUNT = 2;
            static constexpr float PERSISTANCE = 0.1;
            static constexpr float LACUNARITY = 5;
            static constexpr float MAX_HEIGHT = 5;
            static constexpr float MAX_DISTORTION = 7;
            static constexpr float CITY_HEIGHT_CHANGE_SUPRESSION_FACTOR = 0.5;
            static constexpr float SEMI_CITY_AREA_REAL = 15.f;

            static constexpr float DIAMOND_SQUARE_HELP_MAP_SIZE = 256;

            static constexpr float DIAMOND_SQUARE_SPREAD = 0.25; //0.15
            friend void Landscape::GenerateHeightmap_PERLIN();
            friend void Landscape::GenerateHeightmap_DIAMONDSQUARE(bool);
        };

        template<class T>
        glm::vec2 PointRealToMapNormalized(DefaultMap<T>* map, glm::vec2 coords, bool safe = false);
        template<class T>
        glm::vec2 PointMapNormalizedToReal(DefaultMap<T>* map, glm::vec2 coords);
        template<class T>
        glm::vec2 PointMapIntToReal(DefaultMap<T>* map, vec2Int coords);
        template<class T>
        float GetMapCellRealSize(DefaultMap<T>* map);

        float GRASS_SCALE = 0.8;  
        GLuint m_grassTexture;
        float SOIL_SCALE = 1.5;
        GLuint m_soilTexture;
        float SAND_SCALE = 1.5;
        GLuint m_sandTexture;

        HeightMap *m_heightmap, *m_savedHeightmap;
        CGenBiomeMap* m_biomemap;
        SurfaceMap* m_surfacemap;
        RoadScheme* m_roadScheme;
        PopulationMap* m_populationmap = nullptr;;
        std::vector<Building*> m_buildings;
        

        float m_realSize[2];
        float waterLevel;
        HillsGenerator hillsGenerator;  

        class Erosion
        {
            //MAP 400 x 400 and density 0.8 for heightmap
            static constexpr int MAX_STEPS = 30; //500
            //relative to heightmap cell size
            static constexpr float STEP_LENGTH = 0.5; //0.1
            static constexpr float DENSITY = 5; // PER 1x1 real size
            // static constexpr unsigned DROPS_COUNT = 150000;
            static constexpr float INERTIA = 0.2; // 1.0 * STEP_LENGTH;
            // static constexpr float MIN_SLOPE = 0.02 * STEP_LENGTH;
            static constexpr float GRAVITY = 3.0;
            static constexpr float CAPACITY = 1.5; //2
            // потестить второй хейтмап, исправить формулу инерции и исправить багу со скачками высоты
            static constexpr float BASE_EROSION = 0.3; 
            static constexpr float BASE_DEPOSITION = 0.05;
            static constexpr float EROSION_RADIUS = 3;
            // MAX ECityAffectOnHeightmapVAPORATION = 1.0f - float(pow(BASE_EVAPORATION, 1 / MAX_STEPS));
            static constexpr float BASE_EVAPORATION = 0.01;
            static constexpr float NON_EVAPORATE_ANGLE_DEGREES = 30.0;

            // //Map 70 x 70 density 2.8
            // static constexpr int MAX_STEPS = 100; //500
            // //relative to heightmap cell size
            // static constexpr float STEP_LENGTH = 0.5; //0.1
            // static constexpr float DENSITY = 34; // PER 1x1 real size
            // // static constexpr unsigned DROPS_COUNT = 150000;
            // static constexpr float INERTIA = 0.2; // 1.0 * STEP_LENGTH;
            // // static constexpr float MIN_SLOPE = 0.02 * STEP_LENGTH;
            // static constexpr float GRAVITY = 3.0;
            // static constexpr float CAPACITY = 1.5; //2
            // // потестить второй хейтмап, исправить формулу инерции и исправить багу со скачками высоты
            // static constexpr float BASE_EROSION = 0.3; 
            // static constexpr float BASE_DEPOSITION = 0.05;
            // static constexpr float EROSION_RADIUS = 6;
            // // MAX EVAPORATION = 1.0f - float(pow(BASE_EVAPORATION, 1 / MAX_STEPS));
            // static constexpr float BASE_EVAPORATION = 0.01;
            // static constexpr float NON_EVAPORATE_ANGLE_DEGREES = 30.0;

            struct ErosionDroplet
            {
                glm::vec2 positionReal;
                glm::vec2 direction;
                float velocity;
                float sediment;
                float water;
                int stepCounter;

                ErosionDroplet();
                void Reset(glm::vec2 position);
            };


            static void CreateDebugDrop(glm::vec3 pos);
            static void MoveDebugDropTo(glm::vec3 pos);
            static void DebugDropSetColor(glm::vec3 col);

            friend void Landscape::ApplyErosionToHeightmap();
            friend void Landscape::Update();
        };
        bool IsPointErodable(glm::vec2 heightmapTC);
        bool IsBiomeErodable(Biome b);


        void GenerateBiomemap();

        class WaterBiomeGenerator
        {
            static constexpr float COAST_STRAIGHT_SEA_SIZE = 0.5;
            static constexpr float COAST_STRAIGHT_MIDPOINT_MAX_ANGLE_DEGREES = 20;
            static constexpr int COAST_STRAIGHT_MIDPOINT_DIVISION_TIMES = 7;

            friend void Landscape::GenerateBiomemap();
        };

        struct MountainBiomeGenerator
        {
            static constexpr float CENTER_PLAIN_RADIUS = 0.32f;//0.35f; //0.5 is inscribed circle 
            static constexpr float MOUNTAIN_MIN_PROBABILITY = 0.3; //at center radius
            static constexpr float MOUNTAIN_MAX_PROBABILITY = 0.8; //at radius sqrt(2) / 2
            // static constexpr float POINTS_DENSITY = 0.003; //per 1x1 real size square
            static constexpr int PERLIN_OCTAVE_COUNT = 2;
            static constexpr float PERLIN_FREQUENCY_REAL = 0.04;
            static constexpr float PERLIN_PERSISTANCE = 0.5;
            static constexpr float PERLIN_LACUNARITY = 2;
            static constexpr float PERLIN_AFFECT_RATIO = 0.17;
            static constexpr float POINTS_AMOUNT = 23;
        };


        float CalculateWaterLevel();
        void ApplyShoreImmerse();
        void CreateWaterPlate();
        float GetWaterTextureRatio(glm::vec3 pos);
        bool IsPointSedimentDropable(glm::vec2 heightmapTC);
        struct WaterGenerator
        {
            static constexpr float FREQUENCY = 0.05;
            static constexpr int OCTAVE_COUNT = 2;
            static constexpr float PERSISTANCE = 0.05;
            static constexpr float LACUNARITY = 7;
            static constexpr float HEIGHT_DELTA = 9;
            static constexpr float SEA_DEPTH = 5;

            static constexpr float SHORE_REAL_SIZE = 10.0f;
            static constexpr float GROUND_HEIGHT_ABOVE_WATER = 1.0f;
            static constexpr float SHALLOW_WATER_REAL_SIZE = 25.0f;

            static constexpr float SAND_TEXTURE_MIN_TO_SHORE = 0.4f;
            static constexpr float SAND_TEXTURE_MAX_TO_SHORE = 0.7f;
        };       
        std::unique_ptr<WaterPlate> waterPlate; 

        void GenerateRoad();
        RoadSlice CreateSliceFromSchemePoint (glm::vec2, glm::vec2, float);
        std::array<glm::vec2, 2> GetPlanarPointsForCrossroadSection(RoadNode *, RoadSection &);
        RoadSlice CreateSliceFromCrossroad_Uncorrected(RoadSection &, bool);
        std::array<std::vector<RoadSection*>, 2> GetLeftRightTurns(RoadNode* crossroad, RoadSection& section);
        
        RoadModel* m_road;
        public:
            struct RoadModelGenerator
            {
                static constexpr float ROAD_FULL_REAL_WIDTH = 6.0f;
                static constexpr float SEGMENT_STRAIGHT_REAL_LENGTH = 5.0f;
                static constexpr float SEGMENT_ROUND_MAX_ANGLE_DEGREES = 10.0f;
                static constexpr float ROAD_HEIGHT_ABOVE_HEIGHTMAP = 0.4f;
                static constexpr float CROSSROAD_CORRECTION_MIN_ANGLE_DEGREES = 150.f;
                static constexpr float TURN_PREFER_RADIUS_RATIO = 1.5f;
                static float GET_TURN_MIN_RADIUS_RATIO();
                static float GET_MIN_ANGLE_BETWEEN_ROADS_DEGREES();
            };
        private:

        void GeneratePopulationMap();
        float GetPopulationRatioByRealPos(glm::vec2);
        public:
            struct PopulationGenerator
            {
                static bool IsBiomePopulatable(Biome b);
                static float GetPointPopulationByLeftPart(float populationLeft);
                static float GetPointCenterDistancePartByLeftPart(float populationLeft);
                
                static constexpr float POPULATION_VOLUME_REAL = 45000; //FINAL 45000
                
                static constexpr int POPULATION_PERLIN_OCTAVE_COUNT = 2;
                static constexpr float POPULATION_PERLIN_FREQUENCY_REAL = 0.07;
                static constexpr float POPULATION_PERLIN_PERSISTANCE = 0.5;
                static constexpr float POPULATION_PERLIN_LACUNARITY = 2;

                static constexpr float POPULATION_RATIO_DISTANCE_TO_CENTER = 10.f; 
                static constexpr float POPULATION_RATIO_LOWLAND = 40.f;
                static constexpr float POPULATION_RATIO_DISTANCE_TO_SHORE = 330.f;
                static constexpr float POPULATION_RATIO_PERLIN = 150.f;
                static constexpr float POPULATION_RATIO_BIOME_POSSIBILITY = -1000000.f;
                static constexpr float POPULATION_RATIO_POSITIVE_SUM = 
                    POPULATION_RATIO_DISTANCE_TO_CENTER +
                    POPULATION_RATIO_LOWLAND +
                    POPULATION_RATIO_DISTANCE_TO_SHORE +
                    POPULATION_RATIO_PERLIN;
            };

            struct RoadSchemeGenerator
            {
                static constexpr int HIGHWAY_RAYS_PER_SEGMENT = 14;
                static constexpr float MIN_DISTANCE_BETWEEN_CROSSROADS = 25.0f;
                static constexpr float MAX_GENERATION_SEGMENT_LENGTH = 100.0f;
                static constexpr float MIN_GENERATION_SEGMENT_LENGTH = 45.0f;
                static constexpr float RAY_GENERATION_STEP_REAL = 1.5f;
                static constexpr float HIGHWAY_MIN_SEGMENT_SUBDIVISION_PART = 0.3f;
                static constexpr int MAX_CROSSROAD_ROADS = 4;
                static constexpr float MERGE_TO_CROSSROAD_DISTANCE_GROWTH = 15.0f;
                static constexpr float MERGE_TO_CROSSROAD_DISTANCE_LINKING = 20.0f;
                static constexpr float MERGE_TO_MIDROAD_DISTANCE_LINKING = 20.0f;
                static constexpr float INTERSECTING_ROADS_DISTANCE_REAL = std::min(
                        MERGE_TO_CROSSROAD_DISTANCE_GROWTH, 
                        MERGE_TO_MIDROAD_DISTANCE_LINKING
                    );
                static constexpr float RAY_GENERATION_LONG_PROMOTION_RATIO = 0.001f; //0 - no promotion, inf - max promotion
                static constexpr float SECTION_FINAL_MAX_POPULATION_RATIO = 0.015;
                static constexpr float DEADEND_SECTION_BIOME_RADIUS_REAL = 6.f;
                static constexpr float DEADEND_SECTION_BLOCKED_RAYS_MIN_PART = 0.25f;
                static constexpr float EXTEND_DEADEND_ANGLE_DISPERSION_DEGREES = 180; //from 0 to 360

                //maximum ratio of lengths: (extension) / (existing)
                static constexpr float EXTEND_DEADEND_MAX_LENGTH_RATIO = 1.5; 
            };
            bool IsCrossroadJoinable(RoadNode* cross);
            bool IsCrossroadJoinable(RoadNode* cross, glm::vec2 posNew);
            bool IsSectionJoinable(RoadSection* sec, float ratio, glm::vec2 posNew);
            void RoundOffTurns();

            struct BuildingsGenerator
            {
                static constexpr float DIST_BETWEEN_ROAD_AND_BASEMENT_REAL = 1.f;

                static constexpr int OUT_MAX_ATTEMPTS = 200;
                static constexpr float OUT_MIN_AREA_REAL = 200;
                static constexpr float OUT_MAX_AREA_REAL = 500;
                static constexpr float OUT_MIN_SIDES_RATIO = 0.25; // a and (1 - a)
                static constexpr float OUT_MAX_SIDES_RATIO = 0.5;
                static constexpr float OUT_NORMAL_DISP = 0.2;
                static constexpr float OUT_OUT_BUILDINGS_EXTRA_SHIFT_REAL = 0.5;

                static constexpr int INSCRIBE_FIND_SQUARE_POINTS = 1000;
                static constexpr int INSCRIBE_SIZE_ATTEMPTS = 15;
                static constexpr int INSCRIBE_OVERAL_ATTEMPTS_AMOUNT = 10000;
                static constexpr float INSCRIBE_MAX_SIDE_RATIO = 0.3;
                static constexpr float INSCRIBE_MIN_AREA_REAL = 50;
                static constexpr float INSCRIBE_MAX_AREA_START_PART = 0.9;

                static constexpr float SCYSCRAPERS_AMOUNT_PART = 0.2;
                static constexpr int SCYSCRAPERS_MAX_AMOUNT = 4;
                static constexpr float SCYSCRAPER_MIN_HEIGHT_REAL = 100;
                static constexpr float SCYSCRAPER_MAX_HEIGHT_REAL = 150;
                static constexpr float PANEL_MIN_HEIGHT_REAL = 20;
                static constexpr float PANEL_MAX_HEIGHT_REAL = 40;
            };
            void GenerateBuildings();
            std::vector<Basement> GenerateBasements();
            float CalculateBasementHeight(Basement& b);
            void AddBuildingsToTexture();
};