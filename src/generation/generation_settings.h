#pragma once
#include "core/tree.h"
#include "core/body.h"
#include "graphics_utils/impostor.h"

enum GenerationTask
{
    GENERATE = 1,
    CLUSTERIZE = 2,
    BILLBOARDS = 4,
    IMPOSTORS = 8,
    IMPOSTOR_FULL_GROVE = 16,
    SYNTS = 32,
    MODELS = 64
};
#define ALL_GENERATION_TASKS (2*MODELS - 1)
#define MINIMUM_FOR_RENDER (GENERATE|CLUSTERIZE|MODELS)
struct GroveGenerationData
{
    std::vector<TreeTypeData> types;
    int trees_count;
    int synts_count;
    int synts_precision;
    float clustering_max_individual_distance = 0.7;
    ImpostorBaker::ImpostorGenerationParams impostor_generation_params;
    Quality bill_1_quality = Quality::MEDIUM;
    Quality bill_2_quality = Quality::LOW;
    unsigned task = GENERATE | CLUSTERIZE | BILLBOARDS | IMPOSTORS | MODELS;
    glm::vec3 pos;
    glm::vec3 size;
    std::string name;
    std::vector<Body *> obstacles;
};