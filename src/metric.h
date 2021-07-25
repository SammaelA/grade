#pragma once
#include <functional>
#include "parameter.h"
#include "grove.h"
#include <string>
#include "texture_manager.h"

class Metric
{
public:
    virtual double get(GrovePacked &g) = 0;
};

class CompressionMetric : public Metric
{
public:
    virtual double get(GrovePacked &g) override;
};

class ImpostorMetric : public Metric
{
public:
    ImpostorMetric(Texture &reference_image);
    ~ImpostorMetric();
    virtual double get(GrovePacked &g) override;
private:
    Texture reference;
    unsigned char *reference_raw = nullptr;
    int w, h;
    glm::vec3 leaves_color = glm::vec3(0,1,0);
    glm::vec3 empty_color = glm::vec3(0,0,0);
    glm::vec3 wood_color = glm::vec3(75.0/256, 37.5/256, 0);
    int get_type(unsigned char r, unsigned char g);
    float diff(unsigned char *imp, unsigned char *reference, int imp_tex_w, int imp_tex_h,
               int imp_offset_w = 0, int imp_offset_h = 0);
};

class DummyMetric : public Metric
{
public:
    virtual double get(GrovePacked &g) override;
};