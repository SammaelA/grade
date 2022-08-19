#pragma once
#include <map>
#include <string>
#include "tinyEngine/resources.h"

class Timestamp
{
public:
    Timestamp(bool need_debug = true);
    void start(std::string name);
    void end(std::string name);
    void resolve();
private:
    struct Stamp
    {
        GLuint   startQuery = 0, endQuery = 0;
        GLuint64 startTime = 0, endTime = 0, avTime = 0;  
        bool started = false, ended = false, first = true;
    };

    static constexpr float spike_mul = 3.0;
    static constexpr float av_mix = 0.975;
    static constexpr int debug_time = 100;
    bool need_debug = false;
    uint64_t ticks = 0;
    std::map<std::string, Stamp> TS = {};
};