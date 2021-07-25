#pragma once
#include <cstdio>
#include <glm/glm.hpp>
#include <stdarg.h>
#include <string>
#include <map>

#define PI 3.14159265f  
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define CLAMP(a,b,c) (MIN(MAX((a),(b)),(c)))
#define SQR(a) ((a)*(a))


void debugnl();
void debug(const char *__restrict __fmt, glm::vec2 vec);
void debugl(uint level, const char *__restrict __fmt, glm::vec2 vec);
void debug(const char *__restrict __fmt, glm::vec3 vec);
void debugl(uint level, const char *__restrict __fmt, glm::vec3 vec);
void debug(const char *__restrict __fmt, glm::vec4 vec);
void debugl(uint level, const char *__restrict __fmt, glm::vec4 vec);
void debugl(uint level, const char *__restrict __fmt, ...);
void debug(const char *__restrict __fmt, ...);
void logerr(const char *__restrict __fmt, ...);

struct AllocData
{
    long active_allocs = 0;
    long all_allocs = 0;
    long alloc_size = 0;
};
extern std::map<std::string, AllocData> alloc_info;

template <typename T>
T *safe_new(int count, std::string tag)
{   
    auto pos = alloc_info.find(tag);
    if (pos == alloc_info.end())
    {
        AllocData ad;
        ad.active_allocs = 1;
        ad.all_allocs = 1;
        ad.alloc_size = sizeof(T)*count;

        alloc_info.emplace(tag,ad);
    }
    else 
    {
        pos->second.active_allocs++;
        pos->second.all_allocs++;
        pos->second.alloc_size += sizeof(T)*count;
    }
    return new T[count];
}

template <typename T>
T *safe_new(std::string tag)
{
    return safe_new<T>(1, tag);
}

template <typename T>
void safe_delete(T *ptr, std::string tag)
{
    if (!ptr)
        return;
    
    auto pos = alloc_info.find(tag);
    if (pos == alloc_info.end())
    {
        AllocData ad;
        ad.active_allocs = -1;
        ad.all_allocs = 0;

        alloc_info.emplace(tag,ad);
    }
    else 
    {
        pos->second.active_allocs--;
    }
    delete[] ptr;
}

void print_alloc_info();

class Countable
{
public:
    static long count;
    Countable();
    ~Countable();
    static void get_cur_count();
};