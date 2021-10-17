#pragma once
#include <cstdio>
#include <glm/glm.hpp>
#include <stdarg.h>
#include <string>
#include <map>
#include <GL/glew.h>
#include <boost/preprocessor.hpp>
#include <chrono>

#define X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE(r, data, elem)    \
    case elem : return BOOST_PP_STRINGIZE(elem);

#define DEFINE_ENUM_WITH_STRING_CONVERSIONS(name, enumerators)                \
    enum name {                                                               \
        BOOST_PP_SEQ_ENUM(enumerators)                                        \
    };                                                                        \
                                                                              \
    inline const char* ToString(name v)                                       \
    {                                                                         \
        switch (v)                                                            \
        {                                                                     \
            BOOST_PP_SEQ_FOR_EACH(                                            \
                X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE,          \
                name,                                                         \
                enumerators                                                   \
            )                                                                 \
            default: return "[Unknown " BOOST_PP_STRINGIZE(name) "]";         \
        }                                                                     \
    }
    
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
void print_FB_status(GLuint status);
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
    unsigned char countable_type_num = 0;
    Countable(int num = 0);
    Countable(const Countable &c);
    Countable(Countable &&c);
    ~Countable();
    Countable &operator=(Countable &c);
    Countable &operator=(Countable &&c);
    static void get_cur_count();
};

extern int debug_level;
extern int show_progress;
class ProgressBar
{
public:
    ProgressBar(std::string _process, long _count, std::string _type, bool _estimate_time = true);
    ~ProgressBar();
    void iter(long n);
    void finish();
    static void single_message(std::string message)
    {
        if (show_progress)
            debug("[P] %s\n", message.c_str());
    }
private:
    std::string process;
    long count;
    std::string type;
    bool estimate_time;
    std::chrono::steady_clock::time_point t_start;
    std::chrono::steady_clock::time_point t_prev;
    double at_per_iter = 0;
    long prev_iter = 0;
    float t_estimate_q = 0.9;
    bool finished = false;
};
