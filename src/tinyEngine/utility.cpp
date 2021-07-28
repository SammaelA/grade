#include "utility.h"
#include <map>
#include <new>
#include <vector>
#include "../generated_tree.h"

#include "malloc.h"

#define DEBUG_LEVEL 100
void debug(const char *__restrict __fmt, ...)
{
    va_list args;
    va_start(args, __fmt);
    vfprintf(stdout, __fmt, args);
    va_end(args);
}
void debugnl()
{
    debug("\n");
}
void debug(const char *__restrict __fmt, glm::vec2 vec)
{
    debug("[%f %f]",vec.x,vec.y);
}
void debugl(uint level, const char *__restrict __fmt, glm::vec2 vec)
{
    if (level && DEBUG_LEVEL != level)
        return;
    debug(__fmt, vec);
}
void debug(const char *__restrict __fmt, glm::vec3 vec)
{
    debug("[%f %f %f]",vec.x,vec.y,vec.z);
}
void debugl(uint level, const char *__restrict __fmt, glm::vec3 vec)
{
    if (level && DEBUG_LEVEL != level)
        return;
    debug(__fmt, vec);
}
void debug(const char *__restrict __fmt, glm::vec4 vec)
{
    debug("[%f %f %f %f]",vec.x,vec.y,vec.z,vec.w);
}
void debugl(uint level, const char *__restrict __fmt, glm::vec4 vec)
{
    if (level && DEBUG_LEVEL != level)
        return;
    debug(__fmt, vec);
}
void debugl(uint level, const char *__restrict __fmt, ...)
{
    if (DEBUG_LEVEL > 0 && level && DEBUG_LEVEL != level)
        return;
    va_list args;
    va_start(args, __fmt);
    vfprintf(stdout,__fmt, args);
    va_end(args);
}
void logerr(const char *__restrict __fmt, ...)
{
    va_list args;
    va_start(args, __fmt);
    vfprintf(stderr, __fmt, args);
    va_end(args);
    fprintf(stderr,"\n");
}

std::map<std::string, AllocData> alloc_info;

void print_alloc_info()
{
    for (auto &p : alloc_info)
    {
        logerr("%s : %ld %ld %ld", p.first.c_str(), p.second.active_allocs, p.second.all_allocs, p.second.alloc_size);
    }
    Countable::get_cur_count();
}
int full_count = 0;
Countable::Countable() 
{
    count++;
    full_count++;
};
Countable::~Countable() {count--;};
void Countable::get_cur_count()  {logerr("count = %d %d",count, full_count);}
long Countable::count = 0;