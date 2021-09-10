#include "utility.h"
#include <map>
#include <new>
#include <vector>
#include "../generated_tree.h"

#include "malloc.h"

int debug_level = 11;
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
    debug("[%f %f] %s",vec.x,vec.y,__fmt);
}
void debugl(uint level, const char *__restrict __fmt, glm::vec2 vec)
{
    if (level && debug_level != level)
        return;
    debug(__fmt, vec);
}
void debug(const char *__restrict __fmt, glm::vec3 vec)
{
    debug("[%f %f %f] %s",vec.x,vec.y,vec.z, __fmt);
}
void debugl(uint level, const char *__restrict __fmt, glm::vec3 vec)
{
    if (level && debug_level != level)
        return;
    debug(__fmt, vec);
}
void debug(const char *__restrict __fmt, glm::vec4 vec)
{
    debug("[%f %f %f %f] %s",vec.x,vec.y,vec.z,vec.w, __fmt);
}
void debugl(uint level, const char *__restrict __fmt, glm::vec4 vec)
{
    if (level && debug_level != level)
        return;
    debug(__fmt, vec);
}
void debugl(uint level, const char *__restrict __fmt, ...)
{
    if (debug_level > 0 && level && debug_level != level)
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
void print_FB_status(GLuint status)
{
    if (status == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
        debugl(9,"GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT\n");
    else if (status == GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS)
        debugl(9,"GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS\n");
    else if (status == GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT)
        debugl(9,"GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT\n");
    else if (status == GL_FRAMEBUFFER_UNSUPPORTED)
        debugl(9,"GL_FRAMEBUFFER_UNSUPPORTED\n");
    else  debugl(9,"GL_FRAMEBUFFER_INCOMPLETE %#010x\n",status);
}
std::map<std::string, AllocData> alloc_info;
long counts[256];
long full_counts[256];
void print_alloc_info()
{
    for (auto &p : alloc_info)
    {
        logerr("%s : %ld %ld %ld", p.first.c_str(), p.second.active_allocs, p.second.all_allocs, p.second.alloc_size);
    }
    for (int i=0;i<256;i++)
    {
        if (counts[i] || full_counts[i])
            logerr("[%d]=%d %d",i,counts[i], full_counts[i]);
    }
    Countable::get_cur_count();
}
int full_count = 0;
Countable::Countable(int num) 
{
    countable_type_num = num;
    count++;
    full_count++;
    counts[countable_type_num]++;
    full_counts[countable_type_num]++;
};
Countable::Countable(const Countable &c)
{
    countable_type_num = c.countable_type_num;
    count++;
    full_count++;
    counts[countable_type_num]++;
    full_counts[countable_type_num]++;
}
Countable::Countable(Countable &&c)
{
    countable_type_num = c.countable_type_num;
    count++;
    full_count++;
    counts[countable_type_num]++;
    full_counts[countable_type_num]++;
}
Countable &Countable::operator=(Countable &&c)
{
    countable_type_num = c.countable_type_num;
}
Countable &Countable::operator=(Countable &c)
{
    countable_type_num = c.countable_type_num;
}
Countable::~Countable() {count--;counts[countable_type_num]--;};
void Countable::get_cur_count()  {logerr("count = %d %d",count, full_count);}
long Countable::count = 0;