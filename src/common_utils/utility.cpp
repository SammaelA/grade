#include "utility.h"
#include <map>
#include <new>
#include <vector>
#include <boost/filesystem.hpp>
#include <iostream>

#include "malloc.h"

int debug_level = 1000;
int show_progress = 0;

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
typedef std::chrono::steady_clock::time_point TP;
void print_time_interval(double ms)
{
    if (ms < 1)
    {
        debug("%f ms",(float)ms);
        return;
    }
    if (ms < 10)
    {
        debug("%.2f ms",(float)ms);
        return;
    }
    if (ms < 100)
    {
        debug("%.1f ms",(float)ms);
        return;
    }
    float s = ms/1000;
    if (s < 1)
    {
        debug("%.3f seconds",(float)s);
        return;
    }
    if (s < 60)
    {
        debug("%.2f seconds",(float)s);
        return;
    }
    float m = s/60;
    if (m < 60)
    {
        debug("%.1f minutes",(float)m);
        return;
    }
    float h = m/60;
    debug("%.1f hours",(float)h);
    return;
}
ProgressBar::ProgressBar(std::string _process, long _count, std::string _type, bool _estimate_time)
{
    process = _process;
    count = _count;
    type = _type;
    estimate_time = _estimate_time;

    if (show_progress)
    {
        debug("[P] Starting %s: %d %s\n", process.c_str(), count, type.c_str());
        if (estimate_time)
        {
            t_start = std::chrono::steady_clock::now();
            t_prev = t_start;
        }
    //float frame_time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t_prev).count();
    }
}
ProgressBar::~ProgressBar()
{
    finish();
}
void ProgressBar::finish()
{
    if (!finished)
    {
        if (show_progress)
        {
            debug("[P] Finished %s ", process.c_str());
            if (estimate_time)
            {
                auto t_finish = std::chrono::steady_clock::now();
                double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_start).count();
                debug("took ");
                print_time_interval(ms);
            }
            debugnl();
        }
        finished = true;
    }
}
void ProgressBar::iter(long n)
{
    if (show_progress)
    {
        debug("[P] %s %d/%d ", process.c_str(), n, count);
        if (estimate_time && n > prev_iter)
        {
            auto t_now = std::chrono::steady_clock::now();
            double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_now - t_prev).count();
            ms = ms / (n - prev_iter);
            float q = MAX(t_estimate_q, 1 - 10.0/count);
            at_per_iter = at_per_iter > 0 ? (t_estimate_q*at_per_iter + (1 - t_estimate_q)*ms) : ms;
            debug("time left ", at_per_iter);
            print_time_interval(at_per_iter*(count - n));
            prev_iter = n;
            t_prev = t_now;
        }
        else if (estimate_time)
        {
            t_prev = std::chrono::steady_clock::now();
        }
        debugnl();
    }
}

bool prepare_directory(std::string &save_path)
{
  bool status = true;
  try
  {
    if (boost::filesystem::exists(save_path))
    {
      if (boost::filesystem::is_directory(save_path))
      {
        printf("replacing previous save\n");
        boost::filesystem::remove_all(save_path);
      }
      else
      {
        logerr("path %s represents existing file. Can not save here", save_path.c_str());
        status = false;
      }
    }
    if (status)
    {
      boost::filesystem::create_directory(save_path);
      boost::filesystem::permissions(save_path, boost::filesystem::perms::all_all);
    }
  }
  catch (const std::exception &e)
  {
    status = false;
    std::cerr << e.what() << '\n';
  }

  return status;
}