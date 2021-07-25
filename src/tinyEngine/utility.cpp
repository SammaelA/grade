#include "utility.h"
#include <map>
#include <new>
#include <vector>

#include "malloc.h"

#define DEBUG_LEVEL 10
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


#define Mb 1024*1024

/*struct Buffer
{
    std::vector<unsigned char *> chunks;
    int pos = 0;
    int chunk = 0;
};*/
long new_count = 0;

std::map<uint64_t, std::size_t> large_allocs;

/*void * operator new(std::size_t n) //throw(std::bad_alloc)
{
  new_count++;
  //if (new_count % 10000 == 0)
  //  logerr("new count %d",new_count);
  void *ptr = malloc(n);
  if (n > 100000 && large_allocs.size() < 10000)
  {
      large_allocs.emplace((uint64_t)ptr,n);
  }
  return ptr;
  
}
void operator delete(void * p) //throw()
{
  auto iter = large_allocs.find((uint64_t)p);
  if (iter != large_allocs.end())
    large_allocs.erase(iter);
  if (p)
    free(p);
}
void operator delete(void * p, std::size_t n) //throw()
{
  auto iter = large_allocs.find((uint64_t)p);
  if (iter != large_allocs.end())
    large_allocs.erase(iter);
  if (p)
    free(p);
}

class MemoryManager
{
	struct AllocInfo
	{
		const char* source;
		int         line;
		int         size;
		void*       data;
	};

	static MemoryManager   instance;
	std::vector<AllocInfo> allocations;

public:
	static void RegAllocation(void* ptr, int size, const char* source, int line)
	{
		AllocInfo info;
		info.data = ptr;
		info.size = size;
		info.source = source;
		info.line = line;

		instance.allocations.push_back(info);
	}

	static void UnregAllocation(void* ptr)
	{
		for (std::vector<AllocInfo>::iterator it = instance.allocations.begin(); it != instance.allocations.end(); ++it)
		{
			if (it->data == ptr)
			{
				instance.allocations.erase(it);
				return;
			}
		}
	}

	static void PrintInfo()
	{
		for (std::vector<AllocInfo>::iterator it = instance.allocations.begin(); it != instance.allocations.end(); ++it)
			printf("0x%x: %i bytes at %s: %i", it->data, it->size, it->source, it->line);
	}
};

MemoryManager MemoryManager::instance;

void* operator new(size_t size, const char* location, int line)
{
	void* res = ::operator new(size);
	MemoryManager::RegAllocation(res, size, location, line);
	return res;
}

void  operator delete(void* allocMemory, const char* location, int line)
{
	MemoryManager::UnregAllocation(allocMemory);
}

void  operator delete(void* allocMemory)
{
	MemoryManager::UnregAllocation(allocMemory);
}

*/
void print_alloc_info()
{
    for (auto &p : alloc_info)
    {
        logerr("%s : %d %d %d", p.first.c_str(), p.second.active_allocs, p.second.all_allocs, p.second.alloc_size);
    }
    for (auto &p : large_allocs)
    {
        logerr("still in use (%uud %d)",p.first, p.second);
    }
    //MemoryManager::PrintInfo();
    Countable::get_cur_count();
}
int full_count = 0;
Countable::Countable() {count++;full_count++;};
Countable::~Countable() {count--;};
void Countable::get_cur_count()  {logerr("count = %d %d",count, full_count);}
long Countable::count = 0;