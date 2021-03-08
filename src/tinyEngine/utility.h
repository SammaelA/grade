#pragma once
#include <cstdio>
#include <glm/glm.hpp>
#include <stdarg.h>
#define PI 3.14159265f  
#define MIN(a,b) (a<b ? a : b)
#define MAX(a,b) (a>b ? a : b)
#define CLAMP(a,b,c) (MIN(MAX(a,b),c))
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