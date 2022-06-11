#pragma once
#include <GL/glew.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <GLFW/glfw3.h>
#include <unordered_set>
#include <vector>
#define GLM_FORCE_RADIANS 1
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtx/quaternion.hpp"
#include "glm/gtx/euler_angles.hpp"

#define COMPONENT_COUNTER 16
#define BAD_TEXTURE -1
#define DEFAULT_COLOR glm::vec4(0.0, 0.5, 0.0, 1.0)
#define ASSETS_FOLDER "resources/"
#define SHADERS_FOLDER "shaders/cities_gen/"

#define IS_DERIVED(pointer, class) (dynamic_cast<class*>(pointer) != nullptr)


using ComponentType = unsigned int;

class vec2Int
{
    public:
        int x, y;
        vec2Int() = default;
        vec2Int(int _x, int _y) : x{_x}, y{_y} {} 
        explicit vec2Int(const glm::vec2& v);
        int& operator[](int index);
        vec2Int operator+(const vec2Int& v) const;
        vec2Int operator-(const vec2Int& v) const;
        float length() const;
        operator glm::vec2() const;
};

bool operator==(const vec2Int& v1, const vec2Int& v2);

template<>
struct std::hash<vec2Int>
{
    size_t operator() (const vec2Int& v) const
    {
        return v.x * 10000 + v.y;
    }
};

template<class T>
T LERP(T a, T b, float t)
{
    return (a + t * (b - a));
}

template<class T>
T clamp(T a, T _min, T _max)
{
    if (_min > _max)
    {
        std::cout << "WRONG CLAMP ARGUMENTS!" << std::endl;
        throw std::exception{};
    }
    a = (a < _min) ? _min : a;
    a = (a > _max) ? _max : a;
    return a;
}

template<class T>
int sign(T x)
{
    return (x > 0) - (x < 0);
}

template<class T>
int sign_no_zero(T x)
{
    return (x >= 0) - (x < 0);
}

void intersect(std::unordered_set<unsigned>&, const std::unordered_set<unsigned>&);

template<class T>
void debug_part(T s)
{
    std::cout << s << " ";
}

template<>
void debug_part<glm::vec3>(glm::vec3 v);

template<>
void debug_part<glm::vec2>(glm::vec2 v);


template<>
void debug_part<glm::vec4>(glm::vec4 v);

template<>
void debug_part<vec2Int>(vec2Int v);

template<class First, class... Rest>
void debug_part(First f, Rest... r)
{
    debug_part(f);
    debug_part(r...);
}

template<class T>
void debug(T s)
{
    debug_part(s);
    std::cout << std::endl;
}

template<class First, class... Rest>
void debug(First f, Rest... r)
{
    debug_part(f);
    debug(r...);
}


void checkForGlErrors(std::string text, bool clean = false);

std::vector<std::string> Split(std::string);

float rndNormalized();

glm::vec2 LeftNormal(glm::vec2 v);

glm::vec2 XZ(glm::vec3 v);
glm::vec2 XY(glm::vec3 v);
glm::vec3 X0Z(glm::vec2 v);
glm::vec3 XY0(glm::vec2 v);
glm::vec3 _0Y0(float y);

template<class T>
bool Approx(T a, T b, float DELTA = 0.001)
{
    return glm::length(a - b) < DELTA;
}

int modulo(int x, int y);

float rnd();

float SafeAcosf(float cos);

glm::vec2 NearestPointOfSection(glm::vec2, glm::vec2, glm::vec2);

float DistanceToSection(glm::vec2, glm::vec2, glm::vec2);

std::pair<glm::vec2, bool> IntersectSegments(glm::vec2, glm::vec2, glm::vec2, glm::vec2, float = 0.005);