#include "cities_generator/global.h"
#include <vector>
#include <sstream>
#include <algorithm>

template<>
void debug_part<glm::vec3>(glm::vec3 v)
{
    debug_part(v.x, v.y, v.z);
}

template<>
void debug_part<glm::vec2>(glm::vec2 v)
{
    debug_part(v.x, v.y);
}

template<>
void debug_part<vec2Int>(vec2Int v)
{
    debug_part(v.x, v.y);
}

int& vec2Int::operator[](int index)
{
    if (index < 0 || index > 1)
    {
        debug("WROND INDEX! [vec2Int::operator[]]");
        throw std::exception{};
    }
    if (index == 0)
    {
        return x;
    }
    if (index == 1)
    {
        return y;
    }
    return x;
}

vec2Int vec2Int::operator+(const vec2Int& v) const
{
    vec2Int res = *this;
    res.x += v.x;
    res.y += v.y;
    return res;
}

vec2Int vec2Int::operator-(const vec2Int& v) const
{
    vec2Int res = *this;
    res.x -= v.x;
    res.y -= v.y;
    return res;
}

float vec2Int::length() const
{
    return sqrtf(x*x + y*y);
}

bool operator==(const vec2Int& v1, const vec2Int& v2)
{
    return (v1.x == v2.x && v1.y == v2.y);
}

vec2Int::operator glm::vec2() const
{
    return glm::vec2(x, y);
}

template<>
void debug_part<glm::vec4>(glm::vec4 v)
{
    debug_part(v.x, v.y, v.z, v.w);
}

void checkForGlErrors(std::string text, bool clean)
{
    bool noErrorsFlag = true;
    GLenum err;
    while((err = glGetError()) != GL_NO_ERROR)
    {
        debug("ERROR:", text, err);
        noErrorsFlag = false;
    }
    if (!noErrorsFlag && !clean)
        throw std::exception{};
}

std::vector<std::string> Split(std::string s)
{
    std::string SPACE_STRING = " ";
    if (s.size() > 0)
    {
        if (s[s.size() - 1] != SPACE_STRING[0])
            s += SPACE_STRING;
    }

    std::vector<std::string> result;
    size_t pos = 0;
    while ((pos = s.find(SPACE_STRING)) != std::string::npos) {
        result.push_back(s.substr(0, pos));
        s.erase(0, pos + SPACE_STRING.length());
    }
    return result;
}

vec2Int::vec2Int(const glm::vec2& v)
{
    x = int(v.x);
    y = int(v.y);
}

float rndNormalized()
{
    return (float)random() / RAND_MAX;
}

glm::vec2 LeftNormal(glm::vec2 v)
{
    return glm::normalize(glm::vec2{-v.y, v.x});
}

glm::vec2 XZ(glm::vec3 v)
{
    return glm::vec2{v.x, v.z};
}

glm::vec3 X0Z(glm::vec2 v)
{
    return glm::vec3{v.x, 0, v.y};
}

glm::vec2 XY(glm::vec3 v)
{
    return glm::vec2{v.x, v.y};
}

glm::vec3 XY0(glm::vec2 v)
{
    return glm::vec3{v.x, v.y, 0};
}

int modulo(int x, int y)
{
    if (y > 0)
    {
        if (x >= 0)
            return x % y;
        else 
            return (y - ((-x) % y)) % y;
    }
    else
    {
        debug("MODULO ACTION UNDEFINED!");
        throw std::exception{};
    }

}

float rnd()
{
    return (float)rand() / RAND_MAX;
}

float SafeAcosf(float cos)
{
    return acosf(clamp(cos, -1.f, 1.f));
}

glm::vec2 NearestPointOfSection(glm::vec2 point, glm::vec2 s1, glm::vec2 s2)
{
    glm::vec2 pointToS1 = glm::normalize(s1 - point);   
    glm::vec2 pointToS2 = glm::normalize(s2 - point);
    glm::vec2 s1ToS2 = glm::normalize(s2 - s1);

    float s1Cos = glm::dot(-pointToS1, s1ToS2);
    float s2Cos = glm::dot(-pointToS2, -s1ToS2);

    if (s1Cos * s2Cos > 0)
    {
        return s1Cos * glm::length(point - s1) * s1ToS2 + s1;
    }
    else
    {
        return (s1Cos > 0) ? s2 : s1;
    }
}


float DistanceToSection(glm::vec2 point, glm::vec2 s1, glm::vec2 s2)
{
    glm::vec2 pointToS1 = glm::normalize(s1 - point);   
    glm::vec2 pointToS2 = glm::normalize(s2 - point);
    glm::vec2 s1ToS2 = glm::normalize(s2 - s1);

    float s1Cos = glm::dot(-pointToS1, s1ToS2);
    float s2Cos = glm::dot(-pointToS2, -s1ToS2);

    if (s1Cos * s2Cos > 0)
    {
        return std::sqrt(fabsf(1 - s1Cos * s1Cos)) * glm::length(point - s1);
    }
    else
    {
        return (s1Cos > 0) ? glm::length(point - s2) : glm::length(point - s1);
    }
}

std::pair<glm::vec2, bool> IntersectSegments(glm::vec2 p1_1, glm::vec2 p1_2, glm::vec2 p2_1, glm::vec2 p2_2, float DELTA)
{
    float denominator = (p1_1.x - p1_2.x) * (p2_1.y - p2_2.y) - (p1_1.y - p1_2.y) * (p2_1.x - p2_2.x);
    if (Approx(denominator, 0.f))
        return std::pair<glm::vec2, bool>{glm::vec2{0, 0}, false};
    float t = (p1_1.x - p2_1.x) * (p2_1.y - p2_2.y) - (p1_1.y - p2_1.y) * (p2_1.x - p2_2.x);
    t /= denominator;
    float u = (p1_1.x - p2_1.x) * (p1_1.y - p1_2.y) - (p1_1.y - p2_1.y) * (p1_1.x - p1_2.x);
    u /= denominator;
    if (u <= 1.f && u >= 0.f && t <= 1.f && t >= 0.f)
    {
        glm::vec2 result{p1_1.x + t * (p1_2.x - p1_1.x), p1_1.y + t * (p1_2.y - p1_1.y)};
        return std::pair<glm::vec2, bool>{result, true};
    }
    else
    {
        if (DistanceToSection(p1_1, p2_1, p2_2) < DELTA)
            return std::pair<glm::vec2, bool>{p1_1, true};
        if (DistanceToSection(p1_2, p2_1, p2_2) < DELTA)
            return std::pair<glm::vec2, bool>{p1_2, true};
        if (DistanceToSection(p2_1, p1_1, p1_2) < DELTA)
            return std::pair<glm::vec2, bool>{p2_1, true};
        if (DistanceToSection(p2_2, p1_1, p1_2) < DELTA)
            return std::pair<glm::vec2, bool>{p2_2, true};
        return std::pair<glm::vec2, bool>{glm::vec2{0, 0}, false};
    }
}

glm::vec3 _0Y0(float y)
{
    return glm::vec3{0, y, 0};
}