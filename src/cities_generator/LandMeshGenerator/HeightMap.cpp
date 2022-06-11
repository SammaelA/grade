#include "HeightMap.h"
#include "third_party/stb_image.h"
#include "cities_generator/ImageLoader.h"
#include "common_utils/perlin.h"

HeightMap::HeightMap(int width, int height) :
    DefaultMap(width, height)
{
    m_isValidBorders = false;
}

HeightMap::~HeightMap() {}

float HeightMap::GetSafe_EdgeStrectched(int x, int y)
{
    if (x >= -m_width && x <= m_width && y >= -m_height && y <= m_height)
    {
        return Get(x,y);
    }
    if (x < -m_width)
        return (GetSafe_EdgeStrectched(-m_width, y));
    if (y < -m_height)
        return (GetSafe_EdgeStrectched(x, -m_height));
    if (x > m_width)
        return (GetSafe_EdgeStrectched(m_width, y));
    if (y > m_height)
        return (GetSafe_EdgeStrectched(x, m_height));
}

float HeightMap::GetSafe_EdgeStrectched(vec2Int c)
{
    return GetSafe_EdgeStrectched(c.x, c.y);
}

float HeightMap::GetSafe_EdgeRepeated(int x, int y)
{
    if (x >= -m_width && x <= m_width && y >= -m_height && y <= m_height)
    {
        return Get(x,y);
    }
    if (x < -m_width)
        return (GetSafe_EdgeRepeated(2 * m_width + x, y));
    if (y < -m_height)
        return (GetSafe_EdgeRepeated(x, 2 * m_height + y));
    if (x > m_width)
        return (GetSafe_EdgeRepeated(x - 2 * m_width, y));
    if (y > m_height)
        return (GetSafe_EdgeRepeated(x, y - 2 * m_height));
}

float HeightMap::GetSafe_EdgeMirrored(vec2Int v)
{
    vec2Int correctCoords = v;
    
    correctCoords.x = modulo(correctCoords.x + m_width, 2 * m_width) - m_width;
    if (modulo(int(floor((float)(v.x + m_width) / (2 * m_width))), 2) == 1)
        correctCoords.x = -correctCoords.x;
    
    correctCoords.y = modulo(correctCoords.y + m_height, 2 * m_height) - m_height;
    if (modulo(int(floor((float)(v.y + m_height) / (2 * m_height))), 2) == 1)
        correctCoords.y = -correctCoords.y;

    return Get(correctCoords.x, correctCoords.y);
}

float HeightMap::GetSafe_EdgeMirrored(int x, int y)
{
    return GetSafe_EdgeMirrored(vec2Int(x,y));
}

void HeightMap::Set(int x, int y, float val)
{
    Set(vec2Int{x, y}, val);
}

void HeightMap::Set(vec2Int pos, float val)
{
    if (pos.x >= -m_width && pos.x <= m_width && pos.y >= -m_height && pos.y <= m_height)
    {
        m_data[(2*m_width + 1)*(pos.y + m_height) + pos.x + m_width] = val;
        UpdateHeightBorders(val);
    }
    else
    {
        debug("OUT OF BOUNDARIES [HeightMap::Set]");
        throw std::exception{};
    }
}

void HeightMap::FillByPerlin(float min, float max, glm::ivec2 shift)
{
    if (!m_data)
        return;
    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            float height = min + 0.5 * (max - min)
                * (1 + f_perlin(
                    (float)(i + m_width+shift.x)/(2*m_width + 1),
                    (float)(j + m_height+shift.y)/(2*m_height + 1)));
            Set(i,j,height);
        }
    }
    SaveAsTexture();
}

void HeightMap::UpdateHeightBorders(float height)
{
    if (!m_isValidBorders)
    {
        m_isValidBorders = true;
        m_maxHeight = m_minHeight = height;
    }
    m_maxHeight = std::max(m_maxHeight, height);
    m_minHeight = std::min(m_minHeight, height);
}

void HeightMap::RecalculateHeightBorders()
{
    m_isValidBorders = false;
    for (int i = 0; i < m_height * 2 + 1; i++)
    {
        for (int j = 0; j < m_width * 2 + 1; j++)
        {
            int coords[] = {i - m_height, j - m_width};
            float height = Get(coords[0], coords[1]);
            if (!m_isValidBorders)
            {
                m_isValidBorders = true;
                m_maxHeight = height;
                m_minHeight = height;
                continue;
            }
            m_maxHeight = std::max(m_maxHeight, height);
            m_minHeight = std::min(m_minHeight, height);
        } 
    }
}

std::array<std::array<int, 2>, 2> HeightMap::GetSurroundingRect(glm::vec2 TC)
{
    std::array<std::array<int, 2>, 2> surroundingRect = 
    {
        std::array<int, 2> {(int) floor(TC.x), (int) floor(TC.y)},
        std::array<int, 2> {(int) floor(TC.x) + 1, (int) floor(TC.y) + 1}
    }; 
    return surroundingRect;
}

HeightMap::HeightMap(const HeightMap& h)
    : DefaultMap(h)
{
    m_minHeight = h.m_minHeight;
    m_maxHeight = h.m_maxHeight;
    m_isValidBorders = h.m_isValidBorders;
}

void HeightMap::ComponentWiseOperation(const HeightMap& other, const char operation)
{
    if (m_data == nullptr || 
        other.m_data == nullptr || 
        m_width != other.m_width || 
        m_height != other.m_height
    )
    {
        debug("CANT OPERATE THIS HEIGHTMAPS!");
        throw std::exception{};
    }

    for (int i = 0; i < m_height * 2 + 1; i++)
    {
        for (int j = 0; j < m_width * 2 + 1; j++)
        {
            int coords[] = {i - m_height, j - m_width};
            float result = 0;
            switch (operation)
            {
                case '+':
                    result = Get(coords[0], coords[1]) + other.Get(coords[0], coords[1]);
                    break;

                case '-':
                    result = Get(coords[0], coords[1]) - other.Get(coords[0], coords[1]);
                    break;

                default:
                    debug("WRONG OPERATION");
                    throw std::exception{};
                
            }
            Set(coords[0], coords[1], result);
        } 
    }
    RecalculateHeightBorders();
}

void HeightMap::operator-=(const HeightMap& other)
{
    ComponentWiseOperation(other, '-');
}

void HeightMap::operator+=(const HeightMap& other)
{
    ComponentWiseOperation(other, '+');
}

unsigned HeightMap::PointToPixel(float height)
{
    unsigned color = (int)std::round((height - m_minHeight)/(m_maxHeight - m_minHeight) * 255);
    return (color << 24) | (color << 16) | (color << 8) | 255;
            
}

glm::vec2 HeightMap::PointIntToNormalized(vec2Int pos, bool safe)
{
    return PointIntContinuousToNormalized(glm::vec2(pos), safe);
}

glm::vec2 HeightMap::PointIntContinuousToNormalized(glm::vec2 pos, bool safe)
{
    glm::vec2 result{
        (float)(m_width + pos.x) / (float) (2 * m_width), 
        (float)(m_height + pos.y) / (float) (2 * m_height)
    };
    // if (result.x > 1.f || result.x < 0.f || result.y > 1.f || result.y < 0.f)
    // {
    //     debug("WRONG NORMALIZATION COORDS! [HeightMap::PointIntContinuousToNormalized]", pos);
    //     throw std::exception{};
    // }
    return result;
}

glm::vec2 HeightMap::PointNormalizedToIntContinuous(glm::vec2 pos, bool safe)
{
    if (pos.x < 0 || pos.y < 0 || pos.x > 1 || pos.y > 1)
    {
        if (!safe)
        {
            debug("BAD INPUT! [HeightMap_PointNormalizedToIntContinuous]", pos.x, pos.y);
            throw std::exception{};
        }
    }
    return glm::vec2 {
        pos.x * (float) (2 * m_width) -  m_width,
        pos.y * (float) (2 * m_height) -  m_height,
    };
}

void HeightMap::BlurHeightmap()
{
    const float RADIUS = 2.2;
    const float GAUSS_SIGMA = 0.7; // less means narrower blur

    const int radiusInt = int(ceilf(RADIUS));
    const float GAUSS_MULTIPLIER_OUT = 1 / (GAUSS_SIGMA * 2.50662827463);
    const float GAUSS_MULTIPLIER_IN = 1 / (GAUSS_SIGMA * GAUSS_SIGMA * 2);

    float weightSum = 0;
    for (int i = -radiusInt; i <= radiusInt; i++)
    {
        for (int j = -radiusInt; j <= radiusInt; j++)
        {
            float len = sqrtf(i * i + j * j);
            if (len > RADIUS)
                continue;
            float weight = GAUSS_MULTIPLIER_OUT * exp(-(len*len) * GAUSS_MULTIPLIER_IN);
            weightSum += weight;
        }
    }
    const float CORRECTION = 1.f / weightSum;

    for (int i = -radiusInt; i <= radiusInt; i++)
    {
        for (int j = -radiusInt; j <= radiusInt; j++)
        {
            float len = sqrtf(i * i + j * j);
            if (len > RADIUS)
                continue;
            float weight = GAUSS_MULTIPLIER_OUT * exp(-(len*len) * GAUSS_MULTIPLIER_IN);
            weight *= CORRECTION;
        }
    }

    
    HeightMap* temp = new HeightMap(*this);
    for (int i = -m_width; i <= m_width; i++)
    {
        for (int j = -m_height; j <= m_height; j++)
        {
            float currentHeight = 0;
            for (int ki = -radiusInt; ki <= radiusInt; ki++)
            {
                for (int pj = -radiusInt; pj <= radiusInt; pj++)
                {
                    float len = sqrtf(ki * ki + pj * pj);
                    if (len > RADIUS)
                        continue;
                    float weight = exp(-(len * len) * GAUSS_MULTIPLIER_IN);
                    weight *= GAUSS_MULTIPLIER_OUT * CORRECTION;
                    currentHeight += weight * temp->GetSafe_EdgeStrectched(vec2Int{i + ki, j + pj});
                } 
            } 
            Set(i, j, currentHeight);
        }
    }
    delete temp;
}

float HeightMap::GetInterpolated_IntContinuous(glm::vec2 _coords)
{
    float resultHeight = 0;
    std::array<std::array<int, 2>, 2> surroundingRect = 
        HeightMap::GetSurroundingRect(_coords);
    for (int i = 0; i < 4; i++)
    {
        int coords[] = {i / 2, i % 2};
        float currentHeight = Get(
            surroundingRect[coords[0]][0], 
            surroundingRect[coords[1]][1]
        );

        float part = (1 - std::abs(surroundingRect[coords[0]][0] - _coords.x)) * 
            (1 - std::abs(surroundingRect[coords[1]][1] - _coords.y));

        resultHeight += part * currentHeight;
    }
    return resultHeight;
}