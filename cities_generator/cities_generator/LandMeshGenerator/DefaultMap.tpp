#include <algorithm>

template<class T>
DefaultMap<T>::DefaultMap(int width, int height)
{
    if (width <= 0 || height <= 0)
    {
        debug("INCORRECT SIZE");
        throw std::exception{};
    }
    m_width = width;
    m_height = height;
    dataSize = (2*m_width + 1)*(2*m_height + 1) * sizeof(T);
    m_data = (T*) malloc(dataSize);
    m_textureId = BAD_TEXTURE;
}

template<class T>
DefaultMap<T>::~DefaultMap()
{
    free(m_data);
    if (m_textureId != BAD_TEXTURE)
        glDeleteTextures(1, &m_textureId);
}

template<class T>
DefaultMap<T>::DefaultMap(const DefaultMap& h)
{
    m_width = h.m_width;
    m_height = h.m_height;

    dataSize = h.dataSize;
    m_data = (T*) malloc(dataSize);
    std::memcpy(m_data, h.m_data, dataSize);
    SaveAsTexture();
}

template<class T>
T DefaultMap<T>::Get(int x, int y, bool safe, T baseValue) const
{
    if (x >= -m_width && x <= m_width && y >= -m_height && y <= m_height)
    {
        return m_data[(2*m_width + 1)*(y + m_height) + x + m_width];
    }
    if (safe)
    {
        return baseValue;
    }
    debug("OUT OF BOUNDARIES [DefaultMap::Get]");
    throw std::exception{};
}

template<class T>
T DefaultMap<T>::Get(vec2Int pos, bool safe, T baseValue) const
{
    return Get(pos.x, pos.y, safe, baseValue);
}

template<class T>
void DefaultMap<T>::Set(int x, int y, T val)
{
    if (x >= -m_width && x <= m_width && y >= -m_height && y <= m_height)
    {
        m_data[(2*m_width + 1)*(y + m_height) + x + m_width] = val;
    }
    else
    {
        debug("OUT OF BOUNDARIES [DefaultMap::Set]");
        throw std::exception{};
    }
}

template<class T>
void DefaultMap<T>::Set(vec2Int pos, T val)
{
    Set(pos.x, pos.y, val);
}

template<class T>
void DefaultMap<T>::SaveAsTexture()
{
    if (m_data == nullptr || m_width == 0 || m_height == 0)
    {
        debug("CANT SAVE TEXTURE");
        throw std::exception{};
    }

    int TEX_WIDTH = 1024;
    int TEX_HEIGHT = 1024;
    bool isScaled = true;

    if (m_width < 513 && m_width < 513)
    {
        isScaled = false;
        TEX_WIDTH = m_width * 2 + 1;
        TEX_HEIGHT = m_height * 2 + 1;
    }
    
    GLubyte* image = (GLubyte*)malloc(TEX_WIDTH * TEX_HEIGHT * 4 * sizeof(GLubyte));
    for (unsigned j = 0; j < TEX_HEIGHT; j++) {
        for (unsigned i = 0; i < TEX_WIDTH; i++) {
            T value;
            if (isScaled)
            {
                value = Get(
                    (int)((float) i / (TEX_WIDTH - 1) * 2 * m_height) - m_width,
                    (int)((float) j / (TEX_HEIGHT - 1) * 2 * m_width) - m_height
                );
            }
            else
            {
                value = Get(i - m_width, j - m_height);
            }
            unsigned color = PointToPixel(value);
            image[0 + (i + j*(TEX_WIDTH))*4] = (GLubyte) (color >> 24) ;
            image[1 + (i + j*(TEX_WIDTH))*4] = (GLubyte) (color << 8 >> 24);
            image[2 + (i + j*(TEX_WIDTH))*4] = (GLubyte) (color << 16 >> 24);
            image[3 + (i + j*(TEX_WIDTH))*4] = (GLubyte) (color << 24 >> 24);
        }
    }
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &m_textureId);
    glBindTexture(GL_TEXTURE_2D, m_textureId);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, 
                    GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
                    GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEX_WIDTH, 
                TEX_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, 
                image);
    stbi_image_free(image);
    glBindTexture(GL_TEXTURE_2D, 0);
}

template<class T>
bool DefaultMap<T>::IsTextureReady()
{
    return m_textureId != BAD_TEXTURE;
}

template<class T>
bool DefaultMap<T>::IsPointValid(int x, int y)
{
    return (x >= -m_width && x <= m_width && y >= -m_height && y <= m_height);
}

template<class T>
bool DefaultMap<T>::IsPointValid_Normalized(glm::vec2 v)
{
    return (v.x >= 0.f && v.x <= 1.f && v.y >= 0.f && v.y <= 1.f);
}

template<class T>
bool DefaultMap<T>::IsPointValid_IntContinuous(glm::vec2 v)
{
    float x = v.x;
    float y = v.y;
    return (x >= -m_width && x <= m_width && y >= -m_height && y <= m_height);
}

template<class T>
bool DefaultMap<T>::IsPointValid(vec2Int v) 
{
    return IsPointValid(v.x, v.y);
}

template<class T>
GLuint DefaultMap<T>::GetTexture()
{
    if (IsTextureReady())
        return m_textureId;
    else
    {
        debug("TEXTURE IS NOT READY");
        throw std::exception{};
    }
}

template<class T>
glm::vec2 DefaultMap<T>::GetSize()
{
    return glm::vec2{m_width, m_height};
}

template<class T>
unsigned DefaultMap<T>::PointToPixel(T val)
{
    return (255u << 24) | (255u << 8) | 255u;
}

template<class T>
unsigned DefaultMap<T>::ColorLerp(unsigned col1, unsigned col2, float t)
{
    unsigned res = 0;
    for (int i = 0; i < 4; i++)
    {
        res |= (unsigned) LERP((float)((col1 << 8 * i) >> 24), (float)((col2 << 8 * i) >> 24), t) << 8 * (3 - i); 
    }
    return res;
}
