#pragma once
#include "cities_generator/global.h"
#include "cities_generator/shaderClass.h"

enum RENDER_MODE
{
    DEFUALT = 0,
    CLIPPING = 1,
    BLOOM = 2,
};

class VertexDataHolder
{
    friend class Entity;
    public:
        GLuint VBO, VAO, IBO;
        unsigned int vertexesNum;
        unsigned int indicesNum;
        
        bool IsUsingIndexes();
        unsigned Get_IBO_element_type();
        
        virtual void Render(RENDER_MODE mode);
        VertexDataHolder();
        virtual ~VertexDataHolder();

    protected:
        void Set_IBO_element_type(unsigned);

    private:
        unsigned m_IBO_element_type;
        
};
