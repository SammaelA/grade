#include "VertexDataHolder.h"

#include "../LandMeshGenerator/LandGeneratorClass.h"
#include "cities_generator/ECSClass.h"

unsigned VertexDataHolder::Get_IBO_element_type()
{
    if(m_IBO_element_type == 0)
    {
        std::cout << "WRONG IBO ELEMENT TYPE" << std::endl;
        throw std::exception{}; 
    }
    else
    {
        return m_IBO_element_type;
    }
}

void VertexDataHolder::Set_IBO_element_type(unsigned newType)
{
    if (!IsUsingIndexes() && newType != 0)
        glGenBuffers(1, &IBO);
    m_IBO_element_type = newType;
}

VertexDataHolder::VertexDataHolder()
{
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    m_IBO_element_type = 0;
}

VertexDataHolder::~VertexDataHolder()
{
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    if (IsUsingIndexes())
    {
        glDeleteBuffers(1, &IBO);
    }
}

void VertexDataHolder::Render(RENDER_MODE rmode)
{
    glBindVertexArray(VAO);
    if (IsUsingIndexes())
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
        glDrawElements(GL_TRIANGLES, indicesNum, Get_IBO_element_type(), 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
    else
    {
        glDrawArrays(GL_TRIANGLES, 0, vertexesNum);
    }
    glBindVertexArray(0);
}

bool VertexDataHolder::IsUsingIndexes()
{
    return m_IBO_element_type != 0;
}