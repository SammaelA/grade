#pragma once
#include "cities_generator/global.h"
#include <string>
#include <vector>

class DebugSphere;
class Entity;

class Renderer
{
    class DebugSphere
    {
        friend class Renderer;

        static std::vector<DebugSphere*> spheres; 
        static void Clear();

        DebugSphere(float radius, glm::vec3 col, glm::vec3 pos, std::string);
        ~DebugSphere();

        static DebugSphere* GetSphereByName(std::string, bool safe = false);
        
        glm::vec3 col;
        float radius;
        std::string name;
        Entity* entityPointer;
    };
    
    static Renderer* instance;

    public:
        static Renderer* GetInstance();

        void MyBufferClearSystem(GLfloat r, GLfloat g, GLfloat b);
        void RenderTexture(GLuint tex, GLfloat scaleFrom, GLfloat scaleTo, bool keepRatio);
        void CreateDebugSphere(float r, glm::vec3 pos, glm::vec3 col=glm::vec3{-1}, std::string n = "");
        void DebugSphereMove(glm::vec3 pos, std::string);
        void DebugSphereSetColor(glm::vec3 col, std::string);
        bool DebugSphereDoesExist(std::string name);
        void Destroy();
        void RenderDebugSpheres();
        void SetFullViewport();
};