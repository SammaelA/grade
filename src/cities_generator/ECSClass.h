#pragma once
#include <GL/glew.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <array>
#include <unordered_set>
#include <functional>
#include <iostream>
#include <map>

#include "./shaderClass.h"
#include "cities_generator/global.h"

#ifndef GLM_FORCE_RADIANS
#define GLM_FORCE_RADIANS 1
#endif
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "ShaderDataHolders/Archetypes/BasicGLTFModel.h"


class Component
{
    public:
        ComponentType type; 
        static std::array<std::unordered_set<unsigned>, COMPONENT_COUNTER> componentOwnersArray;
        Component(unsigned, unsigned);
        Component() = default;
        virtual ~Component() = default;

        template <typename TComponent>
        static TComponent* AnyComponentTyped(bool=false);
};

class Archetype;

class ArchetypeGraphEdge
{
    public:
        Archetype* mate;
        ComponentType componentLink;        
};

class Archetype
{
    public:
        std::vector<ComponentType> type;
        std::vector<std::vector<Component*>> components;
        std::vector<unsigned int> entitiesId; 
        std::vector<ArchetypeGraphEdge> mates;
        void AddComponentToEntity(unsigned int EntityIndex, Component *component);
        void debugPrintType();
        void debugShowProducedTree();

        static Archetype* basicNode;
        static int AddEmptyEntity(unsigned id);
        static void debugShowTree();
        static void handleNewArchetypes();
        static bool hasUnhandledArchetypes;

    private:
        Archetype(std::vector<ComponentType> type);
};

class Record {
    public:
        Archetype *archetype;
        unsigned int index;
        Record(Archetype *a, unsigned int i) : archetype(a), index(i) {}
        Record() = default;
};

class Entity
{
    public:
        unsigned int id;
        // std::vector<Component*> components;

        static std::vector<Entity*> entities;
        static std::unordered_map<unsigned int, Entity*> entities_map;
        static unsigned int currentId;
        static Entity* EntityById(unsigned,bool=false);
        static Entity* AnyObjectWithComponent(unsigned,bool=false);
        static Entity* AnyObjectWithTag(std::string,bool=false);
        static std::vector<Entity*> AllEntitiesWithComponentsAndTags(unsigned n_components, 
            unsigned n_tags, const std::vector<unsigned>&, const  std::vector<std::string>&, bool=false);
        static std::unordered_map<unsigned int, Record> entitiesArchetype;

        Entity();
        bool HasComponent(unsigned);
        void AddComponent(Component *comp);
        
        template <typename TComponent>
        TComponent* GetComponent(bool=false);

        void giveTag(std::string);
        bool hasTag(std::string);

        void LoadGltfModel(std::string path, std::string obj_name, std::string shaderName);

        ~Entity();
};

class Queue;

class System
{
    public:
        std::vector<ComponentType> activeComponents;
        std::unordered_set<Archetype*> relevantArchetypes;
        std::unordered_set<Queue*> queues;
        bool isSeparateSystem;

        static std::vector<System*> allSystems;
        
        System();
        void UpdateArchetypes();
        virtual void RunSystem(std::vector<Component*>) = 0;
        void AddQueue(std::string, std::vector<ComponentType>);
        Queue* GetQueueByName(std::string);

        template <typename TComponent>
        TComponent* SingleCompFindQueue(std::string);

        virtual ~System() {};
    
        static void UpdateAllSystems(); 
        static void RunAllSystems();
};

class Queue : public System
{
    friend class System;
    
    std::string name;
    virtual void RunSystem(std::vector<Component*> c);

    Queue(std::string);

    public:
        struct Iterator
        {
            using iterator_category = std::input_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = std::vector<Component*>;
            using pointer           = std::vector<Component*>*;
            using reference         = std::vector<Component*>&;
            using const_reference   = const std::vector<Component*>&;

            void UpdatePositionsInArchetype();
            bool ValidateIndexes();
            void UpdateCurrentComponents();
            bool IsEndIterator();

            Iterator(const Iterator&) = default; // нужен итератор в begin
            Iterator(Queue*, bool);
            ~Iterator() = default;

            const_reference operator*() const; 

            // Prefix increment
            Iterator& operator++();  

            // Postfix increment
            Iterator operator++(int);

            friend bool operator== (const Iterator& a, const Iterator& b);
            friend bool operator!= (const Iterator& a, const Iterator& b);

            private:
                value_type m_currentComponents;
                std::unordered_set<Archetype*>::iterator m_archetypesIter;
                unsigned m_entityInArchetypeNum;
                Queue* m_parentQueue;
                std::vector<unsigned> m_componentPositionsInArchetype;
        };

        Iterator begin();
        Iterator end();
};

class SpriteComponent : public Component
{
    public:
        float pos[2];
        float size[2];
        GLuint VBO, VAO, EBO;
        Shader* spriteShader;
        GLuint spriteTexture;
        static unsigned int typeId;
        
        SpriteComponent(float,float,float,float,std::string, unsigned);
        SpriteComponent(float,float,float,float,std::string,std::string, unsigned);
        ~SpriteComponent();
        static unsigned GetCompId() {return typeId;}
};

class ModelInfoComponent : 
    public Component, 
    public BasicGLTFModel
{
    public:
        static unsigned int typeId;
        
        // ModelInfoComponent(std::string shaderName, unsigned EID);
        ModelInfoComponent(unsigned EID);
        virtual ~ModelInfoComponent();
        static unsigned GetCompId() {return typeId;}
};

class ButtonComponent : public Component
{
    public:
        float pos[2], size[2];
        int layer;
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        std::function<void(float, float, int)> func;
        ButtonComponent(int, float, float, float, float, std::function<void(float, float, int)>, unsigned);
};

class ShaderLibrary : public Component
{
    public:
        std::vector<Shader*> shaderVector;
        static unsigned int typeId;
        static ShaderLibrary* library;
        static unsigned GetCompId() {return typeId;}
        Shader* GetShaderPointer(std::string);
        ShaderLibrary(unsigned EID) : Component(EID, ShaderLibrary::typeId) {};
        static ShaderLibrary* GetLibrary();
};

class TransformComponent : public Component
{
    public:
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        glm::vec3 position;
        glm::vec3 scale;
        glm::mat4 rotation;
        glm::mat4 GetMatrix();
        glm::mat4 GetViewMatrixFromThisPoint();
        TransformComponent(unsigned);
        TransformComponent(float, float, float, unsigned);
        TransformComponent(const TransformComponent&, unsigned);
        static glm::vec3 GetDefaultViewDirection() {return glm::vec3(1,0,0);};
};

class CameraComponent : public Component
{
    public:
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        CameraComponent(unsigned EID) : Component(EID, CameraComponent::typeId) {}
        glm::mat4 GetViewProjectionMatrix();
        glm::mat4 GetProjectionMatrix();
        void SetProjectionMatrix(glm::mat4 projM);
        glm::mat4 GetInvProjMatrix();
        glm::mat4 GetViewMatrix();
        void MirrorFromWater(float waterLevel);
        void ConsoleHandler();
        glm::vec3 GetPosition();
        void HandleInput();
        glm::vec2 GetZDists();

    private:
        glm::vec3 position, direction;
        glm::mat4 projectionMatrix, invprojectionMatrix, viewMatrix;
        float zNear, zFar;

    friend Entity* createCameraTemplate(float, float, float, float, glm::vec3, glm::vec3, std::string);
};

class CameraControlComponent : public Component
{
    public:
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        CameraControlComponent(unsigned);
        std::string state;
        int switchStateKey;

        float freeRotateSpeed[2];
        float freeVerticalLimit;
        float freeMovementSpeed;

        glm::vec3 satelliteCenter;
        float satelliteRotateSpeed[2];
        float satelliteVerticalLimit;
        float satelliteRadius;
        float satelliteLatitudeDegrees;
        float satelliteLongitudeDegrees;
        float satelliteZoomSpeed;
        float satelliteZoomBorders[2];
};

class TagComponent : public Component
{
    friend void Entity::giveTag(std::string);
    public:
        std::vector<std::string> tags;
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        static std::unordered_map<std::string,std::unordered_set<unsigned>> tagOwnersMap;
    private:
        TagComponent(unsigned);
};

class SunLightComponent : public Component
{
    public:
        glm::vec3 color;
        glm::mat4 projectionMatrix;
        unsigned int SHADOW_WIDTH, SHADOW_HEIGHT;
        unsigned int depthMapFBO;
        GLuint depthTexture;
        GLuint shadowTexture;

        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        SunLightComponent(glm::vec3 col, unsigned w, unsigned h, float scale, unsigned ID);
};

class ImageVertexInfoComponent : public Component
{
    public:
        GLuint VBO, VAO, EBO;
        static unsigned int typeId;
        Shader* shader;

        static unsigned GetCompId() {return typeId;}
        ImageVertexInfoComponent(unsigned EID);
        ~ImageVertexInfoComponent();
};

class CompA : public Component
{
    public:
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        CompA(unsigned EID)
            :Component(EID, CompA::GetCompId()) {}
};

class CompB : public Component
{
    public:
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        CompB(unsigned EID)
            :Component(EID, CompB::GetCompId()) {}
};

class CompC : public Component
{
    public:
        static unsigned int typeId;
        static unsigned GetCompId() {return typeId;}
        CompC(unsigned EID)
            :Component(EID, CompC::GetCompId()) {}
};

Entity* createImageTemplate(float, float, float, float, std::string);
Entity* createShaderLibraryTemplate();
Entity* createCameraTemplate(float,float,float,float,glm::vec3,glm::vec3,std::string);



template <typename TComponent>
TComponent* Component::AnyComponentTyped(bool is_optional)
{
    auto iterator = Component::componentOwnersArray[TComponent::typeId].begin();
    while (iterator != Component::componentOwnersArray[TComponent::typeId].end())
    {
        Entity* entity = Entity::EntityById(*iterator, true);
        if (entity == nullptr)
        {
            iterator = Component::componentOwnersArray[TComponent::typeId].erase(iterator);
        }
        else
        {
            return entity->GetComponent<TComponent>();
        }
    }
    if (is_optional)
    {
        return nullptr;
    }
    else
    {
        std::cout<<"----------------NO SUCH COMPONENTS [AnyComponentTyped]-----------------\n";
        throw std::exception{};
    }
}

template <typename TComponent>
TComponent* Entity::GetComponent(bool is_optional)
{
    Record rec = Entity::entitiesArchetype[id];
    auto archetypeType = rec.archetype->type;
    auto iter = std::find(archetypeType.begin(), archetypeType.end(), TComponent::typeId);

    if (iter != archetypeType.end())
    {
        return static_cast<TComponent*>(
            rec.archetype->components[iter - archetypeType.begin()][rec.index]);
    }
    if (is_optional)
    {
        return nullptr;
    }
    else
    {
        std::cout<<"----------------NO SUCH COMPONENT [GetComponent]-----------------\n";
        throw std::exception{};
    } 
}

class SleepInterval
{
    private:
        static std::map<std::string, SleepInterval*> intervals;
        bool wait(float milisec);
        float counter;
        SleepInterval() : counter{0} {};
    public:
        // return true once in "waitMilisec" miliseconds
        static bool interval(float waitMilisec, std::string name); 
};

#include "ECSClassTemplates.tpp"