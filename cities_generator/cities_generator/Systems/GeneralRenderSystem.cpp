#include "cities_generator/ECSClass.h"
#include "cities_generator/global.h"

using namespace glm;

#define SYSTEM_NAME GeneralRenderSystem
#define TYPE1 ModelInfoComponent
#define Q1_TYPE1 TransformComponent
#define Q1_TYPE2 SunLightComponent
#define Q2_TYPE1 CameraComponent


class SYSTEM_NAME : public System
{
    static SYSTEM_NAME* inst;
    SYSTEM_NAME() 
    {
        //components
        activeComponents.push_back(TYPE1::GetCompId());

        //Just finishing creation
        UpdateArchetypes();

        //Necessary queues
        AddQueue("get_sun_queue",
                 std::vector<ComponentType>{Q1_TYPE1::GetCompId(),
                                            Q1_TYPE2::GetCompId()});
        
        AddQueue("get_camera_queue",
                 std::vector<ComponentType>{Q2_TYPE1::GetCompId()});

    }
    virtual void RunSystem(std::vector<Component*> components)
    {
        //Preparations
        TYPE1* model = static_cast<TYPE1*>(components[0]);
        
        auto query1_components = (*GetQueueByName("get_sun_queue")->begin());
        Q1_TYPE1* sunTransform = static_cast<Q1_TYPE1*>(query1_components[0]);
        Q1_TYPE2* sunComp = static_cast<Q1_TYPE2*>(query1_components[1]);

        Q2_TYPE1* cam = SingleCompFindQueue<Q2_TYPE1>("get_camera_queue");

        //System body
        
    }
    virtual ~SYSTEM_NAME() = default;
};

SYSTEM_NAME* SYSTEM_NAME::inst = new SYSTEM_NAME();
