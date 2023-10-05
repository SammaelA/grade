#include "cities_generator/ECSClass.h"
#include "cities_generator/global.h"
#include "./ExternFunctions/externFunctions.h"
#include <iostream>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include "./EventClasses.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "third_party/stb_image.h"
#include "GameManager.h"
#include "InputManager.h"
#include "cities_generator/Renderer.h"

ComponentType SpriteComponent::typeId = 1;
ComponentType TagComponent::typeId = 4;
ComponentType ModelInfoComponent::typeId = 5;
ComponentType ButtonComponent::typeId = 6;
ComponentType ShaderLibrary::typeId = 7;
ComponentType TransformComponent::typeId = 8;
ComponentType CameraComponent::typeId = 9;
ComponentType CameraControlComponent::typeId = 10;
ComponentType SunLightComponent::typeId = 11;
ComponentType ImageVertexInfoComponent::typeId = 12;
ComponentType CompA::typeId = 13;
ComponentType CompB::typeId = 14;
ComponentType CompC::typeId = 15;
//НЕ ЗАБЫВАТЬ УВЕЛИЧИВАТЬ СЧЁТЧИК КОМПОНЕНТОВ!


SpriteComponent::SpriteComponent(float x, float y, float w, float h, std::string path, unsigned EID)
    :SpriteComponent(x, y, w, h, path, "spriteShader", EID) {}

SpriteComponent::SpriteComponent(float x, float y, float w, float h, std::string path, std::string name, unsigned EID)
    : Component(EID, SpriteComponent::typeId)
{
    spriteShader = ShaderLibrary::GetLibrary()->GetShaderPointer(name);

    pos[0] = x;
    pos[1] = y;
    size[0] = w;
    size[1] = h;

    int imageW, imageH;
    unsigned char* image = stbi_load(path.c_str(), &imageW, &imageH, 0, 4);


    glGenTextures(1, &spriteTexture);
    glBindTexture(GL_TEXTURE_2D, spriteTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    glGenerateMipmap(GL_TEXTURE_2D);
    stbi_image_free(image);
    glBindTexture(GL_TEXTURE_2D, 0);

    GLfloat vertexes[] = {
     1.0f,  1.0f, 0.0f,   1.0f, 1.0f,   
     1.0f,  0.0f, 0.0f,   1.0f, 0.0f,   
     0.0f,  0.0f, 0.0f,   0.0f, 0.0f,   
     0.0f,  1.0f, 0.0f,   0.0f, 1.0f    
    };
    GLuint indices[] = {
        0, 1, 3,
        1, 2, 3  
    };
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);// Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertexes), vertexes, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, 0); // Note that this is allowed, the call to glVertexAttribPointer registered VBO as the currently bound vertex buffer object so afterwards we can safely unbind
    glBindVertexArray(0); // Unbind VAO (it's always a good thing to unbind any buffer/array to prevent strange bugs), remember: do NOT unbind the EBO, keep it bound to this VAO
}

SpriteComponent::~SpriteComponent()
{
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

TagComponent::TagComponent(unsigned EID)
    : Component(EID, TagComponent::typeId) {};

ButtonComponent::ButtonComponent(int l, float x, float y, float w, float h, std::function<void(float, float, int)> f, unsigned EID)
    : Component(EID, ButtonComponent::typeId), pos{x,y}, size{w,h}
{
    func = f;
    layer = l;
}

Entity* createImageTemplate(float x, float y, float width, float height, std::string imagePath)
{
    auto im = new Entity();
    // im->components.push_back(new SpriteComponent(x,y,width,height,imagePath, im->id));
    im->AddComponent(new SpriteComponent(x,y,width,height,imagePath, im->id));
    return im;
}

// ModelInfoComponent::ModelInfoComponent(std::string name, unsigned EID)
//     : Component(EID, ModelInfoComponent::typeId)
// {
//     shader = ShaderLibrary::GetLibrary()->GetShaderPointer(name);
//     GLfloat vertexes[] = {
//         -0.5f, -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
//          0.5f, -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
//          0.5f,  0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
//          0.5f,  0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
//         -0.5f,  0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
//         -0.5f, -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        
//         -0.5f, -0.5f,  0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 
//          0.5f, -0.5f,  0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
//          0.5f,  0.5f,  0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
//          0.5f,  0.5f,  0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
//         -0.5f,  0.5f,  0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
//         -0.5f, -0.5f,  0.5f, 0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
        
//         -0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.5f, -1.0f, 0.0f, 0.0f,
//         -0.5f,  0.5f, -0.5f, 0.0f, 0.0f, 0.5f, -1.0f, 0.0f, 0.0f,
//         -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.5f, -1.0f, 0.0f, 0.0f,
//         -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.5f, -1.0f, 0.0f, 0.0f,
//         -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.5f, -1.0f, 0.0f, 0.0f,
//         -0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.5f, -1.0f, 0.0f, 0.0f,
        
//          0.5f,  0.5f,  0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
//          0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
//          0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
//          0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
//          0.5f, -0.5f,  0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
//          0.5f,  0.5f,  0.5f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        
//         -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
//          0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
//          0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
//          0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
//         -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
//         -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        
//         -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
//          0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
//          0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
//          0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
//         -0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
//         -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
//     };

//     vertexesNum = 36;
//     IsUsingIndexes() = false;
//     IBO_element_type = 0;
//     diffuseColor = DEFAULT_COLOR;

//     glGenVertexArrays(1, &VAO);
//     glGenBuffers(1, &VBO);
//     glBindVertexArray(VAO);
//     glBindBuffer(GL_ARRAY_BUFFER, VBO);
//     glBufferData(GL_ARRAY_BUFFER, sizeof(vertexes), vertexes, GL_STATIC_DRAW);
//     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)0);
//     glEnableVertexAttribArray(0);
//     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
//     glEnableVertexAttribArray(1);
//     glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)(6 * sizeof(GLfloat)));
//     glEnableVertexAttribArray(2);
//     glBindBuffer(GL_ARRAY_BUFFER, 0);
//     glBindVertexArray(0);
// }

ModelInfoComponent::ModelInfoComponent(unsigned EID) :
    Component(EID, ModelInfoComponent::typeId){}

Entity* createShaderLibraryTemplate()
{
    auto ent = new Entity();
    ShaderLibrary* libra = new ShaderLibrary(ent->id);
    // ent->components.push_back(libra);
    ent->AddComponent(libra);
    const char* VSPathes[] = {"shaders/cities_gen/spriteVertex.vs", 
                              "shaders/cities_gen/lowPolyVertex.vs", 
                              "shaders/cities_gen/screenDisplay.vs",
                              "shaders/cities_gen/test3dVertex.vs",
                              "shaders/cities_gen/landscapeVertex.vs",
                              "shaders/cities_gen/debugSphereVertex.vs",
                              "shaders/cities_gen/waterPlateVertex.vs",
                              "shaders/cities_gen/skyboxVertex.vs",
                              "shaders/cities_gen/roadVertex.vs",
                              "shaders/cities_gen/buildingVertex.vs"};

    const char* FSPathes[] = {"shaders/cities_gen/spriteFragment.fs",
                              "shaders/cities_gen/phongFragment.fs",
                              "shaders/cities_gen/textureDisplay.fs",
                              "shaders/cities_gen/test3dFragment.fs",
                              "shaders/cities_gen/landscapeFragment.fs",
                              "shaders/cities_gen/debugSphereFragment.fs",
                              "shaders/cities_gen/waterPlateFragment.fs",
                              "shaders/cities_gen/skyboxFragment.fs",
                              "shaders/cities_gen/roadFragment.fs",
                              "shaders/cities_gen/buildingFragment.fs"};

    std::vector<std::string> names = {"spriteShader",
                                      "simpleShader",
                                      "textureDisplayShader",
                                      "test3dShader",
                                      "landscapeShader",
                                      "debugSphereShader",
                                      "waterPlateShader",
                                      "skyboxShader",
                                      "roadShader",
                                      "buildingShader"};
    for (int i = 0; i < names.size(); i++)
    {
        libra->shaderVector.push_back(new Shader(VSPathes[i], FSPathes[i], names[i]));
    }
    return ent;
}

Shader* ShaderLibrary::GetShaderPointer(std::string title)
{
    Shader* shader = nullptr;
    for (int i = 0; i < shaderVector.size(); i++)
    {
        if (shaderVector[i]->name == title)
        {
            shader = shaderVector[i];
            break;
        }
    }
    return shader;
}

ShaderLibrary* ShaderLibrary::GetLibrary()
{
    if (library == nullptr)
    {
        library = createShaderLibraryTemplate()->GetComponent<ShaderLibrary>();
    }
    return library;    
}   

TransformComponent::TransformComponent(float x, float y, float z, unsigned EID)
    : TransformComponent(EID)
{
    position = glm::vec3(x, y, z);
}

TransformComponent::TransformComponent(unsigned EID) : Component(EID, TransformComponent::typeId) 
{
    rotation = glm::mat4(1.0);
    scale = glm::vec3(1,1,1);
}

TransformComponent::TransformComponent(const TransformComponent& t, unsigned EID)
    :TransformComponent(t.position.x, t.position.y, t.position.z, EID) 
{}

glm::mat4 TransformComponent::GetMatrix()
{
    glm::mat4 res = glm::mat4(1.0f);
    res = glm::translate(res, position);
    res *= rotation;
    res = glm::scale(res, scale);
    return res;
}

glm::mat4 TransformComponent::GetViewMatrixFromThisPoint()
{
    glm::vec4 direction = glm::vec4(TransformComponent::GetDefaultViewDirection(), 0) * rotation;
    return glm::lookAt(position, position + glm::vec3(direction), glm::vec3(0,1,0));
}

Entity* createCameraTemplate(float angleDegrees, float screenRatio, float near, float far, glm::vec3 pos, glm::vec3 dir, std::string tag = "")
{
    auto ent = new Entity();
    if(tag.size() > 0)
    {
        ent->giveTag(tag);
    }
    auto cam = new CameraComponent(ent->id);
    cam->position = pos;
    cam->direction = dir;
    cam->zNear = near;
    cam->zFar = far;
    cam->SetProjectionMatrix(glm::perspective(angleDegrees, screenRatio, near, far));
    cam->viewMatrix = glm::lookAt(pos, dir, glm::vec3(0,1,0));
    // ent->components.push_back(cam);
    ent->AddComponent(cam);
    return ent;    
}

glm::mat4 CameraComponent::GetViewProjectionMatrix()
{
    return projectionMatrix * viewMatrix;
}

glm::mat4 CameraComponent::GetViewMatrix()
{
    return viewMatrix;
}

void CameraComponent::SetProjectionMatrix(glm::mat4 projM)
{
    projectionMatrix = projM;
    invprojectionMatrix = glm::inverse(projectionMatrix);
}

glm::mat4 CameraComponent::GetInvProjMatrix()
{
    return invprojectionMatrix;
}

glm::mat4 CameraComponent::GetProjectionMatrix()
{
    return projectionMatrix;
}

void CameraComponent::ConsoleHandler()
{
    const std::string path = "../Code/Logs/SavedCameraPositions";
    std::function<void(Event*)> f = [&](Event *p)
    {
        if (static_cast<ConsoleInputEvent*>(p)->command != "camera")
            return;

        if (static_cast<ConsoleInputEvent*>(p)->argument == "save")
        {
            std::ofstream file;
            file.open(path);
            if (file.is_open())
            {
                file.clear();
                for (int i = 0; i < 3; i++)
                {
                    file << position[i] << " ";
                }
                for (int i = 0; i < 3; i++)
                {
                    file << direction[i] << " ";
                }
                file.close();
            }
        } 
        else if (static_cast<ConsoleInputEvent*>(p)->argument == "load")
        {
            std::ifstream file;
            file.open(path);
            if (file.is_open())
            {
                float x;
                file >> position.x >> position.y >> position.z;
                file >> direction.x >> direction.y >> direction.z;
                file.close();
            }
        } 
        else if (static_cast<ConsoleInputEvent*>(p)->argument == "fps")
        {
            debug(GameManager::GetInstance()->GetFps());
        } 
        else
        {
            std::string argument = static_cast<ConsoleInputEvent*>(p)->argument;
            std::vector<std::string> argV = Split(argument);
            
            if (argV[0] == "show_sun_pos")
            {
                if (argV.size() != 4)
                {
                    debug("NEW POSITION NEEDED!");
                    throw std::exception{};
                }
                glm::vec3 sunPos;
                for (int i = 0; i < 3; i++)
                {
                    sunPos[i] = std::stof(argV[i + 1]);
                }
                debug("Check it out!");
                std::string debugSphereName = "DebugSun";
                if (Renderer::GetInstance()->DebugSphereDoesExist(debugSphereName))
                    Renderer::GetInstance()->DebugSphereMove(GetPosition() + sunPos, debugSphereName);
                else
                {
                    Renderer::GetInstance()->CreateDebugSphere(
                        1, GetPosition() + sunPos,  glm::vec3{0.7, 0.7, 0}, debugSphereName
                    );
                }
            }
        } 
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("consoleinputevent", f);
}

void CameraComponent::MirrorFromWater(float waterLevel)
{
    position.y = waterLevel - (position.y - waterLevel);
    direction.y = waterLevel - (direction.y - waterLevel);
    viewMatrix = (glm::lookAt(position, direction, glm::vec3 {0, 1, 0}));
}

glm::vec3 CameraComponent::GetPosition()
{
    return position;
}

void CameraComponent::HandleInput()
{
    CameraControlComponent* control = Component::AnyComponentTyped<CameraControlComponent>();
    float mouseDelta[2];

    for (int j = 0; j < 2; j++)
        mouseDelta[j] = (float) InputManager::GetInstance()->mouseDelta[j];

    if(control->state == "free")
    {
        glm::vec3 shift = direction - position;
        float delta = GameManager::GetInstance()->timeDelta;
        shift = glm::rotateY(shift, -control->freeRotateSpeed[0] * delta * mouseDelta[0]);   
        glm::vec3 right = glm::cross(shift, glm::vec3{0,1,0});
        float angle = atan2((shift.y), sqrt(shift.x * shift.x + shift.z * shift.z)); 
        float delta_rotation = control->freeRotateSpeed[1] * delta * mouseDelta[1];

        if(abs(glm::degrees(angle + delta_rotation)) > control->freeVerticalLimit)
        {
            delta_rotation = (angle > 0) ? glm::radians(control->freeVerticalLimit) : glm::radians(-control->freeVerticalLimit);
            delta_rotation -= angle;
        }

        glm::vec3 new_shift = glm::rotate(shift, delta_rotation, right);
        shift = glm::normalize(new_shift); 
        direction = position + shift;
        
        glm::vec3 moveVectors[] = {glm::normalize(direction - position), glm::vec3(0,1,0),
                        -glm::normalize(glm::cross(direction - position, glm::vec3(0,1,0)))};

        float speed = control->freeMovementSpeed * GameManager::GetInstance()->timeDelta;
        auto& butts = InputManager::GetInstance()->pressedButtons;
        if (std::find(butts.begin(),butts.end(), GLFW_KEY_LEFT_SHIFT) != butts.end())
            speed *= 10;

        for (int i = 0; i < InputManager::GetInstance()->pressedButtons.size(); i++)
        {
            switch (InputManager::GetInstance()->pressedButtons[i])
            {
                case GLFW_KEY_W:
                    position += moveVectors[0] * speed;
                    direction += moveVectors[0] * speed;            
                    break;

                case GLFW_KEY_S:
                    position -= moveVectors[0] * speed;
                    direction -= moveVectors[0] * speed;            
                    break;

                case GLFW_KEY_E:
                    position += moveVectors[1] * speed;
                    direction += moveVectors[1] * speed;            
                    break;

                case GLFW_KEY_Q:
                    position -= moveVectors[1] * speed;
                    direction -= moveVectors[1] * speed;            
                    break;

                case GLFW_KEY_A:
                    position += moveVectors[2] * speed;
                    direction += moveVectors[2] * speed;            
                    break;

                case GLFW_KEY_D:
                    position -= moveVectors[2] * speed;
                    direction -= moveVectors[2] * speed;            
                    break;

                default:
                    break;
            }
        }
    }
    else if (control->state == "satellite")
    {
        direction = control->satelliteCenter;
        for (int i = 0; i < InputManager::GetInstance()->pressedButtons.size(); i++)
        {
            if (InputManager::GetInstance()->pressedButtons[i] == GLFW_MOUSE_BUTTON_MIDDLE)
            {
                float delta_rotation_longitude = control->satelliteRotateSpeed[0] * mouseDelta[0];
                control->satelliteLongitudeDegrees += glm::degrees(delta_rotation_longitude);
                float delta_rotation_latitude = -control->satelliteRotateSpeed[1] * mouseDelta[1];
                float newLatitude = glm::degrees(delta_rotation_latitude) + control->satelliteLatitudeDegrees;
                newLatitude = (abs(newLatitude) <= control->satelliteVerticalLimit)
                    ? newLatitude
                    : (newLatitude > 0)
                        ? control->satelliteVerticalLimit
                        : -control->satelliteVerticalLimit;  
                control->satelliteLatitudeDegrees = newLatitude;
            }
        }
        float yPos = control->satelliteRadius * sinf(glm::radians(control->satelliteLatitudeDegrees)) + control->satelliteCenter.y;
        float realRadius = control->satelliteRadius * cosf(glm::radians(control->satelliteLatitudeDegrees));
        float xPos = cosf(glm::radians(control->satelliteLongitudeDegrees)) * realRadius + control->satelliteCenter.x;
        float zPos = sinf(glm::radians(control->satelliteLongitudeDegrees)) * realRadius + control->satelliteCenter.z;
        position = glm::vec3(xPos,yPos,zPos);
    }
    viewMatrix = (glm::lookAt(position, direction, glm::vec3(0,1,0)));
}

glm::vec2 CameraComponent::GetZDists()
{
    return glm::vec2(zNear, zFar);
}

CameraControlComponent::CameraControlComponent(unsigned EID)
    : Component(EID, CameraControlComponent::typeId), freeRotateSpeed{200,200}, satelliteRotateSpeed{4,1}, satelliteZoomBorders{2,10}
{
    state = "satellite";
    switchStateKey = GLFW_KEY_F;

    freeVerticalLimit = 85;
    freeMovementSpeed = 10;

    satelliteVerticalLimit = 80;
    satelliteRadius = 5;
    satelliteLatitudeDegrees = 0;
    satelliteLongitudeDegrees = 90;
    satelliteZoomSpeed = 0.3;
}

void Entity::LoadGltfModel
    (std::string path, std::string obj_name, std::string shaderName)
{
    ModelInfoComponent* model = new ModelInfoComponent(id);
    AddComponent(model);
    BasicGLTFModel* gltf_model = model;
    *gltf_model = (*BasicGLTFModel::LoadGltfModel(path, obj_name, shaderName));
}

ModelInfoComponent::~ModelInfoComponent()
{
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    if (IsUsingIndexes())
    {
        glDeleteBuffers(1, &IBO);
    }
}

SunLightComponent::SunLightComponent(glm::vec3 col, unsigned width, unsigned height, float orthoScale,  unsigned EID)
    :Component(EID, SunLightComponent::GetCompId()), color{col}, SHADOW_WIDTH{width}, SHADOW_HEIGHT{height}
{
    //angleDegrees, screenRatio, near, far);
    projectionMatrix = 
        glm::ortho(-1 * orthoScale, 1 * orthoScale, -1 * orthoScale, 1 * orthoScale, 0.1f, 80.0f);  

    //CREATE FB TEXTURE
	glGenTextures(1, &depthTexture);
	glBindTexture(GL_TEXTURE_2D, depthTexture);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT24, SHADOW_WIDTH, SHADOW_HEIGHT, 0,GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
    glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenTextures(1, &shadowTexture);
    glBindTexture(GL_TEXTURE_2D, shadowTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    //CREATE FBO AND ATTACH TEXTURE TO IT
	glGenFramebuffers(1, &depthMapFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,GL_TEXTURE_2D, depthTexture,0);
	// glDrawBuffer(GL_NONE);
	// glReadBuffer(GL_NONE);

    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status != GL_FRAMEBUFFER_COMPLETE)
		std::cout<<"FBO NO OK!!"<<std::endl;
}

ImageVertexInfoComponent::ImageVertexInfoComponent(unsigned EID)
    :Component(EID, ImageVertexInfoComponent::GetCompId())
{
    shader = ShaderLibrary::GetLibrary()->GetShaderPointer("textureDisplayShader");
    GLfloat vertexes[] = {
     1.0f,  1.0f,   
     1.0f, -1.0f,   
    -1.0f, -1.0f,   
    -1.0f,  1.0f   
    };
    GLuint indices[] = {
        0, 1, 3,
        1, 2, 3  
    };
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);// Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertexes), vertexes, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0); // Note that this is allowed, the call to glVertexAttribPointer registered VBO as the currently bound vertex buffer object so afterwards we can safely unbind
    glBindVertexArray(0); // Unbind VAO (it's always a good thing to unbind any buffer/array to prevent strange bugs), remember: do NOT unbind the EBO, keep it bound to this VAO
}

ImageVertexInfoComponent::~ImageVertexInfoComponent()
{
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

bool SleepInterval::wait(float milisec)
{
    float delta = GameManager::GetInstance()->timeDelta * 1000;
    if (counter <= 0)
    {
        counter = milisec + counter - delta;
        return true;
    }
    counter -= delta;
    return false;
}

bool SleepInterval::interval(float waitMilisec, std::string name)
{
    SleepInterval* sleep = nullptr;
    if (intervals.find(name) == intervals.end())
    {
        sleep =  new SleepInterval();
        intervals.emplace(name, sleep);
    }
    else
    {
        sleep = intervals.find(name)->second;
    }
    return sleep->wait(waitMilisec);
}
