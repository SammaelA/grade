#include <iostream>
#define GLEW_STATIC
#include <cmath>
#include <chrono>
#include <csignal>
#include <algorithm>
#include <unordered_map>
#include <sstream>

#include "./EventClasses.h"
#include "./shaderClass.h"
#include "./ECSClass.h"
#include "./LandMeshGenerator/LandGeneratorClass.h"
#include "./ExternFunctions/externFunctions.h"
#include "cities_generator/global.h"
#include "./static_init.h"
#include "third_party/stb_image.h"
#include "cities_generator/Renderer.h"
#include "InputManager.h"
#include "GameManager.h"
#include "Skybox.h"

using namespace glm;

//Init and destroy functions
GLFWwindow* MyInit(GLuint&,GLuint&);
void initializeGame();
void setCustomCursor();
void DeleteSingletones();

//Input handlers
void WindowCloseEventHandler(GLFWwindow*);
void WindowButtonHandler(GLFWwindow*);
void CorrectHeldButtonsHandler();
void ChangeCameraTypeHandler();
void ControlSatelliteCameraHandler();
void DebugConsoleInputHandler();
void SavedHeightmapSwapHandler(Landscape*);

//Systems
void SpriteWithTagRenderSystem(std::string);
void SunShadowsRenderSystem();
void GeneralRenderSystem();
void clearEventsSystem();

//UNKNOWN
void BindTextureToFrame(GLint frameTex = -1, GLint depthT = -1);

//Callback functions for handling user input
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode); //My handmade key enets handler
void mouse_callback(GLFWwindow*, int, int, int); 
void scroll_callback(GLFWwindow*, double ,double);

//Global window variables
GLuint SCREEN_WIDTH, SCREEN_HEIGHT;
GLFWwindow* globalWindow;


int cities_generator_main()
{
    globalWindow = MyInit(SCREEN_WIDTH, SCREEN_HEIGHT);
    glfwSetKeyCallback(globalWindow, key_callback); // Set the required callback functions
    glfwSetMouseButtonCallback(globalWindow, mouse_callback);
    glfwSetScrollCallback(globalWindow, scroll_callback);

    //Generation object
    Landscape landscape;
    landscape.Init();

    //skybox object
    Skybox skybox(std::string(ASSETS_FOLDER)+"Textures/Skybox/", glm::normalize(glm::vec3{22, 17, 4}), M_PI_2 * -1.0/3);

    //Sun Entity
    Entity* sun = new Entity();
    sun->AddComponent(new SunLightComponent({1,1,1}, 8000, 8000, 16.0f, sun->id));
    TransformComponent* sunTransform = new TransformComponent(sun->id);
    sunTransform->position = skybox.GetSunPosition();
    sunTransform->rotation = glm::eulerAngleZY(glm::pi<float>() / 4 * 0.85f, glm::pi<float>() / 4 * 5.0f);
    sun->AddComponent(sunTransform);

    //Initialize game
    initializeGame();
    
    SunShadowsRenderSystem(); // For static objects

    // Game loop
    while (!glfwWindowShouldClose(globalWindow))
    {   
            
        // EVENT HANDLING
        clearEventsSystem();
        glfwPollEvents();
        WindowCloseEventHandler(globalWindow);
        WindowButtonHandler(globalWindow);
        CorrectHeldButtonsHandler();
        ChangeCameraTypeHandler();
        ControlSatelliteCameraHandler();
        DebugConsoleInputHandler();
        Component::AnyComponentTyped<CameraComponent>()->ConsoleHandler();
        SavedHeightmapSwapHandler(&landscape);
        Component::AnyComponentTyped<CameraComponent>()->HandleInput();
        
        // SYSTEMS
        Renderer::GetInstance()->MyBufferClearSystem(1, 1, 1);
        InputManager::GetInstance()->Update();
        GameManager::GetInstance()->Update();
        landscape.Update();
        
        
        //RENDER

        // SpriteWithTagRenderSystem("BackgroundImage"); 
        GeneralRenderSystem();

        if (landscape.IsRenderingWater())
        {
            WaterPlate* water = landscape.GetWaterPlate();

            glBindFramebuffer(GL_FRAMEBUFFER, water->FBO);
            BindTextureToFrame(water->bloomTexture, water->depthTexture);
            Renderer::GetInstance()->MyBufferClearSystem(0,0,0);
            water->Render(RENDER_MODE::BLOOM);

            CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
            cam->MirrorFromWater(landscape.GetWaterLevel());
            glBindFramebuffer(GL_FRAMEBUFFER, water->FBO);
            BindTextureToFrame(water->aboveWaterTexture, water->depthTexture);
            Renderer::GetInstance()->MyBufferClearSystem(1,1,1);
            landscape.SetClippingSettings(
                glm::vec4 {
                    0, 1, 0, -landscape.GetWaterLevel() + landscape.GetWaterLevelSafeDelta()
                }
            );
            landscape.Render(RENDER_MODE::CLIPPING);
            landscape.RenderBuildings(RENDER_MODE::DEFUALT);
            skybox.Render();
            cam->MirrorFromWater(landscape.GetWaterLevel());

            glBindFramebuffer(GL_FRAMEBUFFER, water->FBO);
            BindTextureToFrame(water->underWaterTexture, water->depthTexture);
            Renderer::GetInstance()->MyBufferClearSystem(1,1,1);
            landscape.SetClippingSettings(
                glm::vec4 {
                    0, -1, 0, landscape.GetWaterLevel() + landscape.GetWaterLevelSafeDelta()
                }
            );
            landscape.Render(RENDER_MODE::CLIPPING);
            
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            landscape.Render(RENDER_MODE::DEFUALT);
            water->Render(RENDER_MODE::DEFUALT);
        }
        else
        {
            landscape.Render(RENDER_MODE::DEFUALT);
            
        }
        skybox.Render();
        landscape.RenderRoad(RENDER_MODE::DEFUALT);
        landscape.RenderBuildings(RENDER_MODE::DEFUALT);
        Renderer::GetInstance()->RenderDebugSpheres();
        
        static MAP currentMap = MAP::NONE;
        std::function<void(Event*)> f = [&](Event *p)
        {
            if (static_cast<KeyBoardEvent*>(p)->action == GLFW_PRESS)
            {
                if (static_cast<KeyBoardEvent*>(p)->key == GLFW_KEY_M)
                {
                    ++currentMap;
                    if (currentMap == MAP::ENUM_SIZE)
                    {
                        currentMap = MAP::NONE;
                    }
                }
            } 
        };
        EventHandler::GetHandler()->ApplyFunctionToEvents("keyboardevent", f);
        landscape.RenderMap(currentMap);
        

        //All ECS Systems
        System::RunAllSystems();

        //Suppotrs dynamic archetypes
        //MUST BE AFTER ALL ARCHETYPE MANIPULATION
        Archetype::handleNewArchetypes();

        // SWAP THE SCREEN BUFFERS
        glfwSwapBuffers(globalWindow);
    }

    DeleteSingletones();

    glfwTerminate();
    
    return 0;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    EventHandler::GetHandler()->NewEvent("keyboardevent", key, action);
} 

void mouse_callback(GLFWwindow* window, int button, int action, int mods)
{
    double xpos, ypos;
    glfwGetCursorPos(globalWindow, &xpos, &ypos);
    EventHandler::GetHandler()->NewEvent("mouseevent", button, action, xpos, SCREEN_HEIGHT - ypos);
} 

void scroll_callback(GLFWwindow* window, double xoff, double yoff)
{
    EventHandler::GetHandler()->NewEvent("mouseevent", 3, GLFW_PRESS, xoff, yoff);
} 

void WindowCloseEventHandler(GLFWwindow *window)
{
    std::function<void(Event*)> f = [&](Event *p)
    {
        if (static_cast<KeyBoardEvent*>(p)->key == GLFW_KEY_ESCAPE)
            glfwSetWindowShouldClose(window, GL_TRUE);
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("keyboardevent", f);
};

void WindowButtonHandler(GLFWwindow *window)
{
    std::function<void(Event*)> f = [&](Event *p)
    {
        MouseEvent* event = static_cast<MouseEvent*>(p);
        if (event->action == GLFW_PRESS)
        {
            
            std::vector<ButtonComponent*> buts;
            std::vector<Entity*> entitiesWithButs = 
                std::move(Entity::AllEntitiesWithComponentsAndTags(1, 0, 
                    std::vector<unsigned>{ButtonComponent::GetCompId()},
                    std::vector<std::string>{}, true));

            for (size_t i = 0; i < entitiesWithButs.size(); i++)
            {
                ButtonComponent* but = entitiesWithButs[i]->GetComponent<ButtonComponent>(true);
                if (but != nullptr)
                {
                    if (event->pos[0] < but->pos[0] || event->pos[0] > but->pos[0] + but->size[0] ||
                        event->pos[1] < but->pos[1] || event->pos[1] > but->pos[1] + but->size[1])
                        continue;
                    buts.push_back(but);
                }
            }
            if (buts.size() > 0)
            {
                std::sort(buts.begin(), buts.end(),
                    [] (ButtonComponent* const& a, ButtonComponent* const& b) { return a->layer > b->layer; });
                
                buts[0]->func(event->pos[0], event->pos[1], event->button);
            }
        }
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("mouseevent", f);
}

void clearEventsSystem()
{
    EventHandler::GetHandler()->clearEvents();
}

void SpriteWithTagRenderSystem(std::string tag)
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    GLint TransformUniformPos = 0;
    bool firstSprite = true;
    std::vector<Entity*> spriteEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
                                            1, 0, std::vector<unsigned>{SpriteComponent::GetCompId()},
                                            std::vector<std::string>{}, true));
    for (size_t i = 0; i < spriteEntities.size(); i++)
    {
        if (!spriteEntities[i]->hasTag(tag))
            continue;

        SpriteComponent *sprite = spriteEntities[i]->GetComponent<SpriteComponent>();

        if (firstSprite)
        {
            TransformUniformPos = glGetUniformLocation(sprite->spriteShader->Program, "Transform");
            firstSprite = false;
        }

        sprite->spriteShader->Use();
        glUniform4f(TransformUniformPos, sprite->size[0] / SCREEN_WIDTH, sprite->size[1] / SCREEN_HEIGHT, sprite->pos[0] / SCREEN_WIDTH, sprite->pos[1] / SCREEN_HEIGHT);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sprite->EBO);
        glBindVertexArray(sprite->VAO);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, sprite->spriteTexture);
        glUniform1i(glGetUniformLocation(sprite->spriteShader->Program, "Texture"), 0);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    }
    glDisable(GL_BLEND);
}

void SavedHeightmapSwapHandler(Landscape* land)
{
    std::function<void(Event*)> f = [&](Event *p)
    {
        if (static_cast<KeyBoardEvent*>(p)->action != GLFW_PRESS)
            return;
        if (static_cast<KeyBoardEvent*>(p)->key == GLFW_KEY_UP)
            land->SwapHeightMaps();
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("keyboardevent", f);
};

GLFWwindow* MyInit(GLuint &w, GLuint &h)
{
    std::cout << "Starting GLFW context, OpenGL 3.3" << std::endl;
    // Init GLFW
    glfwInit();
    // Set all the required options for GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    // Create a GLFWwindow object that we can use for GLFW's functions
    const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    w = mode->width;
    h = mode->height;
    GLFWwindow* window = glfwCreateWindow(w, h, "NumericDictators_prerender", glfwGetPrimaryMonitor(), nullptr);
    glfwMakeContextCurrent(window);
    // Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
    glewExperimental = GL_TRUE;
    // Initialize GLEW to setup the OpenGL Function pointers
    glewInit();
    checkForGlErrors("JUST CLEAN glewInit errors", true);
    // Define the viewport dimensions
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);  
    glViewport(0, 0, width, height);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GLFW_TRUE);
    return window;
}

void initializeGame()
{
    setCustomCursor();

    Entity* BackgroungImage = createImageTemplate(0,0,SCREEN_WIDTH,SCREEN_HEIGHT, "./Assets/Common/background.jpg");
    BackgroungImage->giveTag("BackgroundImage");
    // BackgroungImage->components.push_back(new ButtonComponent)

    // Entity* TestModel = new Entity();
    // TestModel->LoadGltfModel
    //     ("./Assets/Models/Full_Desk/Full_Desk.gltf", "BLP_Env_Globe", "simpleShader");

    // Entity* TestModel1 = new Entity();
    // TestModel1->LoadGltfModel
    //     ("./Assets/Models/Full_Desk/Full_Desk.gltf", "BLP_Env_DeskDrawer", "simpleShader");

    // Entity* TestModel2 = new Entity();
    // TestModel2->LoadGltfModel
    //     ("./Assets/Models/Full_Desk/Full_Desk.gltf", "BLP_Env_Book", "simpleShader");
    
    // Entity* TestModel3 = new Entity();
    // TestModel3->LoadGltfModel
    //     ("./Assets/Models/Full_Desk/Full_Desk.gltf", "BLP_Env_Desk", "simpleShader");

    // Entity* TestModel4 = new Entity();
    // TestModel4->LoadGltfModel
    //     ("../Assets/Models/Sign/sign.gltf", "Plane", "simpleShader");
    
    // Entity* TestModel5 = new Entity();
    // TestModel5->LoadGltfModel
    //     ("../Assets/Models/Floor/floor.gltf", "Cube", "simpleShader");

    // Entity* TestModel7 = new Entity();
    // TestModel7->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.001", "simpleShader");

    // Entity* TestModel8 = new Entity();
    // TestModel8->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.002", "simpleShader");

    // Entity* TestModel9 = new Entity();
    // TestModel9->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.003", "simpleShader");

    // Entity* TestModel10 = new Entity();
    // TestModel10->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.004", "simpleShader");

    // Entity* TestModel11 = new Entity();
    // TestModel11->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.005", "simpleShader");

    // Entity* TestModel12 = new Entity();
    // TestModel12->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.006", "simpleShader");

    // Entity* TestModel13 = new Entity();
    // TestModel13->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cube.007", "simpleShader");

    // Entity* TestModel14 = new Entity();
    // TestModel14->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cone", "simpleShader");

    // Entity* TestModel15 = new Entity();
    // TestModel15->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cone.001", "simpleShader");

    // Entity* TestModel16 = new Entity();
    // TestModel16->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Cone.002", "simpleShader");

    // Entity* TestModel17 = new Entity();
    // TestModel17->LoadGltfModel
    //     ("../Assets/Models/Road/road.gltf", "Cube", "simpleShader");

    // Entity* TestModel18 = new Entity();
    // TestModel18->LoadGltfModel
    //     ("../Assets/Models/LowPolyScene/Scene.gltf", "Sphere", "simpleShader");

    // Entity* one = new Entity();
    // one->AddComponent(new CompA(one->id));

    // Entity* two = new Entity();
    // two->AddComponent(new CompB(two->id));

    // Entity* three = new Entity();
    // three->AddComponent(new CompA(three->id));

    // debug("FFFF");
    // debug(Entity::entitiesArchetype[two->id].archetype->components[0].size());
    // debug(Entity::entitiesArchetype[one->id].archetype->components[0].size()); 
    // debug(Entity::entitiesArchetype[three->id].archetype->components[0].size()); 
    // debug(Archetype::basicNode->entitiesId.size());
    // Archetype::debugShowTree();
    // debug("LLLLLLL");
    // unsigned u;
    // debug(three->HasComponent(CompA::GetCompId()));
    // debug(three->HasComponent(CompB::GetCompId()));


    // Entity* TestCube = new Entity();
    // TestCube->components.push_back(new ModelInfoComponent("simpleShader", TestCube->id));
    // TestCube->components.push_back(new TransformComponent(1,0,-1, TestCube->id));
    // Entity* TestCube1 = new Entity();
    // TestCube1->components.push_back(new ModelInfoComponent("simpleShader", TestCube1->id));
    // TestCube1->components.push_back(new TransformComponent(-1,0,-1, TestCube1->id));
    // Entity* TestCube2 = new Entity();
    // TestCube2->components.push_back(new ModelInfoComponent("simpleShader", TestCube2->id));
    // TestCube2->components.push_back(new TransformComponent(1,0,1, TestCube2->id));
    // Entity* TestCube3 = new Entity();
    // TestCube3->components.push_back(new ModelInfoComponent("simpleShader", TestCube3->id));
    // TestCube3->components.push_back(new TransformComponent(-1,0,1, TestCube3->id));
    auto cam = createCameraTemplate(45.0f, (float)SCREEN_WIDTH/(float)SCREEN_HEIGHT, 0.1f, 1000.0f, glm::vec3(0,0,0), glm::vec3(1,0,0), "mainCamera");
    // cam->components.push_back(new CameraControlComponent(cam->id));
    cam->AddComponent(new CameraControlComponent(cam->id));

    // Renderer::GetInstance()->CreateDebugSphere(3, glm::vec3{0,4,0}, glm::vec3{0.8, 0.8, 0});
    // Renderer::GetInstance()->CreateDebugSphere(1, glm::vec3{0,1,0}, glm::vec3{0.1, 0.8, 0});
}

void GeneralRenderSystem()
{
    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    Renderer::GetInstance()->MyBufferClearSystem(0,0,0);

    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // GLint TransformUniformPos;
    // bool firstSprite = true;
    std::vector<Entity*> modelEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
                                            1, 0, std::vector<unsigned>{ModelInfoComponent::GetCompId()},
                                            std::vector<std::string>{}, true));
                                            
    glm::mat4 viewProj = glm::mat4(1.0f);
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    
    viewProj = cam->GetViewProjectionMatrix();
             
    for (size_t i = 0; i < modelEntities.size(); i++)
    {
        ModelInfoComponent *model = modelEntities[i]->GetComponent<ModelInfoComponent>();
        
        Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
                2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];

        TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
        SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();

        model->UseShader();

        GLint ViewProjMatrixUniformPos, GlobalIlluminationPointUniformPos,
                GlobalIlluminationColorUniformPos, diffuseColorPos, isIsotropicColorPos,
                isDepthOnlyPos, GlobalIlluminationViewProjUniformPos, CameraPositionPos;
        
        ViewProjMatrixUniformPos = glGetUniformLocation(model->m_shader->Program, "view_proj_matrix");
        GlobalIlluminationPointUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_pos");
        GlobalIlluminationColorUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_color"); 
        GlobalIlluminationViewProjUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_view_proj_matrix");
        isDepthOnlyPos = glGetUniformLocation(model->m_shader->Program, "is_depth_only"); 
        CameraPositionPos = glGetUniformLocation(model->m_shader->Program, "camera_position"); 

        // std::vector<Entity*> sunLightEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
        //     1, 0, std::vector<unsigned>{SunLightComponent::GetCompId()}, std::vector<std::string>{}, true));
        // SunLightComponent* sunS = sunLightEntities[0]->GetComponent<SunLightComponent>();
        // TransformComponent* sunT = sunLightEntities[0]->GetComponent<TransformComponent>();
        glm::mat4 sunProjM = sunComp->projectionMatrix;
        glm::mat4 sunVeiwM = sunTransform->GetViewMatrixFromThisPoint();
        glm::mat4 sunViewProj = sunProjM * sunVeiwM;

        model->LoadInShader_Scene3DObject();
        glUniformMatrix4fv(ViewProjMatrixUniformPos, 1, GL_FALSE, glm::value_ptr(viewProj));
        glUniformMatrix4fv(GlobalIlluminationViewProjUniformPos, 1, GL_FALSE, glm::value_ptr(sunViewProj));
        glUniform3f(GlobalIlluminationPointUniformPos, sunTransform->position.x, sunTransform->position.y, sunTransform->position.z);
        glUniform3f(GlobalIlluminationColorUniformPos, sunComp->color.x, sunComp->color.y, sunComp->color.z);
        glUniform3fv(CameraPositionPos, 1, glm::value_ptr(cam->GetPosition()));
        glUniform1i(isDepthOnlyPos, false);
        
        //TEXTURES
        int TEX_COUNT = 0;
        model->LoadInShader_BasicTexturedObject(TEX_COUNT);
        glActiveTexture(GL_TEXTURE0 + TEX_COUNT);
        glBindTexture(GL_TEXTURE_2D, sunComp->shadowTexture);
        glUniform1i(glGetUniformLocation(model->m_shader->Program, "shadow_texture"), TEX_COUNT);
        
        //DRAW CALL
        model->Render(RENDER_MODE::DEFUALT);
    }
}

void SunShadowsRenderSystem()
{
    // glEnable(GL_CULL_FACE);
    // glCullFace(GL_FRONT); 
    glEnable(GL_DEPTH_TEST);
    
    std::vector<Entity*> sunLightEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
        1, 0, std::vector<unsigned>{SunLightComponent::GetCompId()}, std::vector<std::string>{}, true));
    std::vector<Entity*> modelEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
        1, 0, std::vector<unsigned>{ModelInfoComponent::GetCompId()}, std::vector<std::string>{}, true));
    

    for (int j = 0; j < sunLightEntities.size(); j++)
    {
        SunLightComponent* sunComp = sunLightEntities[j]->GetComponent<SunLightComponent>();
        TransformComponent* sunT = sunLightEntities[j]->GetComponent<TransformComponent>();

        glm::mat4 projM = sunComp->projectionMatrix;
        glm::mat4 veiwM = sunT->GetViewMatrixFromThisPoint();
        glm::mat4 viewProj = projM * veiwM;

        glViewport(0, 0, sunComp->SHADOW_WIDTH, sunComp->SHADOW_HEIGHT);
        glBindFramebuffer(GL_FRAMEBUFFER, sunComp->depthMapFBO);
        BindTextureToFrame(sunComp->shadowTexture, sunComp->depthTexture);
        Renderer::GetInstance()->MyBufferClearSystem(1,0,0);
        
        for (size_t i = 0; i < modelEntities.size(); i++)
        {
            ModelInfoComponent *model = modelEntities[i]->GetComponent<ModelInfoComponent>();
            
            glm::mat4 transform;
            if (modelEntities[i]->HasComponent(TransformComponent::GetCompId()))
            {
                transform = modelEntities[i]->GetComponent<TransformComponent>()->GetMatrix();
            }
            else
            {
                transform = glm::mat4(1.0f);
            }
            
            model->UseShader();

            GLint TransformUniformPos, ViewProjMatrixUniformPos, isDepthOnlyPos;
            
            TransformUniformPos = glGetUniformLocation(model->m_shader->Program, "transform_matrix");
            ViewProjMatrixUniformPos = glGetUniformLocation(model->m_shader->Program, "view_proj_matrix");
            isDepthOnlyPos = glGetUniformLocation(model->m_shader->Program, "is_depth_only"); 

            glUniformMatrix4fv(TransformUniformPos, 1, GL_FALSE, glm::value_ptr(transform));
            glUniformMatrix4fv(ViewProjMatrixUniformPos, 1, GL_FALSE, glm::value_ptr(viewProj));
            glUniform1i(isDepthOnlyPos, true);
            
            //DRAW CALL
            glBindVertexArray(model->VAO);
            if (model->IsUsingIndexes())
            {
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model->IBO);
                glDrawElements(GL_TRIANGLES, model->indicesNum, model->Get_IBO_element_type(), 0);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            }
            else
            {
                glDrawArrays(GL_TRIANGLES, 0, model->vertexesNum);
            }
            glBindVertexArray(0);
        }
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // glCullFace(GL_BACK); 
    // glDisable(GL_CULL_FACE);
    
	// GLsync fence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, NULL); //Create the fence
    // GLenum state = glClientWaitSync(fence, GL_SYNC_FLUSH_COMMANDS_BIT, 1000000000); //Timeout after 1 second
	// glDeleteSync(fence); //Delete the fence object
	// if (state != GL_CONDITION_SATISFIED) //If anything other than satisfied. Print it out
	// 	std::cout << "Fence blocked! CODE: " << state << std::endl;
}


void CorrectHeldButtonsHandler()
{
    std::function<void(Event*)> f_key = [&](Event *p)
    {
        if (static_cast<KeyBoardEvent*>(p)->action == GLFW_PRESS)
        {
            int k = static_cast<KeyBoardEvent*>(p)->key;
            for (int i = 0; i < InputManager::GetInstance()->pressedButtons.size(); i++)
            {
                if(k == InputManager::GetInstance()->pressedButtons[i])
                    return;
            }
            InputManager::GetInstance()->pressedButtons.push_back(k);
        } 
        else if (static_cast<KeyBoardEvent*>(p)->action == GLFW_RELEASE)
        {
            int k = static_cast<KeyBoardEvent*>(p)->key;
            for (int i = 0; i < InputManager::GetInstance()->pressedButtons.size(); i++)
            {
                if(k == InputManager::GetInstance()->pressedButtons[i])
                {
                    InputManager::GetInstance()->pressedButtons[i] = InputManager::GetInstance()->pressedButtons[InputManager::GetInstance()->pressedButtons.size() - 1];
                    InputManager::GetInstance()->pressedButtons.resize(InputManager::GetInstance()->pressedButtons.size() - 1);
                    return;
                }
            }   
        }  
    };
    std::function<void(Event*)> f_mouse = [&](Event *p)
    {
        if (static_cast<MouseEvent*>(p)->action == GLFW_PRESS)
        {
            int k = static_cast<MouseEvent*>(p)->button;
            for (int i = 0; i < InputManager::GetInstance()->pressedButtons.size(); i++)
            {
                if(k == InputManager::GetInstance()->pressedButtons[i])
                    return;
            }
            InputManager::GetInstance()->pressedButtons.push_back(k);
        } 
        else if (static_cast<MouseEvent*>(p)->action == GLFW_RELEASE)
        {
            int k = static_cast<MouseEvent*>(p)->button;
            for (int i = 0; i < InputManager::GetInstance()->pressedButtons.size(); i++)
            {
                if(k == InputManager::GetInstance()->pressedButtons[i])
                {
                    InputManager::GetInstance()->pressedButtons[i] = InputManager::GetInstance()->pressedButtons[InputManager::GetInstance()->pressedButtons.size() - 1];
                    InputManager::GetInstance()->pressedButtons.resize(InputManager::GetInstance()->pressedButtons.size() - 1);
                    return;
                }
            }   
        }  
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("keyboardevent", f_key);
    EventHandler::GetHandler()->ApplyFunctionToEvents("mouseevent", f_mouse);
}

void ChangeCameraTypeHandler()
{
    Entity* mainCameraObjet = Entity::AnyObjectWithTag("mainCamera");
    if (!mainCameraObjet->HasComponent(CameraControlComponent::GetCompId()) ||
        !mainCameraObjet->HasComponent(CameraComponent::GetCompId()))
    {
        debug("MAINCAMERA ERROR");
        throw std::exception{};
    }
    CameraControlComponent* control = mainCameraObjet->GetComponent<CameraControlComponent>();

    std::function<void(Event*)> f = [&](Event *p)
    {
        if (static_cast<KeyBoardEvent*>(p)->action == GLFW_PRESS)
        {
            if (static_cast<KeyBoardEvent*>(p)->key == control->switchStateKey)
            {
                if (control->state == "free")
                {
                    control->state = "satellite";
                    glfwSetInputMode(globalWindow, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
                }
                else if (control->state == "satellite")
                {
                    control->state = "free";
                    glfwSetInputMode(globalWindow, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
                }
            }
        } 
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("keyboardevent", f);
}

void ControlSatelliteCameraHandler()
{
    CameraControlComponent* control = Component::AnyComponentTyped<CameraControlComponent>();
    
    std::function<void(Event*)> f_mouse = [&](Event *p)
    {
        if (static_cast<MouseEvent*>(p)->button == 3 && static_cast<MouseEvent*>(p)->action == GLFW_PRESS)
        {
            control->satelliteRadius -= static_cast<MouseEvent*>(p)->pos[1] * control->satelliteZoomSpeed;
            if (control->satelliteRadius > control->satelliteZoomBorders[1])
                control->satelliteRadius = control->satelliteZoomBorders[1];
            if (control->satelliteRadius < control->satelliteZoomBorders[0])
                control->satelliteRadius = control->satelliteZoomBorders[0];
            
            // std::cout << static_cast<MouseEvent*>(p)->button << " " << static_cast<MouseEvent*>(p)->pos[1] << std::endl;
        }   
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("mouseevent", f_mouse);
}

void DebugConsoleInputHandler()
{
    std::function<void(Event*)> f = [&](Event *p)
    {
        if (static_cast<KeyBoardEvent*>(p)->action == GLFW_PRESS)
        {
            if (static_cast<KeyBoardEvent*>(p)->key == GLFW_KEY_ENTER)
            {
                glfwIconifyWindow(globalWindow);

                debug("ENTER COMMAND!");

                char buffer[256];
                std::cin.getline(buffer, 256);
                std::string input = buffer;
                std::vector<std::string> inputV = Split(input);

                std::string commandName, argument;
                commandName = inputV[0];
                for (int i = 1; i < inputV.size(); i++)
                {
                    argument += inputV[i];
                    if (i != inputV.size() -1)
                        argument += " ";
                }
                
                EventHandler::GetHandler()->NewEvent(
                    "consoleinputevent", 
                    commandName.c_str(),
                    argument.c_str()
                );
            }
        } 
    };
    EventHandler::GetHandler()->ApplyFunctionToEvents("keyboardevent", f);
}

void intersect(std::unordered_set<unsigned>& reduced, const std::unordered_set<unsigned>& substracted)
{
    for (auto it = reduced.begin(); it != reduced.end(); )
    {
        if (substracted.count(*it) == 0)
        {
            it = reduced.erase(it);
        }
        else
        {
            it++;
        }   
    }
}

void setCustomCursor()
{
    GLFWimage image;
    std::string path = std::string(ASSETS_FOLDER)+"Common/cursor.png";
    image.pixels = stbi_load(path.c_str(), &image.width, &image.height, 0, 4);
    GLFWcursor * cursor = glfwCreateCursor ( &image, 0, 0);
    glfwSetCursor(globalWindow, cursor);
}

void BindTextureToFrame(GLint frameTex, GLint depthT)
{
    if (frameTex != -1)
    {
        glBindTexture(GL_TEXTURE_2D, frameTex); // Need to generate texture object
        glBindTexture(GL_TEXTURE_2D, 0);

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, (GLuint)frameTex, 0); 
    }
    if (depthT != -1)
    {
        glBindTexture(GL_TEXTURE_2D, depthT); // Need to generate texture object
        glBindTexture(GL_TEXTURE_2D, 0);

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,GL_TEXTURE_2D, (GLuint)depthT,0);
        // glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, (GLuint)depthT, 0); 
        // glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, (GLuint)depthT, 0); 
    }
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
	    std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_UNDEFINED)
            debug("1");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
            debug("2");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT)
            debug("3");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER)
            debug("4");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER)
            debug("5");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_UNSUPPORTED)
            debug("6");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE)
            debug("7");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE)
            debug("8");
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS)
            debug("9");
        throw std::exception{};
    }
}

void DeleteSingletones()
{
    Renderer::GetInstance()->Destroy();
}