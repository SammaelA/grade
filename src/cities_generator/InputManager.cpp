#include "InputManager.h"
#include "cities_generator/global.h"

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;
extern GLFWwindow* globalWindow;

InputManager* InputManager::m_instance = nullptr;

InputManager::InputManager() :
    mousePos{0, 0}, 
    mouseDelta{0, 0}
{}

InputManager* InputManager::GetInstance()
{
    if (m_instance == nullptr)
        Instantiate();
    return m_instance;
}

void InputManager::Instantiate()
{
    if (m_instance == nullptr)
        m_instance = new InputManager();
}

void InputManager::Update()
{
    glm::vec2 prevDelta {mouseDelta[0], mouseDelta[1]};
    double xpos, ypos;
    glfwGetCursorPos(globalWindow, &xpos, &ypos);
    for (int j = 0; j < 2; j++)
        mouseDelta[j] = -mousePos[j];
                
    mousePos[0] = xpos / SCREEN_WIDTH;
    mousePos[1] = 1 - ypos / SCREEN_HEIGHT;
    
    for (int j = 0; j < 2; j++)
        mouseDelta[j] += mousePos[j];

    if (glm::length(prevDelta) < 0.00001 && 
        glm::length(glm::vec2{mouseDelta[0], mouseDelta[1]}) > MAX_AMPLITUDE)
    {
        mouseDelta[0] = mouseDelta[1] = 0;
    }
}