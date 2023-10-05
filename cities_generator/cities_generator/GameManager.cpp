#include "GameManager.h"
#include "cities_generator/global.h"

GameManager* GameManager::m_instance = nullptr;

GameManager::GameManager()
{
    old_time = std::chrono::high_resolution_clock::now();
}

GameManager* GameManager::GetInstance()
{
    if (m_instance == nullptr)
        Instantiate();
    return m_instance;
}

void GameManager::Instantiate()
{
    if (m_instance == nullptr)
        m_instance = new GameManager();
}

void GameManager::Update()
{
    cur_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> delta = cur_time - old_time;
    timeDelta = (float)delta.count() * 0.001;
    old_time = cur_time;

    static float fpsAccumulator = 0;
    static float frameCounter = 0;
    fpsAccumulator += timeDelta;
    frameCounter++;
    if (fpsAccumulator >= 1.0f)
    {
        fps = frameCounter;
        frameCounter = 0.f;
        fpsAccumulator = 0.f;
    }
}

float GameManager::GetFps()
{
    if (fps < 0)
        return 45;
    else
        return fps;
}