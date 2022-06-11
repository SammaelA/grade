#include <chrono>

class GameManager
{
    GameManager();
    static GameManager* m_instance;
    std::chrono::_V2::system_clock::time_point old_time, cur_time;
    float fps = 1;
    
    public:
        static void Instantiate();
        static GameManager* GetInstance();
        void Update();
        float GetFps();
        
        float timeDelta;
};