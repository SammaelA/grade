
#include <vector>
class InputManager
{
    InputManager();
    static InputManager* m_instance;
    
    public:
        static void Instantiate();
        static InputManager* GetInstance();
        void Update();

        double mousePos[2];
        double mouseDelta[2];
        const float MAX_AMPLITUDE = 0.2;
        std::vector<int> pressedButtons;
};