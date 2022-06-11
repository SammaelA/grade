#include "externFunctions.h"

/*
unsigned char* SOIL_load_image_safe(const char* path, int* w, int* h, 
                                    int* channels, int force_channels,
                                    std::string error_message)
{
    auto res = SOIL_load_image(path, w, h, channels, force_channels);
    if (res == 0)
    {
        std::cout << "______SOIL_load_image failed______\n" << error_message << std::endl;
        throw std::exception{};
    }
    return res;
}
*/
std::chrono::milliseconds MyClock::previousTime;

void MyClock::Tick()
{
    previousTime = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch());
}

int64_t MyClock::Tock()
{
    std::chrono::milliseconds curTime = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch());
    return (curTime - previousTime).count();
}
