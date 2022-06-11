#include <iostream>
#include <string>
#include <chrono>

/*
unsigned char* SOIL_load_image_safe(const char* path, int* w, int* h, 
                                    int* channels, int force_channels,
                                    std::string error_message);*/

class MyClock
{
    static std::chrono::milliseconds previousTime;
    public:
        static void Tick();
        static int64_t Tock();
};
