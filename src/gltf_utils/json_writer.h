#pragma once

#include <vector>
#include <string>

class JsonWriter
{
public:
    std::string &get_text() { return text;}
    void start_container();
    void close_container();
    void start_array();
    void close_array();
    
private:
    std::string text;
};