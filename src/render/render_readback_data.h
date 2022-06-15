#pragma once
struct RenderReadbackInputData
{
    glm::vec2 cursor_screen_pos;
};
struct RenderReadbackData
{
    bool valid = false;
    bool cursor_on_geometry = false;
    glm::vec3 cursor_world_pos; 
};