#pragma once
struct RenderReadbackInputData
{
    float2 cursor_screen_pos = float2(0,0);
};
struct RenderReadbackData
{
    bool valid = false;
    bool cursor_on_geometry = false;
    float4 cursor_world_pos_type; 
};