#include "EventClasses.h"
#include <iostream>
#include <stdarg.h>

EventHandler* EventHandler::GetHandler()
{
    if(handler == nullptr)
        handler = new EventHandler();
    return handler;
}

EventHandler::EventHandler()
{
    nameIds = {{"event", 0},
               {"keyboardevent", 1},
               {"mouseevent", 2},
               {"consoleinputevent", 3}};
               
    relativeTable[0] = -1;
    relativeTable[1] = 0;
    relativeTable[2] = 0;
    relativeTable[3] = 0;
}

Event::Event() : typeId(0) {};

KeyBoardEvent::KeyBoardEvent(int k, int a) : key(k), action(a)
{
    typeId = 1;
}

MouseEvent::MouseEvent(int b, int a, float x, float y) : button(b), action(a), pos{x,y}
{
    typeId = 2;
}

ConsoleInputEvent::ConsoleInputEvent(std::string name, std::string arg) : command{name}, argument{arg}
{
    typeId = 3;
}

Event* EventHandler::NewEvent(std::string name, ...)
{
    Event* result;
    auto it = nameIds.find(name);
    int type;
    if (it == nameIds.end())
    {
        std::cout << "BAD EVENT NAME!\n";
        return nullptr;
    }
    type = it->second;
    va_list vl;
    va_start(vl, name);
    switch(type)
    {
        case 0:
            result = new Event();
            break;
        case 1:
            result = new KeyBoardEvent(va_arg(vl, int), va_arg(vl, int));
            break;
        case 2:
            result = new MouseEvent(va_arg(vl, int), va_arg(vl, int), (float)va_arg(vl, double), (float)va_arg(vl, double));
            break;
        case 3:
            result = new ConsoleInputEvent(va_arg(vl, const char*), va_arg(vl, const char*));
            break;
        default:
            break;
    }
    va_end(vl);
    eventContainer.push_back(static_cast<Event*>(result));
    return result;
}

bool EventHandler::IsConvertableTo(Event* from, Event* to)
{
    int curId = from->typeId;
    while(curId != -1)
    {
        if(curId == to->typeId)
            return true;
        curId = relativeTable[curId];
    }
    return false;
}

bool EventHandler::IsConvertableTo(Event* from, int to)
{
    int curId = from->typeId;
    while(curId != -1)
    {
        if(curId == to)
            return true;
        curId = relativeTable[curId];
    }
    return false;
}

void EventHandler::clearEvents()
{
    for (size_t i = 0; i < eventContainer.size(); i++)
    {
        delete eventContainer[i];
    }
    eventContainer.resize(0);
    
}

void EventHandler::ApplyFunctionToEvents(std::string name, std::function<void(Event* e)> func)
{
    auto it = nameIds.find(name);
    int type;
    if (it == nameIds.end())
    {
        std::cout << "BAD EVENT NAME!\n";
        return;
    }
    type = it->second;
    for (size_t i = 0; i < eventContainer.size(); i++)
    {
        if (IsConvertableTo(eventContainer[i], type))
            func(eventContainer[i]);
    }
    
}