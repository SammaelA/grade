#include <vector>
#include <map>
#include <string>
#include <functional>

class Event
{
    friend class EventHandler;
    public:
        Event();
        int typeId;
};

class EventHandler
{
    public:
        static EventHandler* GetHandler();
        Event* NewEvent(std::string, ...);
        void ApplyFunctionToEvents(std::string, std::function<void(Event*)>);
        void clearEvents();
    private:
        bool IsConvertableTo(Event*, Event*);
        bool IsConvertableTo(Event*, int);
        std::vector<Event*> eventContainer;
        static EventHandler* handler;
        std::map<std::string, int> nameIds;
        int relativeTable[128];
        EventHandler();
        
        
};

class KeyBoardEvent : public Event
{
    friend class EventHandler;
    public:
        int key, action;
    private:
        KeyBoardEvent(int, int);
};

class MouseEvent : public Event
{
    friend  class EventHandler;
    public:
        float pos[2];
        int button, action;
    private:
        MouseEvent(int, int, float, float);
};

class ConsoleInputEvent : public Event
{
    friend class EventHandler;
    public:
        std::string command, argument;
    private:
        ConsoleInputEvent(std::string, std::string);
};