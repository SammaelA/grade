#include "common_utils/timestamp.h"
#include "common_utils/utility.h"

Timestamp::Timestamp(bool _need_debug)
{
    need_debug = _need_debug;
}
void Timestamp::start(std::string name)
{
    auto it = TS.find(name);
    if (it == TS.end())
    {
        Stamp s;
        s.startQuery = create_query();
        s.endQuery = create_query();
        s.started = false;
        s.ended = false;

        TS.emplace(name, s);
        it = TS.find(name);
    }
    if (it != TS.end())
    {
      glQueryCounter(it->second.startQuery, GL_TIMESTAMP);
      it->second.started = true;
    }
}
void Timestamp::end(std::string name)
{
    auto it = TS.find(name);
    if (it == TS.end())
    {
        Stamp s;
        s.startQuery = create_query();
        s.endQuery = create_query();
        s.started = false;
        s.ended = false;

        TS.emplace(name, s);
        it = TS.find(name);
    }
    if (it != TS.end())
    {
      glQueryCounter(it->second.endQuery, GL_TIMESTAMP);
      it->second.ended = true;
    }
}
void Timestamp::resolve()
{
    ticks++;
    for (auto &ts : TS)
    {
        if (ts.second.started)
        {
            glGetQueryObjectui64v ( ts.second.startQuery, GL_QUERY_RESULT, &ts.second.startTime );
            if (ts.second.ended)
            {
                            glGetQueryObjectui64v( ts.second.endQuery, GL_QUERY_RESULT, &ts.second.endTime );
                            auto time = ts.second.endTime - ts.second.startTime;
                            if (ts.second.first)
                                ts.second.avTime = time;
                            else if ((ticks > debug_time - 1 && ticks < 2*debug_time)||(time > 0))
                                ts.second.avTime = av_mix*ts.second.avTime + (1 - av_mix)*time;
                            ts.second.first = false;
                            if (need_debug && (ticks % debug_time == 0))
                            {
                                debug("%s %d usec\n",ts.first.c_str(), ts.second.avTime/1000);
                            }
            }
            else
            {
                logerr("timestamp %s started but not ended",ts.first.c_str());
            }
        }
        else if (ts.second.ended)
        {
            glGetQueryObjectui64v( ts.second.endQuery, GL_QUERY_RESULT, &ts.second.endTime );
            logerr("timestamp %s ended but not started",ts.first.c_str());
        }
        ts.second.started = false;
        ts.second.ended = false;
    }
}