#pragma once
#include <vector>

struct BitVector
{
public:
    int size()
    {
        return values.size();
    }
    bool get(int n)
    {
        return (n >=0 && n < cur_size) ? get_unsafe(n) : default_value;
    }
    bool get_unsafe(int n)
    {
        return values[n];
    }
    void set(int n, bool val)
    {
        if (n>=0 && n < cur_size)
            set_unsafe(n, val);
    }
    void set_unsafe(int n, bool val)
    {
        values[n] = val;
    }
    void set_defaut(bool b)
    {
        default_value = b;
    }
    void push_back(bool value = false)
    {
        values.push_back(value);
        cur_size++;
    }
    void resize(int new_size, bool value = false)
    {
        cur_size = new_size;
        values = std::vector<bool>(cur_size, value);
    }

    std::vector<bool> values;
    bool default_value = false;
    int cur_size = 0;
};