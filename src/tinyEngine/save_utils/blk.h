#pragma once
#include <vector>
#include <string>
#include <glm/glm.hpp>

struct Block;
struct DataArray;

enum ValueType
{
    EMPTY,
    BOOL,
    INT,
    DOUBLE,
    VEC2,
    VEC3,
    VEC4,
    STRING,
    BLOCK,
    ARRAY
};

struct Value
{
    ValueType type;
    union
    {
        bool b;
        long i;
        double d;
        glm::vec2 v2;
        glm::vec3 v3;
        glm::vec4 v4;
        std::string *s;
        Block *bl;
        DataArray *a;
    };
    Value()
    {
        type = EMPTY;
    }
    Value(const Value &v)
    {
        type = v.type;
        if (type == ValueType::INT)
            i = v.i;
        else if (type == ValueType::BOOL)
            b = v.b;
        else if (type == ValueType::DOUBLE)
            d = v.d;
        else if (type == ValueType::VEC2)
            v2 = v.v2;
        else if (type == ValueType::VEC3)
            v3 = v.v3;
        else if (type == ValueType::VEC4)
            v4 = v.v4;
        else if (type == ValueType::STRING)
            s = v.s;
        else if (type == ValueType::BLOCK)
            bl = v.bl;
        else if (type == ValueType::ARRAY)
            a = v.a;
    }
    ~Value()
    {
        //if (type == ValueType::BLOCK && b)
        //    delete b;
        //else if (type == ValueType::ARRAY && a)
        //    delete a;   
    }
    void clear();
    /*{
        if (type == ValueType::BLOCK && b)
            delete b;
        else if (type == ValueType::ARRAY && a)
            delete a;  
    }*/
};

struct DataArray
{
    ValueType type = EMPTY;
    std::vector<Value> values;
};

struct Block
{
    int size();
    void clear();

    bool has_tag(std::string name);
    int get_id(std::string name);
    int get_next_id(std::string name, int pos);
    std::string  get_name(int id);
    ValueType get_type(int id);
    ValueType get_type(std::string name);
    
    int get_bool(int id, bool base_val = false);
    int get_int(int id, int base_val = 0);
    double get_double(int id, double base_val = 0);
    glm::vec2 get_vec2(int id, glm::vec2 base_val = glm::vec2(0,0));
    glm::vec3 get_vec3(int id, glm::vec3 base_val = glm::vec3(0,0,0));
    glm::vec4 get_vec4(int id, glm::vec4 base_val = glm::vec4(0,0,0,0));
    std::string get_string(int id, std::string base_val = "");
    Block *get_block(int id);
    bool get_arr(int id, std::vector<double> &values, bool replace = false);
    bool get_arr(int id, std::vector<float> &values, bool replace = false);
    bool get_arr(int id, std::vector<int> &values, bool replace = false);

    int get_bool(std::string name, bool base_val = false);
    int get_int(std::string name, int base_val = 0);
    double get_double(std::string name, double base_val = 0);
    glm::vec2 get_vec2(std::string name, glm::vec2 base_val = glm::vec2(0,0));
    glm::vec3 get_vec3(std::string name, glm::vec3 base_val = glm::vec3(0,0,0));
    glm::vec4 get_vec4(std::string name, glm::vec4 base_val = glm::vec4(0,0,0,0));
    std::string get_string(std::string name, std::string base_val = "");
    Block *get_block(std::string name);
    bool get_arr(std::string name, std::vector<double> &values, bool replace = false);
    bool get_arr(std::string name, std::vector<float> &values, bool replace = false);
    bool get_arr(std::string name, std::vector<int> &values, bool replace = false);

    void add_bool(std::string name, bool base_val = false);
    void add_int(std::string name, int base_val = 0);
    void add_double(std::string name, double base_val = 0);
    void add_vec2(std::string name, glm::vec2 base_val = glm::vec2(0,0));
    void add_vec3(std::string name, glm::vec3 base_val = glm::vec3(0,0,0));
    void add_vec4(std::string name, glm::vec4 base_val = glm::vec4(0,0,0,0));
    void add_string(std::string name, std::string base_val = "");
    void add_block(std::string name, Block *bl);
    void add_arr(std::string name, std::vector<double> &values);
    void add_arr(std::string name, std::vector<float> &values);
    void add_arr(std::string name, std::vector<int> &values);

    std::vector<std::string> names;
    std::vector<Value> values;
};

class BlkManager
{
    public:
    bool load_block_from_file(std::string path, Block &b);
    void save_block_to_file(std::string path, Block &b);
};