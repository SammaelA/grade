#pragma once
#include <vector>
#include <string>
#include <glm/glm.hpp>

struct Block;
class BlkManager;
struct Block
{
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
        MAT4,
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
            glm::mat4 m4;
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
            else if (type == ValueType::MAT4)
                m4 = v.m4;
            else if (type == ValueType::STRING)
                s = v.s;
            else if (type == ValueType::BLOCK)
                bl = v.bl;
            else if (type == ValueType::ARRAY)
                a = v.a;
        }
        void copy(const Value &v)
        {
            clear();
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
            else if (type == ValueType::MAT4)
                m4 = v.m4;
            else if (type == ValueType::STRING)
            {
                s = new std::string();
                if (v.s)
                    *s = *v.s;
            }
            else if (type == ValueType::BLOCK)
            {
                bl = new Block();
                if (v.bl)
                    bl->copy(v.bl);
            }
            else if (type == ValueType::ARRAY)
            {
                a = new DataArray();
                if (v.a)
                    *a = *v.a;
            }
        }
        ~Value()
        {

        }
        void clear();
    };

    struct DataArray
    {
        ValueType type = EMPTY;
        std::vector<Value> values;
    };
    friend class BlkManager;

    int size();
    void clear();
    void copy(Block *b);
    ~Block() {clear();}
    bool has_tag(const std::string &name);
    int get_id(const std::string &name);
    int get_next_id(const std::string &name, int pos);
    std::string get_name(int id);
    ValueType get_type(int id);
    ValueType get_type(const std::string &name);

    int get_bool(int id, bool base_val = false);
    int get_int(int id, int base_val = 0);
    double get_double(int id, double base_val = 0);
    glm::vec2 get_vec2(int id, glm::vec2 base_val = glm::vec2(0, 0));
    glm::vec3 get_vec3(int id, glm::vec3 base_val = glm::vec3(0, 0, 0));
    glm::vec4 get_vec4(int id, glm::vec4 base_val = glm::vec4(0, 0, 0, 0));
    glm::mat4 get_mat4(int id, glm::mat4 base_val = glm::mat4(1.0f));
    std::string get_string(int id, std::string base_val = "");
    Block *get_block(int id);
    bool get_arr(int id, std::vector<double> &values, bool replace = false);
    bool get_arr(int id, std::vector<float> &values, bool replace = false);
    bool get_arr(int id, std::vector<int> &values, bool replace = false);

    int get_bool(const std::string name, bool base_val = false);
    int get_int(const std::string name, int base_val = 0);
    double get_double(const std::string name, double base_val = 0);
    glm::vec2 get_vec2(const std::string name, glm::vec2 base_val = glm::vec2(0, 0));
    glm::vec3 get_vec3(const std::string name, glm::vec3 base_val = glm::vec3(0, 0, 0));
    glm::vec4 get_vec4(const std::string name, glm::vec4 base_val = glm::vec4(0, 0, 0, 0));
    glm::mat4 get_mat4(const std::string name, glm::mat4 base_val = glm::mat4(1.0f));
    std::string get_string(const std::string name, std::string base_val = "");
    Block *get_block(std::string name);
    bool get_arr(const std::string name, std::vector<double> &values, bool replace = false);
    bool get_arr(const std::string name, std::vector<float> &values, bool replace = false);
    bool get_arr(const std::string name, std::vector<int> &values, bool replace = false);

    void add_bool(const std::string name, bool base_val = false);
    void add_int(const std::string name, int base_val = 0);
    void add_double(const std::string name, double base_val = 0);
    void add_vec2(const std::string name, glm::vec2 base_val = glm::vec2(0, 0));
    void add_vec3(const std::string name, glm::vec3 base_val = glm::vec3(0, 0, 0));
    void add_vec4(const std::string name, glm::vec4 base_val = glm::vec4(0, 0, 0, 0));
    void add_mat4(const std::string name, glm::mat4 base_val = glm::mat4(1.0f));
    void add_string(const std::string name, std::string base_val = "");
    void add_block(const std::string name, Block *bl);
    void add_arr(const std::string name, std::vector<double> &values);
    void add_arr(const std::string name, std::vector<float> &values);
    void add_arr(const std::string name, std::vector<int> &values);

    void set_bool(const std::string name, bool base_val = false);
    void set_int(const std::string name, int base_val = 0);
    void set_double(const std::string name, double base_val = 0);
    void set_vec2(const std::string name, glm::vec2 base_val = glm::vec2(0, 0));
    void set_vec3(const std::string name, glm::vec3 base_val = glm::vec3(0, 0, 0));
    void set_vec4(const std::string name, glm::vec4 base_val = glm::vec4(0, 0, 0, 0));
    void set_mat4(const std::string name, glm::mat4 base_val = glm::mat4(1.0f));
    void set_string(const std::string name, std::string base_val = "");
    void set_block(const std::string name, Block *bl);
    void set_arr(const std::string name, std::vector<double> &values);
    void set_arr(const std::string name, std::vector<float> &values);
    void set_arr(const std::string name, std::vector<int> &values);

    void add_value(const std::string &name, const Value &value);
    void set_value(const std::string &name, const Value &value);

    void add_detalization(Block &det);

    std::vector<std::string> names;
    std::vector<Value> values;
};

class BlkManager
{
public:
    bool load_block_from_string(std::string &str, Block &b);
    bool load_block_from_file(std::string path, Block &b);
    void save_block_to_string(std::string &str, Block &b);
    void save_block_to_file(std::string path, Block &b);
};

extern std::string base_blk_path;