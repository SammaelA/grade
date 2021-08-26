#include "blk.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "../utility.h"
#include <boost/algorithm/string.hpp>

int cur_line = 0;

bool is_empty(const char c)
{
    if (c=='\n')
        cur_line++;
    return (c==' ' || c=='\n' || c=='\t');
}
bool is_div(const char c)
{
    return (c==',' || c==';' || c==':' || c=='=' || c=='{' || c=='}'|| c=='\'' || c=='\"');
}
std::string next_token(const char *data, int &pos)
{
    std::string res;
    if (!data || data[pos] == 0)
        res = "";
    else
    {
        while (is_empty(data[pos]))
            pos++;
        if (data[pos] == 0)
            res = "";
        else
        {
            if (is_div(data[pos]))
            {
                pos++;
                res = std::string(1, data[pos-1]);
            }
            else
            {
                const char *start = data + pos;
                int sz = 0;
                while (!is_div(data[pos]) && !is_empty(data[pos]) && data[pos] != 0)
                {
                    pos++;
                    sz++;
                }
                res = std::string(start, sz);
            }
        }
    }
    return res;
}

std::vector<std::pair<ValueType, std::string>> descs = {
    {EMPTY,"em"},
    {INT,"i"},
    {DOUBLE,"r"},
    {VEC2,"p2"},
    {VEC3,"p3"},
    {VEC4,"p4"},
    {STRING,"s"},
    {ARRAY,"arr"}
};
bool load_block(const char *data, int &cur_pos, Block &b);
bool read_array(const char *data, int &cur_pos, DataArray &a);
bool read_value(const char *data, int &cur_pos, Value &v)
{
    std::string token = next_token(data, cur_pos);
    //:<type> = <description> or { <block> }
    if (token == "{")
    {
        v.bl = new Block();
        v.type = ValueType::BLOCK;
        return load_block(data,cur_pos,*(v.bl));
    }
    else if (token == ":")
    {//simple value or array
        std::string type = next_token(data, cur_pos);
        if (type == "tag")
        {
            v.type = ValueType::EMPTY;
            return true;
        }
        std::string eq = next_token(data, cur_pos);
        if (eq != "=")
        {
            logerr("line %d expected = after value type",cur_line);
            v.type = ValueType::EMPTY;
            return false;
        }
        if (type == "b")
        {
            std::string val = next_token(data, cur_pos);
            v.type = ValueType::BOOL;
            v.b = boost::iequals(val,"true");
        }
        else if (type == "i")
        {
            std::string val = next_token(data, cur_pos);
            v.type = ValueType::INT;
            v.i = std::stol(val);
        }
        else if (type == "r")
        {
            std::string val = next_token(data, cur_pos);
            v.type = ValueType::DOUBLE;
            v.d = std::stod(val);
        }
        else if (type == "p2")
        {
            v.type = ValueType::VEC2;
            std::string val;
            bool ok = true;
            v.v2 = glm::vec2(0,0);

            val = next_token(data, cur_pos);
            v.v2.x = std::stod(val);

            val = next_token(data, cur_pos);
            ok = ok && (val == ",");

            if (ok)
            {
               val = next_token(data, cur_pos);
               v.v2.y = std::stod(val);
            }
            if (!ok)
            {
                logerr("line %d wrong description of vector",cur_line);
                v.type = ValueType::EMPTY;
                return false;
            }
        }
        else if (type == "p3")
        {
            v.type = ValueType::VEC3;
            std::string val;
            bool ok = true;
            v.v3 = glm::vec3(0,0,0);

            val = next_token(data, cur_pos);
            v.v3.x = std::stod(val);

            val = next_token(data, cur_pos);
            ok = ok && (val == ",");
            if (ok)
            {
               val = next_token(data, cur_pos);
               v.v3.y = std::stod(val);
            }

            val = next_token(data, cur_pos);
            ok = ok && (val == ",");
            if (ok)
            {
               val = next_token(data, cur_pos);
               v.v3.z = std::stod(val);
            }
            if (!ok)
            {
                logerr("line %d wrong description of vector",cur_line);
                v.type = ValueType::EMPTY;
                return false;
            }
        }
        else if (type == "p4")
        {
            v.type = ValueType::VEC4;
            std::string val;
            bool ok = true;
            v.v4 = glm::vec4(0,0,0,0);

            val = next_token(data, cur_pos);
            v.v4.x = std::stod(val);

            val = next_token(data, cur_pos);
            ok = ok && (val == ",");
            if (ok)
            {
               val = next_token(data, cur_pos);
               v.v4.y = std::stod(val);
            }

            val = next_token(data, cur_pos);
            ok = ok && (val == ",");
            if (ok)
            {
               val = next_token(data, cur_pos);
               v.v4.z = std::stod(val);
            }

            val = next_token(data, cur_pos);
            ok = ok && (val == ",");
            if (ok)
            {
               val = next_token(data, cur_pos);
               v.v4.w = std::stod(val);
            }
            if (!ok)
            {
                logerr("line %d wrong description of vector",cur_line);
                v.type = ValueType::EMPTY;
                return false;
            }
        }
        else if (type == "s")
        {
            std::string par = next_token(data, cur_pos);
            if (par == "\"")
            {
                int len = 0;
                const char *start = data+cur_pos;
                while (data[cur_pos] != 0 && data[cur_pos] != '\"')
                {
                    if (data[cur_pos] == '\n')
                        cur_line++;
                    len++;
                    cur_pos++;
                }
                if (data[cur_pos] == 0)
                {
                    v.type = EMPTY;
                    logerr("line %d expected \" at the end of a string", cur_line);
                    return false;
                }
                else if (data[cur_pos] == '\"')
                {
                    cur_pos++;
                    v.type = ValueType::STRING;
                    v.s =  new std::string(start, len);
                }
            }
        }
        else if (type == "arr")
        {
            v.type = ValueType::ARRAY;
            v.a = new DataArray();
            return read_array(data, cur_pos, *(v.a));
        }

        return true;
    }
    else
    {
        logerr("line %d expected : or { after value/block name", cur_line);
        v.type = ValueType::EMPTY;
        return false;
    }
}
bool read_array(const char *data, int &cur_pos, DataArray &a)
{
    std::string token = next_token(data, cur_pos);
    //{ <value>, <value>, ... <value>}
    if (token == "{")
    {
        bool ok = true;
        while (ok)
        {
            Value val;
            std::string tok = next_token(data, cur_pos);
            if (tok == "}")
            {
                //empty array
                return true;
            }
            val.type = ValueType::DOUBLE;
            val.d = std::stod(tok);
            a.values.push_back(val);

            tok = next_token(data, cur_pos);
            ok = ok && (tok == ",");
            if (tok == "}")
            {
                a.type = ValueType::DOUBLE;
                return true;
            }
        }
        logerr("line %d expected } at the end of array", cur_line);
        return false;
    }
    else
    {
        logerr("line %d expected { at the start of array", cur_line);
        return false;
    }
}
bool load_block(const char *data, int &cur_pos, Block &b)
{
    bool correct = true;
    while (correct)
    {
        std::string token = next_token(data, cur_pos);
        if (token == "}")
        {
            //block closed correctly 
            return correct;
        }
        else if (token == "")
        {
            //end of file
            logerr("line %d block loader reached end of file, } expected",cur_line);
            return false;
        }
        else
        {
            //next value
            b.names.push_back(token);
            b.values.emplace_back();
            correct = correct && read_value(data, cur_pos, b.values.back());
        }
    }
}
bool BlkManager::load_block_from_file(std::string path, Block &b)
{
    std::fstream f(path);
    std::stringstream iss;
    if (f.fail())
    {
        logerr("unable to load file %s",path.c_str());
        return false;
    }
    iss << f.rdbuf();
    std::string entireFile = iss.str();
    cur_line = 0;
    int cur_pos = 0;
    const char *data = entireFile.c_str();
    std::string token = next_token(data, cur_pos);
    if (token == "{")
    {
        return load_block(data,cur_pos,b);
    }
    else
    {
        logerr("global block should start with {");
        return false;
    }
}

int Block::size()
{
    return names.size();
}

int Block::get_id(std::string name)
{
    return get_next_id(name,0);
}
int Block::get_next_id(std::string name, int pos)
{
    for (int i=pos;i<names.size();i++)
    {
        if (names[i] == name)
            return i;
    }
    return -1;
}
ValueType Block::get_type(int id)
{
    return (id >=0 && id <size()) ? values[id].type : ValueType::EMPTY;
}
ValueType Block::get_type(std::string name)
{
    return get_type(get_id(name));
}

int Block::get_bool(int id, bool base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::BOOL) ? values[id].b : base_val;
}
int Block::get_int(int id, int base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::INT) ? values[id].i : base_val;
}
double Block::get_double(int id, double base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::DOUBLE) ? values[id].d : base_val;
}
glm::vec2 Block::get_vec2(int id, glm::vec2 base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::VEC2) ? values[id].v2 : base_val;
}
glm::vec3 Block::get_vec3(int id, glm::vec3 base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::VEC3) ? values[id].v3 : base_val;
}
glm::vec4 Block::get_vec4(int id, glm::vec4 base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::VEC4) ? values[id].v4 : base_val;
}
std::string Block::get_string(int id, std::string base_val)
{
    return (id >=0 && id <size() && values[id].type == ValueType::STRING && values[id].s) ? *(values[id].s) : base_val;
}
Block *Block::get_block(int id)
{
    return (id >=0 && id <size() && values[id].type == ValueType::BLOCK) ? values[id].bl : nullptr;
}
    bool Block::get_arr(int id, std::vector<double> &_values, bool replace)
    {
        if (id >=0 && id <size() && values[id].type == ValueType::ARRAY && values[id].a && 
           (values[id].a->type == DOUBLE))
        {
            if (replace)
                _values.clear();
            for (Value &v : values[id].a->values)
            {
                _values.push_back(v.d);
            }
            return true;
        }
        return false;
    }
    bool Block::get_arr(int id, std::vector<float> &_values, bool replace)
    {
        if (id >=0 && id <size() && values[id].type == ValueType::ARRAY && values[id].a && 
           (values[id].a->type == DOUBLE))
        {
            if (replace)
                _values.clear();
            for (Value &v : values[id].a->values)
            {
                _values.push_back(v.d);
            }
            return true;
        }
        return false;
    }
    bool Block::get_arr(int id, std::vector<int> &_values, bool replace)
    {
        if (id >=0 && id <size() && values[id].type == ValueType::ARRAY && values[id].a && 
           (values[id].a->type == DOUBLE))
        {
            if (replace)
                _values.clear();
            for (Value &v : values[id].a->values)
            {
                _values.push_back(v.d);
            }
            return true;
        }
        return false;
    }


int Block::get_bool(std::string name, bool base_val)
{
    return get_bool(get_id(name), base_val);
}
int Block::get_int(std::string name, int base_val)
{
    return get_int(get_id(name), base_val);
}
double Block::get_double(std::string name, double base_val)
{
    return get_double(get_id(name), base_val);
}
glm::vec2 Block::get_vec2(std::string name, glm::vec2 base_val)
{
    return get_vec2(get_id(name), base_val);
}
glm::vec3 Block::get_vec3(std::string name, glm::vec3 base_val)
{
    return get_vec3(get_id(name), base_val);
}
glm::vec4 Block::get_vec4(std::string name, glm::vec4 base_val)
{
    return get_vec4(get_id(name), base_val);
}
std::string Block::get_string(std::string name, std::string base_val)
{
    return get_string(get_id(name), base_val);
}
Block *Block::get_block(std::string name)
{
    return get_block(get_id(name));
}
    bool Block::get_arr(std::string name, std::vector<double> &_values, bool replace)
    {
        return get_arr(get_id(name), _values, replace);
    }
    bool Block::get_arr(std::string name, std::vector<float> &_values, bool replace)
    {
        return get_arr(get_id(name), _values, replace);
    }
    bool Block::get_arr(std::string name, std::vector<int> &_values, bool replace)
    {
        return get_arr(get_id(name), _values, replace);
    }

void save_value(std::string &str, Value &v);
void save_block(std::string &str, Block &b)
{
    str += "{\n";
    for (int i=0;i<b.size();i++)
    {
        str += b.names[i];
        save_value(str,b.values[i]);
        str += "\n";
    }
    str += "}";
}
void save_arr(std::string &str, DataArray &a)
{
    str += "{ ";
    if (a.type == ValueType::DOUBLE)
    {
        for (int i=0;i<a.values.size();i++)
        {
            str += std::to_string(a.values[i].d);
            if (i < a.values.size() -1)
                str += ", ";
        }
    }
    str += " }";
}
void save_value(std::string &str, Value &v)
{
    if (v.type == EMPTY)
    {
        str +=":tag";
    }
    else if (v.type == BOOL)
    {
        str+=":b = ";
        str+= v.b ? "true" : "false";
    }
    else if (v.type == INT)
    {
        str+=":i = ";
        str+= std::to_string(v.i);
    }
    else if (v.type == DOUBLE)
    {
        str+=":r = ";
        str+= std::to_string(v.d);
    }
    else if (v.type == VEC2)
    {
        str+=":p2 = ";
        str+= std::to_string(v.v2.x) + ", " + std::to_string(v.v2.y);
    }
    else if (v.type == VEC3)
    {
        str+=":p3 = ";
        str+= std::to_string(v.v3.x) + ", " + std::to_string(v.v3.y) + ", " + std::to_string(v.v3.z);
    }
    else if (v.type == VEC4)
    {
        str+=":p4 = ";
        str+= std::to_string(v.v4.x) + ", " + std::to_string(v.v4.y) + ", " + std::to_string(v.v4.z) +
              ", " + std::to_string(v.v4.w);
    }
    else if (v.type == STRING && v.s)
    {
        str+=":s = \"" + *(v.s) + "\"";
    }
    else if (v.type == ARRAY && v.a)
    {
        str+=":arr = ";
        save_arr(str,*(v.a));
    }
    else if (v.type == BLOCK && v.bl)
    {
        str+=" ";
        save_block(str,*(v.bl));
    }
}
void BlkManager::save_block_to_file(std::string path, Block &b)
{
    std::string input;

    save_block(input,b);

    std::ofstream out(path);
    out << input;
    out.close();
}

void Value::clear()
{
    if (type == ValueType::BLOCK && bl)
        delete bl;
    else if (type == ValueType::ARRAY && a)
        delete a;  
    else if (type == ValueType::STRING && s)
        delete s;
    
    type = EMPTY;
}
void Block::clear()
{
    for (int i = 0;i<size();i++)
    {
        values[i].clear();
    }
    values.clear();
    names.clear();
}
bool Block::has_tag(std::string name)
{
    int id = get_id(name);
    return id >=0 && (get_type(id) == ValueType::EMPTY);
}

    void Block::add_bool(std::string name, bool base_val)
    {
        Value val;
        val.type = ValueType::BOOL;
        val.b = base_val;
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_int(std::string name, int base_val)
    {
        Value val;
        val.type = ValueType::INT;
        val.i = base_val;
        values.push_back(val);
        names.push_back(name);    
    }
    void Block::add_double(std::string name, double base_val)
    {
        Value val;
        val.type = ValueType::DOUBLE;
        val.d = base_val;
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_vec2(std::string name, glm::vec2 base_val)
    {
        Value val;
        val.type = ValueType::VEC2;
        val.v2 = base_val;
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_vec3(std::string name, glm::vec3 base_val)
    {
        Value val;
        val.type = ValueType::VEC3;
        val.v3 = base_val;
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_vec4(std::string name, glm::vec4 base_val)
    {
        Value val;
        val.type = ValueType::VEC4;
        val.v4 = base_val;
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_string(std::string name, std::string base_val)
    {
        Value val;
        val.type = ValueType::STRING;
        val.s = new std::string(base_val);
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_block(std::string name, Block *bl)
    {
        Value val;
        val.type = ValueType::BLOCK;
        val.bl = bl;
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_arr(std::string name, std::vector<double> &_values)
    {
        Value val;
        val.type = ValueType::ARRAY;
        val.a = new DataArray();
        val.a->type = ValueType::DOUBLE;
        for (double &d : _values)
        {
            Value av;
            av.type = ValueType::DOUBLE;
            av.d = d;
            val.a->values.push_back(av);
        }
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_arr(std::string name, std::vector<float> &_values)
    {
        Value val;
        val.type = ValueType::ARRAY;
        val.a = new DataArray();
        val.a->type = ValueType::DOUBLE;
        for (float &d : _values)
        {
            Value av;
            av.type = ValueType::DOUBLE;
            av.d = d;
            val.a->values.push_back(av);
        }
        values.push_back(val);
        names.push_back(name);
    }
    void Block::add_arr(std::string name, std::vector<int> &_values)
    {
        Value val;
        val.type = ValueType::ARRAY;
        val.a = new DataArray();
        val.a->type = ValueType::DOUBLE;
        for (int &d : _values)
        {
            Value av;
            av.type = ValueType::DOUBLE;
            av.d = d;
            val.a->values.push_back(av);
        }
        values.push_back(val);
        names.push_back(name);
    }
    std::string  Block::get_name(int id)
    {
        return (id >= 0 && id < names.size()) ? names[id] : "";
    }