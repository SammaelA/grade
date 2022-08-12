#pragma once
#include <string>
#include <cstdio>
#include <vector>
#include <list>
#include <exception>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include "common_utils/utility.h"

struct PackedBranch;
struct PackedLeaf;
struct PackedJoint;
struct BranchCatalogue;
struct InstancedBranch;
class Texture;
class TextureAtlas;
struct InstanceDataArrays;
struct Billboard;
struct BillboardData;
struct BillboardCloudData;
struct GrovePacked;
struct InstancingReadyModel;
struct Impostor;
struct ImpostorsData;
namespace saver
{
    extern std::string textures_path;
    void set_textures_path (std::string &tp);
    bool save(FILE *f, std::string &t);
    bool load(FILE *f, std::string &t);

    bool save(FILE *f, std::vector<GLfloat> &t);
    bool load(FILE *f, std::vector<GLfloat> &t);

    bool save(FILE *f, std::vector<PackedBranch> &t);
    bool load(FILE *f, std::vector<PackedBranch> &t);

    bool save(FILE *f, PackedLeaf &t);
    bool load(FILE *f, PackedLeaf &t);

    bool save(FILE *f, PackedJoint &t);
    bool load(FILE *f, PackedJoint &t);

    bool save(FILE *f, PackedBranch &t);
    bool load(FILE *f, PackedBranch &t);

    bool save(FILE *f, BranchCatalogue &t);
    bool load(FILE *f, BranchCatalogue &t);

    bool save(FILE *f, InstancedBranch &t);
    bool load(FILE *f, InstancedBranch &t);

    bool save(FILE *f, Texture &t);
    bool load(FILE *f, Texture &t);

    bool save(FILE *f, TextureAtlas &t);
    bool load(FILE *f, TextureAtlas &t);

    bool save(FILE *f, InstanceDataArrays &t);
    bool load(FILE *f, InstanceDataArrays &t);

    bool save(FILE *f, Billboard &t);
    bool load(FILE *f, Billboard &t);

    bool save(FILE *f, BillboardData &t);
    bool load(FILE *f, BillboardData &t);

    bool save(FILE *f, BillboardCloudData &t);
    bool load(FILE *f, BillboardCloudData &t);

    bool save(FILE *f, GrovePacked &t);
    bool load(FILE *f, GrovePacked &t);

    bool save(FILE *f, ImpostorsData &t);
    bool load(FILE *f, ImpostorsData &t);

    bool save(FILE *f, Impostor &t);
    bool load(FILE *f, Impostor &t);

    template <typename T>
    bool save(FILE *f, T &t)
    {
        try
        {
            int c = fwrite(&t, sizeof(T), 1, f);
            return c == 1;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool load(FILE *f, T &t)
    {
        try
        {
            int c = fread(&t, sizeof(T), 1, f);
            return c == 1;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T1,typename T2>
    bool save(FILE *f, std::pair<T1,T2> &t)
    {
        return save(f,t.first) && save(f,t.second);
    }
    template <typename T1,typename T2>
    bool load(FILE *f, std::pair<T1,T2> &t)
    {
        return load(f,t.first) && load(f,t.second);
    }
    template <typename T>
    bool savev(FILE *f, std::vector<T> &t)
    {
        try
        {
            if (t.size() > USHRT_MAX)
                throw std::out_of_range("vector is too long. Use savelv() instead of savev()");
            unsigned short k = t.size();
            bool status = true;
            status = status && save(f, k);
            for (int i = 0; i < k; i++)
            {
                status = status && save(f, t[i]);
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool loadv(FILE *f, std::vector<T> &t)
    {
        try
        {
            unsigned short k = 0;
            bool status = true;
            t.clear();
            status = status && load(f, k);
            t.resize(k);
            for (int i = 0; i < k; i++)
            {
                status = status && load(f, t[i]);
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool savelv(FILE *f, std::vector<T> &t)
    {
        try
        {
            if (t.size() > UINT_MAX)
                throw std::out_of_range("vector is too long. You must be joking...");
            unsigned k = t.size();
            bool status = true;
            status = status && save(f, k);
            for (int i = 0; i < k; i++)
            {
                status = status && save(f, t[i]);
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool loadlv(FILE *f, std::vector<T> &t)
    {
        try
        {
            unsigned k = 0;
            bool status = true;
            t.clear();
            status = status && load(f, k);
            t.resize(k);
            for (int i = 0; i < k; i++)
            {
                status = status && load(f, t[i]);
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool savesv(FILE *f, std::vector<T> &t)
    {
        try
        {
            short k = t.size();
            int c = fwrite(&k, sizeof(short), 1, f);
            bool ok = (c == 1);
            if (k > 0)
            {
                c = fwrite(t.data(), sizeof(T), k, f);
                ok = ok && (c == k);
            }
            return ok;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool loadsv(FILE *f, std::vector<T> &t)
    {
        try
        {
            short k;
            int c = fread(&k, sizeof(short), 1, f);
            bool ok = (c == 1);
            if (k > 0 && ok)
            {
                T *ptr = new T[k];
                c = fread(ptr, sizeof(T), k, f);
                ok = ok && (c == k);
                t.clear();
                t.assign(ptr, ptr + c);
                delete[] ptr;
            }
            return ok;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
        template <typename T>
    bool savelst(FILE *f, std::list<T> &t)
    {
        try
        {
            if (t.size() > USHRT_MAX)
                throw std::out_of_range("list is too long. Use savellst() instead of savelst()");
            unsigned short k = t.size();
            bool status = true;
            status = status && save(f, k);
            for (T &val : t)
            {
                status = status && save(f, val);
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool loadlst(FILE *f, std::list<T> &t)
    {
        try
        {
            unsigned short k = 0;
            bool status = true;
            t.clear();
            status = status && load(f, k);
            for (int i = 0; i < k; i++)
            {
                t.emplace_back();
                status = status && load(f, t.back());
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool savellst(FILE *f, std::vector<T> &t)
    {
        try
        {
            if (t.size() > UINT_MAX)
                throw std::out_of_range("list is too long. You must be joking...");
            unsigned k = t.size();
            bool status = true;
            status = status && save(f, k);
            for (T &val : t)
            {
                status = status && save(f, val);
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
    template <typename T>
    bool loadllst(FILE *f, std::list<T> &t)
    {
        try
        {
            unsigned k = 0;
            bool status = true;
            t.clear();
            status = status && load(f, k);
            for (int i = 0; i < k; i++)
            {
                t.emplace_back();
                status = status && load(f, t.back());
            }
            return status;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
    }
} // namespace saver