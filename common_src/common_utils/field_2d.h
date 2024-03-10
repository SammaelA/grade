#pragma once
#include "common_utils/LiteMath_ext.h"
#include <functional>

class Field_2d
{
public:
    Field_2d(glm::vec3 pos, glm::vec2 size, float cell_size){ create(pos, size, cell_size);}
    Field_2d(glm::vec3 pos, int w, int h){ create(pos, w, h);}
    Field_2d() {data = nullptr;}
    Field_2d& operator=(Field_2d &&s)
    {
      w = s.w;
      h = s.h;
      min_val = s.min_val;
      max_val = s.max_val;
      base_val = s.base_val;
      pos = s.pos;
      size = s.size;
      cell_size = s.cell_size;

      if (data)
        delete[] data;
      data = s.data;
      s.data = nullptr;
      return *this;
    }
    ~Field_2d();
    void create(glm::vec3 pos, glm::vec2 size, float cell_size);
    void create(glm::vec3 pos, int w, int h);
    float get_bilinear(glm::vec3 pos) const;
    float get_bilinear(glm::vec2 pos) const;
    void set(glm::vec3 pos, float val);
    void add(Field_2d &field, bool same_size_expected = false);
    void sub(Field_2d &field, bool same_size_expected = false);
    void mul(Field_2d &field, bool same_size_expected = false);
    void div(Field_2d &field, bool same_size_expected = false);
    void fill_const(float val);
    void fill_func(std::function<float(glm::vec2 &)> filler);
    void fill_func(std::function<float(glm::vec2 &, float )> filler);
    void read_func(std::function<void(glm::vec2 &, float )> reader) const;
    void fill_perlin(float base, float min, float max, glm::ivec2 sh = glm::ivec2(0,0));
    glm::vec2 get_range() const{return glm::vec2(min_val,max_val);}
    glm::vec2 get_grad_bilinear(glm::vec3 pos) const;
    glm::vec4 get_borders() const;//[x0,y0] - [x1,y1]
    void get_min_max_imprecise(glm::vec2 from, glm::vec2 to, float *min_v, float *max_v, 
                               glm::vec2 *min_pos = nullptr,
                               glm::vec2 *max_pos = nullptr) const;//real_min >= min_v, real_max <= max_v 
    void print() const;
    glm::vec3 get_pos() const {return pos;}
    glm::vec2 get_size() const {return size;}
    glm::ivec2 get_grid_size() const {return glm::ivec2(w,h);}
    float get_cell_size() const {return cell_size;}

    float get(int x, int y) const;
    glm::vec2 get_grad(int x, int y) const;
    void set(int x, int y, float val);
    void set_safe(int x, int y, float val);
    void add(int x, int y, float val);
    void add_safe(int x, int y, float val);
    int w = 0, h = 0;
    float min_val = 1e10, max_val = -1e10;
    glm::vec3 pos = glm::vec3(0,0,0);
    glm::vec2 size = glm::vec2(0,0);
    float cell_size = 1;
    float *data = nullptr;
    float base_val = 0;
};