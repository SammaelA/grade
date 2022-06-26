#include "common_utils/field_2d.h"
#include "common_utils/utility.h"
#include "graphics_utils/texture_manager.h"
#include "perlin.h"
#include <vector>
    void Field_2d::create(glm::vec3 pos, glm::vec2 size, float cell_size)
    {
      create(pos, ceil(size.x/cell_size), ceil(size.y/cell_size));
      this->cell_size = cell_size;
      this->pos = pos;
      this->size = size;
    }
    void Field_2d::create(glm::vec3 _pos, int _w, int _h)
    {
      w = _w;
      h = _h;
      pos = _pos;
      if (data)
        delete[] data;
      data = new float[(2*w + 1)*(2*h + 1)];
    }
    Field_2d::~Field_2d()
    {
      if (data)
        delete[] data;
    }
    float Field_2d::get(int x, int y) const
    {
        if (x >= -w && x <= w && y >= -h && y <= h)
        {
            return data[(2*w + 1)*(y + h) + x + w];
        }
        else
            return base_val;
    }
    float Field_2d::get_bilinear(glm::vec2 position) const
    {
        if (!data)
            return base_val;
        glm::vec2 rp = glm::vec2(position.x - pos.x, position.y - pos.z)/cell_size;
        glm::ivec2 ps = rp;
        float dx = rp.x - ps.x;
        float dy = rp.y - ps.y;
        return (dx*get(ps.x, ps.y) + (1 - dx)*get(ps.x + 1, ps.y))*(1 - dy) + 
               (dx*get(ps.x, ps.y + 1) + (1 - dx)*get(ps.x + 1, ps.y + 1))*dy;
        
    }

    float Field_2d::get_bilinear(glm::vec3 position) const
    {
        return get_bilinear(glm::vec2(position.x, position.z));
    }

    void Field_2d::set(glm::vec3 position, float val)
    {
        if (!data)
            return;
        glm::vec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;
        glm::ivec2 ps = rp;
        set(rp.x,rp.y,val);
    }
    void Field_2d::fill_const(float val)
    {
        base_val = val;
        if (!data)
            return;
        for (int i=0;i<(2*w + 1)*(2*h + 1);i++)
            data[i] = val;
    }
    void Field_2d::fill_func(std::function<float(glm::vec2 &)> filler)
    {
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                set(i,j,filler(ps));
            }
        }
    }
    void Field_2d::fill_func(std::function<float(glm::vec2 &, float )> filler)
    {
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                set(i,j,filler(ps, get(i,j)));
            }
        }
    }
    void Field_2d::read_func(std::function<void(glm::vec2 &, float )> reader) const
    {
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                reader(ps, get(i,j));
            }
        }
    }
    void Field_2d::fill_perlin(float base, float min, float max, glm::ivec2 sh)
    {
        base_val = base;
        min_val = base;
        max_val = base;

        if (!data)
            return;
        for (int i = -w; i <= w; i++)
        {
            for (int j = -h; j <= h; j++)
            {
                float height = min + 0.5*(max - min)*(1 + f_perlin((float)(i + w+sh.x)/(2*w + 1), (float)(j + h+sh.y)/(2*h + 1)));
                set(i,j,height);
            }
        }
    }
    glm::vec2 Field_2d::get_grad(int x, int y) const
    {
        std::vector<std::vector<float>> sobel_x = {{-1,0,1},
                                                   {-2,0,2},
                                                   {-1,0,1}};
        std::vector<std::vector<float>> sobel_y = {{-1,-2,-1},
                                                   {0,0,0},
                                                   {1,2,1}};
        glm::vec2 g = glm::vec2(0,0);
        for (int i=-1;i<=1;i++)
        {
            for (int j=-1;j<=1;j++)
            {
                float h = get(x+i,y+j);
                g.x += sobel_x[i+1][j+1]*h/cell_size;
                g.y += sobel_y[i+1][j+1]*h/cell_size;
            }
        }
        return g;
    }
    glm::vec2 Field_2d::get_grad_bilinear(glm::vec3 position) const
    {
        if (!data)
            return glm::vec2(0,0);
        glm::vec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;
        glm::ivec2 ps = rp;
        float dx = rp.x - ps.x;
        float dy = rp.y - ps.y;
        return (dx*get_grad(ps.x, ps.y) + (1 - dx)*get_grad(ps.x + 1, ps.y))*(1 - dy) + 
               (dx*get_grad(ps.x, ps.y + 1) + (1 - dx)*get_grad(ps.x + 1, ps.y + 1))*dy;
    }
    void Field_2d::set(int x, int y, float val)
    {
        data[(2*w + 1)*(y + h) + x + w] = val;
        min_val = MIN(min_val,val);
        max_val = MAX(max_val,val);
    }
    void Field_2d::set_safe(int x, int y, float val)
    {
        if (x >= -w && x <= w && y >= -h && y <= h)
        {
            set(x,y,val);
        }
    }
        void Field_2d::add(int x, int y, float val)
    {
        data[(2*w + 1)*(y + h) + x + w] += val;
        min_val = MIN(min_val,val);
        max_val = MAX(max_val,val);
    }
    void Field_2d::add_safe(int x, int y, float val)
    {
        if (x >= -w && x <= w && y >= -h && y <= h)
        {
            add(x,y,val);
        }
    }
    glm::vec4 Field_2d::get_borders() const
    {
        return glm::vec4(pos.x - size.x, pos.z - size.y, pos.x + size.x, pos.z + size.y);
    }
    void Field_2d::print() const
    {
        debug("field 2d\n-------\n");
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                debug("(%f %f) = %f %d %d|",ps.x,ps.y,get(i,j),i,j);
            }
            debugnl();
        }
    }
    void Field_2d::save_as_image(std::string name, float mnv, float mxv) const
    {
        if (mnv == mxv)
        {
            mnv = min_val;
            mxv = max_val;
        }
        if (mnv <= mxv)
        {
            mnv = 0;
            mxv = 1;
        }

        if (!data)
            return;
        unsigned char *image_data = new unsigned char[3*(2*w + 1)*(2*h+1)];
        for (int i=-h;i<=h;i++)
        {
            for (int j=-w;j<=w;j++)
            {
                float val = (get(i,j) - mnv)/(mxv - mnv);
                int pos = (i+h)*(2*w + 1) + j + w;
                image_data[3*pos]   = 255*val;
                image_data[3*pos+1] = 255*val;
                image_data[3*pos+2] = 255*val;
            }
        }
        textureManager.save_bmp_raw(image_data, 2*w+1, 2*h+1, 3, name);

        delete[] image_data;
    }
    
    Texture Field_2d::save_as_texture_RGBA8(float mnv, float mxv) const
    {
      if (mnv == mxv)
        {
            mnv = min_val;
            mxv = max_val;
        }
        if (mnv <= mxv)
        {
            mnv = 0;
            mxv = 1;
        }

        if (!data)
            return textureManager.empty();
        unsigned char *image_data = new unsigned char[4*(2*w + 1)*(2*h+1)];
        for (int i=-h;i<=h;i++)
        {
            for (int j=-w;j<=w;j++)
            {
                float val = (get(i,j) - mnv)/(mxv - mnv);
                int pos = (j+w)*(2*h + 1) + i + h;
                image_data[4*pos]   = 255*val;
                image_data[4*pos+1] = 255*val;
                image_data[4*pos+2] = 255*val;
                image_data[4*pos+3] = 255;
            }
        }

        Texture t = textureManager.create_texture(2*w+1, 2*h+1);
        t = textureManager.load_unnamed(t, image_data);

        delete[] image_data;
        return t;
    }

    void Field_2d::get_min_max_imprecise(glm::vec2 from, glm::vec2 to, float *min_v, float *max_v, 
                                         glm::vec2 *min_pos, glm::vec2 *max_pos) const
    {
        if (!data)
        {
            if (min_v)
                *min_v = base_val;
            if (max_v)
                *max_v = base_val;
            if (min_pos)
                *min_pos = from;
            if (max_pos)
                *max_pos = from;
            return;
        }
        glm::vec2 mnp = (from - glm::vec2(pos.x,pos.y))/cell_size;
        glm::vec2 mxp = (to - glm::vec2(pos.x,pos.y))/cell_size;
        
        float mn = FLT_MAX;
        float mx = FLT_MIN;
        glm::ivec2 mn_pos = glm::ivec2(-1,-1);
        glm::ivec2 mx_pos = glm::ivec2(-1,-1);
        if (mnp.x < 0 || mnp.y < 0)
        {
            mn = base_val;
            mx = base_val;
            mn_pos = from;
            mx_pos = from;
        }
        else if (mxp.x >= w || mxp.y >= h)
        {
            mn = base_val;
            mx = base_val;
            mn_pos = glm::ivec2(w+1,h+1);
            mx_pos = glm::ivec2(w+1,h+1);
        }
        for (int x = MAX(floor(mnp.x), 0);x<MIN(ceil(mxp.x),w);x++)
        {
            for (int y = MAX(floor(mnp.y), 0);y<MIN(ceil(mxp.y),h);y++)
            {
                float val = get(x,y);
                if (val > mx)
                {
                    mx = val;
                    mx_pos = glm::ivec2(x,y);
                }
                if (val < mn)
                {
                    mn = val;
                    mn_pos = glm::ivec2(x,y);
                }
            }
        }

        if (min_v)
            *min_v = mn;
        if (max_v)
            *max_v = mx;
        if (min_pos)
        {
            if (mn_pos == glm::ivec2(-1,-1))
                *min_pos = from;
            else if (mn_pos == glm::ivec2(w + 1,h + 1))
                *min_pos = to;
            else 
                *min_pos = glm::vec2(pos.x,pos.y) + cell_size*glm::vec2(mn_pos);
        }
        if (max_pos)
        {
            if (mx_pos == glm::ivec2(-1,-1))
                *max_pos = from;
            else if (mx_pos == glm::ivec2(w + 1,h + 1))
                *max_pos = to;
            else 
                *max_pos = glm::vec2(pos.x,pos.y) + cell_size*glm::vec2(mx_pos);
        }
    }

    void Field_2d::add(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos == field.pos && size == field.size && w == field.w && h == field.h);
        if (!same_size && same_size_expected)
        {
            logerr("warning: Field_2d add same size and position for fields expected");
        }
        if (same_size)
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    set(i,j,get(i,j) + field.get(i,j));
                }
            }
        }
        else
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) + field.get_bilinear(ps));
                }
            }
        }
    }
    void Field_2d::sub(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos == field.pos && size == field.size && w == field.w && h == field.h);
        if (!same_size && same_size_expected)
        {
            logerr("warning: Field_2d sub same size and position for fields expected");
        }
        if (same_size)
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    set(i,j,get(i,j) - field.get(i,j));
                }
            }
        }
        else
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) - field.get_bilinear(ps));
                }
            }
        }
    }
    void Field_2d::mul(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos == field.pos && size == field.size && w == field.w && h == field.h);
        if (!same_size && same_size_expected)
        {
            logerr("warning: Field_2d mul same size and position for fields expected");
        }
        if (same_size)
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    set(i,j,get(i,j) * field.get(i,j));
                }
            }
        }
        else
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) * field.get_bilinear(ps));
                }
            }
        }
    }
    void Field_2d::div(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos == field.pos && size == field.size && w == field.w && h == field.h);
        if (!same_size && same_size_expected)
        {
            logerr("warning: Field_2d div same size and position for fields expected");
        }
        if (same_size)
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    set(i,j,get(i,j) / field.get(i,j));
                }
            }
        }
        else
        {
            for (int j=-h;j<=h;j++)
            {
                for (int i = -w;i<=w;i++)
                {
                    glm::vec2 ps = glm::vec2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) / field.get_bilinear(ps));
                }
            }
        }
    }