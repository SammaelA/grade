#include "common_utils/field_2d.h"
#include "common_utils/utility.h"
#include "perlin.h"
#include <vector>
    void Field_2d::create(float3 pos, float2 size, float cell_size)
    {
      create(pos, ceil(size.x/cell_size), ceil(size.y/cell_size));
      this->cell_size = cell_size;
      this->pos = pos;
      this->size = size;
    }
    void Field_2d::create(float3 _pos, int _w, int _h)
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
    float Field_2d::get_bilinear(float2 position) const
    {
        if (!data)
            return base_val;
        float2 rp = float2(position.x - pos.x, position.y - pos.z)/cell_size;
        int2 ps = to_int2(rp);
        float dx = rp.x - ps.x;
        float dy = rp.y - ps.y;
        return (dx*get(ps.x, ps.y) + (1 - dx)*get(ps.x + 1, ps.y))*(1 - dy) + 
               (dx*get(ps.x, ps.y + 1) + (1 - dx)*get(ps.x + 1, ps.y + 1))*dy;
        
    }

    float Field_2d::get_bilinear(float3 position) const
    {
        return get_bilinear(float2(position.x, position.z));
    }

    void Field_2d::set(float3 position, float val)
    {
        if (!data)
            return;
        float2 rp = float2(position.x - pos.x, position.z - pos.z)/cell_size;
        int2 ps = to_int2(rp);
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
    void Field_2d::fill_func(std::function<float(float2 &)> filler)
    {
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                set(i,j,filler(ps));
            }
        }
    }
    void Field_2d::fill_func(std::function<float(float2 &, float )> filler)
    {
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                set(i,j,filler(ps, get(i,j)));
            }
        }
    }
    void Field_2d::read_func(std::function<void(float2 &, float )> reader) const
    {
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                reader(ps, get(i,j));
            }
        }
    }
    void Field_2d::fill_perlin(float base, float min, float max, int2 sh)
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
    float2 Field_2d::get_grad(int x, int y) const
    {
        std::vector<std::vector<float>> sobel_x = {{-1,0,1},
                                                   {-2,0,2},
                                                   {-1,0,1}};
        std::vector<std::vector<float>> sobel_y = {{-1,-2,-1},
                                                   {0,0,0},
                                                   {1,2,1}};
        float2 g = float2(0,0);
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
    float2 Field_2d::get_grad_bilinear(float3 position) const
    {
        if (!data)
            return float2(0,0);
        float2 rp = float2(position.x - pos.x, position.z - pos.z)/cell_size;
        int2 ps = to_int2(rp);
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
    float4 Field_2d::get_borders() const
    {
        return float4(pos.x - size.x, pos.z - size.y, pos.x + size.x, pos.z + size.y);
    }
    void Field_2d::print() const
    {
        debug("field 2d\n-------\n");
        for (int i = -w;i<=w;i++)
        {
            for (int j=-h;j<=h;j++)
            {
                float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                debug("(%f %f) = %f %d %d|",ps.x,ps.y,get(i,j),i,j);
            }
            debugnl();
        }
    }

    void Field_2d::get_min_max_imprecise(float2 from, float2 to, float *min_v, float *max_v, 
                                         float2 *min_pos, float2 *max_pos) const
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
        float2 mnp = (from - float2(pos.x,pos.y))/cell_size;
        float2 mxp = (to - float2(pos.x,pos.y))/cell_size;
        
        float mn =  1e8;
        float mx =  -1e8;
        int2 mn_pos = int2(-1,-1);
        int2 mx_pos = int2(-1,-1);
        if (mnp.x < -w || mnp.y < -h)
        {
            mn = base_val;
            mx = base_val;
            mn_pos = to_int2(from);
            mx_pos = to_int2(from);
        }
        else if (mxp.x >= w || mxp.y >= h)
        {
            mn = base_val;
            mx = base_val;
            mn_pos = int2(w+1,h+1);
            mx_pos = int2(w+1,h+1);
        }
        for (int x = MAX(floor(mnp.x), -w);x<MIN(ceil(mxp.x),w);x++)
        {
            for (int y = MAX(floor(mnp.y), -h);y<MIN(ceil(mxp.y),h);y++)
            {
                float val = get(x,y);
                if (val > mx)
                {
                    mx = val;
                    mx_pos = int2(x,y);
                }
                if (val < mn)
                {
                    mn = val;
                    mn_pos = int2(x,y);
                }
            }
        }

        if (min_v)
            *min_v = mn;
        if (max_v)
            *max_v = mx;
        if (min_pos)
        {
            if (mn_pos.x == -1 && mn_pos.y == -1)
                *min_pos = from;
            else if (mn_pos.x == w+1 && mn_pos.y == h+1)
                *min_pos = to;
            else 
                *min_pos = float2(pos.x,pos.y) + cell_size*float2(mn_pos);
        }
        if (max_pos)
        {
            if (mx_pos.x == -1 && mx_pos.y == -1)
                *max_pos = from;
            else if (mx_pos.x == w+1 && mx_pos.y == h+1)
                *max_pos = to;
            else 
                *max_pos = float2(pos.x,pos.y) + cell_size*float2(mx_pos);
        }
    }

    void Field_2d::add(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos.x == field.pos.x && pos.y == field.pos.y && pos.z == field.pos.z &&
                          size.x == field.size.x && size.y == field.size.y && 
                          w == field.w && h == field.h);
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
                    float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) + field.get_bilinear(ps));
                }
            }
        }
    }
    void Field_2d::sub(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos.x == field.pos.x && pos.y == field.pos.y && pos.z == field.pos.z &&
                          size.x == field.size.x && size.y == field.size.y && 
                          w == field.w && h == field.h);
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
                    float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) - field.get_bilinear(ps));
                }
            }
        }
    }
    void Field_2d::mul(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos.x == field.pos.x && pos.y == field.pos.y && pos.z == field.pos.z &&
                          size.x == field.size.x && size.y == field.size.y && 
                          w == field.w && h == field.h);
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
                    float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) * field.get_bilinear(ps));
                }
            }
        }
    }
    void Field_2d::div(Field_2d &field, bool same_size_expected)
    {
        if (!field.data)
            return;
        bool same_size = (pos.x == field.pos.x && pos.y == field.pos.y && pos.z == field.pos.z &&
                          size.x == field.size.x && size.y == field.size.y && 
                          w == field.w && h == field.h);
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
                    float2 ps = float2(pos.x + cell_size*i, pos.z + cell_size*j);
                    set(i,j,get(i,j) / field.get_bilinear(ps));
                }
            }
        }
    }