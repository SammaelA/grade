#include "field_2d.h"
#include "tinyEngine/utility.h"
#include <vector>
/* Function to linearly interpolate between a0 and a1
 * Weight w should be in the range [0.0, 1.0]
 */
float interpolate(float a0, float a1, float w) {
    /* // You may want clamping by inserting:
     * if (0.0 > w) return a0;
     * if (1.0 < w) return a1;
     */
    return (a1 - a0) * w + a0;
    /* // Use this cubic interpolation [[Smoothstep]] instead, for a smooth appearance:
     * return (a1 - a0) * (3.0 - w * 2.0) * w * w + a0;
     *
     * // Use [[Smootherstep]] for an even smoother result with a second derivative equal to zero on boundaries:
     * return (a1 - a0) * ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) + a0;
     */
}

typedef struct {
    float x, y;
} vector2;

/* Create random direction vector
 */
vector2 randomGradient(int ix, int iy) {
    // Random float. No precomputed gradients mean this works for any number of grid coordinates
    float random = 2920.f * sin(ix * 21942.f + iy * 171324.f + 8912.f) * cos(ix * 23157.f * iy * 217832.f + 9758.f);
    return (vector2) { .x = cosf(random), .y = sinf(random) };
}

// Computes the dot product of the distance and gradient vectors.
float dotGridGradient(int ix, int iy, float x, float y) {
    // Get gradient from integer coordinates
    vector2 gradient = randomGradient(ix, iy);

    // Compute the distance vector
    float dx = x - (float)ix;
    float dy = y - (float)iy;

    // Compute the dot-product
    return (dx*gradient.x + dy*gradient.y);
}

// Compute Perlin noise at coordinates x, y
float perlin(float x, float y) {
    // Determine grid cell coordinates
    int x0 = (int)x;
    int x1 = x0 + 1;
    int y0 = (int)y;
    int y1 = y0 + 1;

    // Determine interpolation weights
    // Could also use higher order polynomial/s-curve here
    float sx = x - (float)x0;
    float sy = y - (float)y0;

    // Interpolate between grid point gradients
    float n0, n1, ix0, ix1, value;

    n0 = dotGridGradient(x0, y0, x, y);
    n1 = dotGridGradient(x1, y0, x, y);
    ix0 = interpolate(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);
    ix1 = interpolate(n0, n1, sx);

    value = interpolate(ix0, ix1, sy);
    return value;
}
float f_perlin(float x, float y) 
{
    std::vector<float> x_mul = {1,3,7,11,21};
    std::vector<float> y_mul = {1,5,9,13,23};
    float res = 0.0;
    int it = 0;
    for (float mx : x_mul)
    {
        for (float my : y_mul)
        {
            res += perlin(mx*x,my*y);
            it++;
        }
    }
    return res;
}
    Field_2d::Field_2d(glm::vec3 pos, glm::vec2 size, float cell_size):
    Field_2d(pos, ceil(size.x/cell_size), ceil(size.y/cell_size))
    {
        this->cell_size = cell_size;
    }
    Field_2d::Field_2d(glm::vec3 _pos, int _w, int _h)
    {
        w = _w;
        h = _h;
        pos = _pos;

        data = safe_new<float>((2*w + 1)*(2*h + 1),"field_2d_data");
    }
    Field_2d::~Field_2d()
    {
        safe_delete<float>(data,"field_2d_data");
    }
    float Field_2d::get(int x, int y)
    {
        if (x >= -w && x <= w && y >= -h && y <= h)
        {
            return data[(2*w + 1)*(y + h) + x + w];
        }
        else
            return base_val;
    }
    float Field_2d::get_bilinear(glm::vec3 position)
    {
        if (!data)
            return base_val;
        glm::vec2 rp = glm::vec2(position.x - pos.x, position.z - pos.z)/cell_size;
        glm::ivec2 ps = rp;
        float dx = rp.x - ps.x;
        float dy = rp.y - ps.y;
        return (dx*get(ps.x, ps.y) + (1 - dx)*get(ps.x + 1, ps.y))*(1 - dy) + 
               (dx*get(ps.x, ps.y + 1) + (1 - dx)*get(ps.x + 1, ps.y + 1))*dy;
        
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
    glm::vec2 Field_2d::get_grad(int x, int y)
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
    glm::vec2 Field_2d::get_grad_bilinear(glm::vec3 position)
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
    glm::vec4 Field_2d::get_borders()
    {
        return glm::vec4(pos.x - size.x, pos.z - size.y, pos.x + size.x, pos.z + size.y);
    }
    void Field_2d::print()
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
    void Field_2d::get_min_max_imprecise(glm::vec2 from, glm::vec2 to, float *min_v, float *max_v, 
                                         glm::vec2 *min_pos, glm::vec2 *max_pos)
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