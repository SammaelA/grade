#include "modeling.h"
#include "common_utils/utility.h"
#include "tinyEngine/engine.h"

#define PI 3.14159265f

namespace visualizer
{
void box_to_model(Box *box, Mesh *m)
{
    std::vector<float3> pos;
    std::vector<float2> colors;
    for (int i = 0; i <= 1; i++)
    {
        for (int j = 0; j <= 1; j++)
        {
            for (int k = 0; k <= 1; k++)
            {
                pos.push_back(box->pos + (float)i * box->a + (float)j * box->b + (float)k * box->c);
                colors.push_back(float2(k, j));
            }
        }
    }
    int indicies[36] = {0, 2, 1, 3, 1, 2, 2, 3, 6, 7, 6, 3, 7, 5, 6, 4, 6, 5, 0, 1, 4, 5, 1, 4, 0, 4, 2, 6, 2, 4, 1, 5, 3, 7, 3, 5};
    int _b = m->positions.size() / 3;
    for (int i = 0; i < pos.size(); i++)
    {
        m->positions.push_back(pos[i].x);
        m->positions.push_back(pos[i].y);
        m->positions.push_back(pos[i].z);

        m->colors.push_back(colors[i].x);
        m->colors.push_back(colors[i].y);
        m->colors.push_back(0.0);
        m->colors.push_back(1.0);
    }
    for (int i = 0; i < 36; i++)
    {
        m->indices.push_back(_b + indicies[i]);
    }
}
void ellipsoid_to_model(Ellipsoid *b, Mesh *m, int sectors, int stacks, bool smooth)
{
    float x, y, z, xy;                              // vertex position
    float nx, ny, nz, lengthInv = 1.0f;             // normal
    float s, t;                                     // texCoord

    float sectorStep = 2 * PI / sectors;
    float stackStep = PI / stacks;
    float sectorAngle, stackAngle;
    int _b = m->positions.size() / 3;

    for(int i = 0; i <= stacks; ++i)
    {
        stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = cosf(stackAngle);             // r * cos(u)
        z = sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for(int j = 0; j <= sectors; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            float3 pos = to_float3(LiteMath::inverse4x4(b->transform)*float4(x,y,z,1)); 
            m->positions.push_back(pos.x);
            m->positions.push_back(pos.y);
            m->positions.push_back(pos.z);

            // normalized vertex normal
            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;
            m->normals.push_back(nx);
            m->normals.push_back(ny);
            m->normals.push_back(nz);

            // vertex tex coord between [0, 1]
            s = (float)j / sectors;
            t = (float)i / stacks;
            m->colors.push_back(s);
            m->colors.push_back(t);
            m->colors.push_back(0.0);
            m->colors.push_back(1.0);
        }
    }

    // indices
    //  k1--k1+1
    //  |  / |
    //  | /  |
    //  k2--k2+1
    unsigned int k1, k2;

    for(int i = 0; i < stacks; ++i)
    {
        k1 = i * (sectors + 1);     // beginning of current stack
        k2 = k1 + sectors + 1;      // beginning of next stack

        for(int j = 0; j < sectors; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding 1st and last stacks
            if(i != 0)
            {
                m->indices.push_back(_b + k1);
                m->indices.push_back(_b + k2);
                m->indices.push_back(_b + k1 + 1);
            }

            if(i != (stacks-1))
            {
                m->indices.push_back(_b + k1 + 1);
                m->indices.push_back(_b + k2);
                m->indices.push_back(_b + k2 + 1);
            }
        }
    }
}
void cylinder_to_model(Cylinder *b, Mesh *m, int sectors)
{
    int v0 = m->positions.size()/3;
    for (int i=0;i<sectors;i++)
    {
        float3 vert = b->pos - b->c + cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b;
        m->positions.push_back(vert.x);
        m->positions.push_back(vert.y);
        m->positions.push_back(vert.z);

        float3 n = normalize(cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b);
        m->normals.push_back(n.x);
        m->normals.push_back(n.y);
        m->normals.push_back(n.z);

        float col_x = 1 - abs(2.0f*i/sectors - 1);
        m->colors.push_back(col_x);
        m->colors.push_back(0.0);
        m->colors.push_back(0.0);
        m->colors.push_back(1.0);
    }
    int v1 = m->positions.size()/3;
    for (int i=0;i<sectors;i++)
    {
        float3 vert = b->pos + b->c + cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b;
        m->positions.push_back(vert.x);
        m->positions.push_back(vert.y);
        m->positions.push_back(vert.z);

        float3 n = normalize(cos(2*PI*i/sectors)*b->a + sin(2*PI*i/sectors)*b->b);
        m->normals.push_back(n.x);
        m->normals.push_back(n.y);
        m->normals.push_back(n.z);

        float col_x = 1 - abs(2.0f*i/sectors - 1);
        m->colors.push_back(col_x);
        m->colors.push_back(1.0);
        m->colors.push_back(0.0);
        m->colors.push_back(1.0);
    }
    int v_down = m->positions.size()/3;

    float3 vert = b->pos - b->c;
    m->positions.push_back(vert.x);
    m->positions.push_back(vert.y);
    m->positions.push_back(vert.z);

    float3 n = normalize(-b->c);
    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);
    
    m->colors.push_back(0.0);
    m->colors.push_back(1.0);
    m->colors.push_back(0.0);
    m->colors.push_back(1.0);

    int v_up = m->positions.size()/3;

    vert = b->pos + b->c;
    m->positions.push_back(vert.x);
    m->positions.push_back(vert.y);
    m->positions.push_back(vert.z);

    n = normalize(b->c);
    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);
    
    m->colors.push_back(0.0);
    m->colors.push_back(0.0);
    m->colors.push_back(0.0);
    m->colors.push_back(1.0);

    for (int i=0;i<sectors;i++)
    {
        int i1 = (i+1)%sectors;
        m->indices.push_back(v0 + i);
        m->indices.push_back(v1 + i1);
        m->indices.push_back(v0 + i1);

        m->indices.push_back(v1 + i);
        m->indices.push_back(v1 + i1);
        m->indices.push_back(v0 + i);

        m->indices.push_back(v_down);
        m->indices.push_back(v0 + i1);
        m->indices.push_back(v0 + i);

        m->indices.push_back(v_up);
        m->indices.push_back(v1 + i);
        m->indices.push_back(v1 + i1);
    }
}
void body_to_model(Body *b, Mesh *m, bool fixed_tc, float4 tc)
{
    int last_tc = m->colors.size();
    Box *box = dynamic_cast<Box *>(b);
    if (box)
        box_to_model(box,m);
    Ellipsoid *el = dynamic_cast<Ellipsoid *>(b);
    if (el)
        ellipsoid_to_model(el,m,20,20);
    Cylinder *cyl = dynamic_cast<Cylinder *>(b);
    if (cyl)
        cylinder_to_model(cyl,m,20);
    if (fixed_tc)
    {
        for (int i=last_tc;i<m->colors.size() - 3;i+=4)
        {
            m->colors[i] = tc.x;
            m->colors[i+1] = tc.y;
            m->colors[i+2] = tc.z;
            m->colors[i+3] = tc.w;
        }
    }
}

int HMLOD = 0;
float3 get_n(Heightmap &h, float3 terr_pos, float2 step)
{
    terr_pos.y = 0;
    float3 terr_pos1 = terr_pos + float3(step.x,0,0);
    float3 terr_pos2 = terr_pos + float3(0,0,step.y);
    terr_pos.y = h.get_height(terr_pos);
    terr_pos1.y = h.get_height(terr_pos1);
    terr_pos2.y = h.get_height(terr_pos2);
    float3 n = normalize(cross(terr_pos1 - terr_pos, terr_pos2 - terr_pos));
    if (dot(n, float3(0,1,0)) < 0)
        n = -n;
    return n;
}

void add_vertex(Heightmap &h, float3 &terr_pos, Mesh *m, float2 step)
{
    terr_pos.y = h.get_height(terr_pos);
    float3 n = get_n(h, terr_pos, step);

    m->positions.push_back(terr_pos.x);
    m->positions.push_back(terr_pos.y);
    m->positions.push_back(terr_pos.z);

    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);

    float l = 0.1 * length(h.get_grad(terr_pos));
    float2 tc = float2(2.0f*abs(fract(float2(terr_pos.x, terr_pos.z)/47.1f) - 0.5f));
    tc = float2(terr_pos.x + urand(-HMLOD,HMLOD), terr_pos.z + urand(-HMLOD,HMLOD))/27.1f;
    m->colors.push_back(tc.x);
    m->colors.push_back(tc.y);
    m->colors.push_back(0);
    m->colors.push_back(0);
}

void heightmap_to_model(Heightmap &h, Mesh *m, float2 detailed_size, float2 full_size, 
                                    float precision, int LODs)
{
    int x = MAX(pow(2, round(log2(2*detailed_size.x / precision))), 4);
    int y = MAX(pow(2, round(log2(2*detailed_size.y / precision))), 4);
    x = h.get_grid_size().x;
    y = h.get_grid_size().y;
    //x = 4;
    //y = 4;
    float2 step = 2.0f * detailed_size / float2(x, y);
    float3 center_pos = h.get_pos();
    float2 total_size = detailed_size;
    HMLOD = 0;
    //detailed part of hmap
    for (int i = 0; i <= x; i++)
    {
        for (int j = 0; j <= y; j++)
        {
            int ind = m->positions.size() / 3;
            float3 terr_pos = float3(center_pos.x - detailed_size.x + step.x * i, 0, 
                                           center_pos.z - detailed_size.y + step.y * j);
            add_vertex(h, terr_pos, m, step);

            if (i != x && j != y)
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 1);
                m->indices.push_back(ind + (y+1) + 1);
                m->indices.push_back(ind);
                m->indices.push_back(ind + (y+1) + 1);
                m->indices.push_back(ind + (y+1));
            }
        }
    }
    bool sides_spec = y % 2;
    x = x/2 + 2;//we can safly do it, as x and y are powers of 2
    y = y/2;
    step *= 2.0f;
    HMLOD = 1;
    float3 pos;
    float3 start_pos_top = float3(center_pos.x - detailed_size.x - step.x, 0, center_pos.z + detailed_size.y + step.y);
    float3 start_pos_right = float3(center_pos.x + detailed_size.x + step.x, 0, center_pos.z - detailed_size.y);
    float3 start_pos_bottom = float3(center_pos.x - detailed_size.x - step.x, 0, center_pos.z - detailed_size.y - step.y);
    float3 start_pos_left = float3(center_pos.x - detailed_size.x - step.x, 0, center_pos.z - detailed_size.y);
    float3 end_pos_bottom = float3(center_pos.x + detailed_size.x + step.x, 0, center_pos.z - detailed_size.y);
    float3 end_pos_top = float3(center_pos.x + detailed_size.x + step.x, 0, center_pos.z + detailed_size.y);

    while (total_size.x < full_size.x || total_size.y < full_size.y)
    {
       
        float3 start_pos = start_pos_top;
        for (int i = 0;i<x;i++)
        {
            bool less = (i == 0) || (i == x-1);
            int ind = m->positions.size() / 3;

            pos = start_pos + float3(i*step.x,0,0);
            add_vertex(h, pos, m, step);
            pos = start_pos + float3(i*step.x,0,-step.y);
            add_vertex(h, pos, m, step);
            if (!less)
            {
              pos = start_pos + float3((i+0.5)*step.x,0,-step.y);
              add_vertex(h, pos, m, step);
            }
            if (i == x - 1)
                pos = end_pos_top;
            else
                pos = start_pos + float3((i+1)*step.x,0,-step.y);
            add_vertex(h, pos, m, step);
            if (i == x - 1)
                pos = end_pos_top + float3(0,0,step.y);
            else
                pos = start_pos + float3((i+1)*step.x,0,0);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);

            if (less)
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
            else
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 4);
                
                m->indices.push_back(ind + 4);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
        }

        start_pos = start_pos_bottom;
        for (int i = 0;i<x;i++)
        {
            bool less = (i == 0) || (i == x-1);
            int ind = m->positions.size() / 3;

            pos = start_pos + float3(i*step.x,0,0);
            add_vertex(h, pos, m, step);
            pos = start_pos + float3(i*step.x,0,step.y);
            add_vertex(h, pos, m, step);
            if (!less)
            {
                pos = start_pos + float3((i+0.5)*step.x,0,step.y);
                add_vertex(h, pos, m, step);
            }
            
            if (i == x - 1)
                pos = end_pos_bottom;
            else
                pos = start_pos + float3((i+1)*step.x,0,step.y);
            add_vertex(h, pos, m, step);
            if (i == x - 1)
                pos = end_pos_bottom + float3(0,0,-step.y);
            else
                pos = start_pos + float3((i+1)*step.x,0,0);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            if (less)
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
            else
            {
                m->indices.push_back(ind);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 4);
                
                m->indices.push_back(ind + 4);
                m->indices.push_back(ind + 2);
                m->indices.push_back(ind + 3);
            }
        }

    
        
        start_pos = start_pos_left;
        for (int i=0;i<y - sides_spec;i++)
        {
            int ind = m->positions.size() / 3;

            pos = start_pos + float3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + float3(step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + float3(step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + float3(step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + float3(0,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);
        }
        if (sides_spec)
        {
            int i = y -1;
            int ind = m->positions.size() / 3;

            pos = start_pos + float3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + float3(step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + float3(step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + float3(step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos_top + float3(0,0,-step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos_top + float3(step.x,0,-step.y);
            add_vertex(h, pos, m, step);

            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);

            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 3);
            m->indices.push_back(ind + 5);
        }

          start_pos = start_pos_right;
        for (int i=0;i<y - sides_spec;i++)
        {
            int ind = m->positions.size() / 3;

            pos = start_pos + float3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + float3(-step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + float3(-step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + float3(-step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + float3(0,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);
        }
        if (sides_spec)
        {
            int i = y -1;
            int ind = m->positions.size() / 3;

            pos = start_pos + float3(0,0,i*step.y);
            add_vertex(h, pos, m, step);
            pos = start_pos + float3(-step.x,0,i*step.y);
            add_vertex(h, pos, m, step);

            pos = start_pos + float3(-step.x,0,(i+0.5)*step.y);
            add_vertex(h, pos, m, step);
            
            pos = start_pos + float3(-step.x,0,(i+1)*step.y);
            add_vertex(h, pos, m, step);

            pos = end_pos_top;
            add_vertex(h, pos, m, step);
            
            pos = end_pos_top + float3(-step.x,0,0);
            add_vertex(h, pos, m, step);

            m->indices.push_back(ind);
            m->indices.push_back(ind + 1);
            m->indices.push_back(ind + 2);
            
            m->indices.push_back(ind);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 4);
            
            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 2);
            m->indices.push_back(ind + 3);

            m->indices.push_back(ind + 4);
            m->indices.push_back(ind + 3);
            m->indices.push_back(ind + 5);
        }

        //break;
        start_pos_top += float3(-2.0f*step.x, 0, 2.0f*step.y);
        start_pos_bottom += float3(-2.0f*step.x, 0, -2.0f*step.y);
        end_pos_top += float3(2.0f*step.x, 0, step.y);
        start_pos_right += float3(2.0f*step.x, 0, -step.y);
        start_pos_left += float3(-2.0f*step.x, 0, -step.y);
        end_pos_bottom += float3(2.0f*step.x, 0, -step.y);

        sides_spec = sides_spec || y % 2;

        x = x/2 + 2;
        y = y/2 + 1;
        step *= 2.0f;
        HMLOD*=2;
        total_size += step;
    }
}

void visualize_light_voxels(LightVoxelsCube *voxels, Mesh *m, float3 shift, float3 scale)
{
  visualize_light_voxels(voxels, m, 
                         voxels->get_center() - voxels->get_voxel_size() * float3(voxels->get_vox_sizes()),
                         voxels->get_voxel_size() * (2.0f * float3(voxels->get_vox_sizes()) + float3(1)),
                         2.0f * float3(voxels->get_voxel_size()),
                         0.5f * voxels->get_voxel_size(),
                         0.1,
                         shift,
                         scale);
}
void visualize_light_voxels(LightVoxelsCube *voxels, Mesh *m, float3 pos, float3 size, float3 step,
                            float dot_size, float threshold, float3 shift, float3 scale, int mip)
{
  int count = ((int)(size.x / step.x)) * ((int)(size.y / step.y)) * ((int)(size.z / step.z));
  for (float x = pos.x; x < pos.x + size.x; x += step.x)
  {
    for (float y = pos.y; y < pos.y + size.y; y += step.y)
    {
      for (float z = pos.z; z < pos.z + size.z; z += step.z)
      {
        float occ = 1e8;
        int bord = (x < pos.x + 0.5f*step.x) + (x > pos.x + size.x - step.x) +
                   (y < pos.y + 0.5f*step.y) + (y > pos.y + size.y - step.y) +
                   (z < pos.z + 0.5f*step.z) + (z > pos.z + size.z - step.z);
        occ = voxels->get_occlusion_simple_mip(float3(x, y, z), mip);
        if (bord < 2 && (occ < threshold || occ > 1e8))
          continue;
        float4 tex;
        tex.w = 1;
        tex.z = MIN(1, occ / (10 * threshold));
        tex.y = MIN(1, occ / (100 * threshold));
        tex.x = MIN(1, occ / (1000 * threshold));
        if (bord >= 2)
          tex = float4(1,0,0,1);
        Box b = Box(shift + float3(scale.x * x, scale.y * y, scale.z * z), float3(dot_size, 0, 0), float3(0, dot_size, 0), float3(0, 0, dot_size));
        body_to_model(&b, m, true, tex);
      }
    }
  }
}

void visualize_aabb(AABB &box, Mesh *m, float3 &color)
{
  ::std::vector<AABB> boxes = {box};
  ::std::vector<float3> colors = {color};
  visualize_aabb(boxes, m, colors);
}

void visualize_aabb(::std::vector<AABB> &boxes, Mesh *m, ::std::vector<float3> &colors)
{
  if (colors.size() == boxes.size())
  {
    for (int i = 0; i < colors.size(); i++)
    {
      float3 sz = boxes[i].max_pos - boxes[i].min_pos;
      Box *box_ptr = new Box(boxes[i].min_pos, float3(sz.x, 0, 0), float3(0, sz.y, 0), float3(0, 0, sz.z));
      body_to_model(box_ptr, m, true, to_float4(colors[i], 1));
      delete box_ptr;
    }
  }
  else
  {
    logerr("visualize_aabb colors and boxes arrays should have the same size");
  }
}

void simple_mesh_to_model_332(const std::vector<float> &verts, Mesh *m)
{
  const int FPV = 8;
  if (verts.size() == 0 || verts.size() % FPV*3 != 0)
  {
    logerr("simple_mesh_to_model_332: invalid input data. Should have %d*3*n floats", FPV);
    return;
  }
  int verts_cnt = verts.size()/FPV;
  m->positions.resize(3*verts_cnt);
  m->normals.resize(3*verts_cnt);
  m->colors.resize(4*verts_cnt);
  m->indices.resize(verts_cnt);

  for (int i=0;i<verts_cnt;i++)
  {
    m->positions[3*i] = verts[FPV*i];
    m->positions[3*i+1] = verts[FPV*i+1];
    m->positions[3*i+2] = verts[FPV*i+2];

    m->normals[3*i] = verts[FPV*i+3];
    m->normals[3*i+1] = verts[FPV*i+4];
    m->normals[3*i+2] = verts[FPV*i+5];

    m->colors[4*i] = verts[FPV*i+6];
    m->colors[4*i+1] = verts[FPV*i+7];
    m->colors[4*i+2] = 0;
    m->colors[4*i+3] = 1;

    m->indices[i] = i;
  }
}
}