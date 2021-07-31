#include "metric.h"
#include "tinyEngine/utility.h"

double CompressionMetric::get(GrovePacked &g)
{
    int origins = 0;
    int instances = g.instancedBranches.size();
    for (int j = 0; j < g.instancedCatalogue.levels(); j++)
    {
        origins += g.instancedCatalogue.get_level(j).size();
    }
    float comp = (float)instances / origins;
    logerr("instanced branches %d/%d = %f", instances, origins, comp);
    return comp;
}

ImpostorMetric::ImpostorMetric(Texture &_reference_image):
reference(_reference_image)
{
    if (reference.type != GL_TEXTURE_2D || !reference.is_valid())
    {
        logerr("only a valid 2d RGBA texture could be a reference for impostor metric");
        reference_raw = nullptr;
        w = 0;
        h = 0;
    }
    else
    {
        w = reference.get_W();
        h = reference.get_H();
        reference_raw = safe_new<unsigned char>(4*w*h, "metric_reference_raw");

        glBindTexture(GL_TEXTURE_2D, reference.texture);

        glGetTexImage(GL_TEXTURE_2D,
                    0,
                    GL_RGBA,
                    GL_UNSIGNED_BYTE,
                    reference_raw);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
}

ImpostorMetric::~ImpostorMetric()
{
    safe_delete<unsigned char>(reference_raw, "metric_reference_raw");
}

double ImpostorMetric::get(GrovePacked &g)
{
    if (g.impostors.size() < 1 || !g.impostors[1].valid)
        return 0;
    Texture &imp = g.impostors[1].atlas.tex(0);
    int iw, ih;
    unsigned char *imp_raw = nullptr;
    if (imp.type != GL_TEXTURE_2D_ARRAY || !imp.is_valid())
    {
        logerr("only a valid 2d RGBA texture array with impostors could be used for impostor metric");
        imp_raw = nullptr;
        iw = 0;
        ih = 0;
        return 0;
    }
    else
    {
        iw = imp.get_W();
        ih = imp.get_H();

        imp_raw = safe_new<unsigned char>(4*iw*ih, "metric_impostor_raw");
        glBindTexture(GL_TEXTURE_2D_ARRAY, imp.texture);

        glGetTexImage(GL_TEXTURE_2D_ARRAY,
                    0,
                    GL_RGBA,
                    GL_UNSIGNED_BYTE,
                    imp_raw);
        glBindTexture(GL_TEXTURE_2D_ARRAY, 0);
        
        float av_m = 0;
        textureManager.save_bmp(imp,"imp");
        
        /*for (int i=0;i<iw;i+=19)
        {
            for (int j = 0;j<ih;j+=19)
            {
                uint imp_index = 4*( (iw - 1 - i)*(ih) + (ih - 1 - j) );
                int t1 = get_type(imp_raw[imp_index], imp_raw[imp_index+1]);
                if (t1 != 3)
                    debug("%d|",t1);
                else 
                    debug(" |");
            }
            debugnl();
        }*/
        //tard.clear();
        for (int i=0;i<3;i++)
        {
            for (int j=0;j<3;j++)
            {
                if (i == 2 && j == 2)
                {
                    //top view
                    continue;
                }
                float m = 1 - diff(imp_raw,reference_raw,iw,ih,i*w,j*h);
                av_m += m;
            }
        }
        av_m /= 8;

        safe_delete<unsigned char>(imp_raw, "metric_impostor_raw");
        
        return av_m;
    }
    return 0;
}
int ImpostorMetric::get_type(unsigned char r, unsigned char g)
{
    #define D(r,g,v) (abs((r/256.0) - v.x) + abs((g/256.0) - v.y))
    if (r + g < 10)
        return 3;
    else if (r < 10 || g/r > 2.5)
        return 2;
    else
        return 1;
    float d1 = D(r,g,leaves_color);
    float d2 = D(r,g,wood_color);
    float d3 = D(r,g,empty_color);
    if (d1 < d2)
    {
        if (d1 < d3)
            return 1;
        else 
            return 3;
    }
    else
    {
        if (d2 < d3)
            return 2;
        else 
            return 3;
    }
}
float ImpostorMetric::diff(unsigned char *imp, unsigned char *reference, int imp_tex_w, int imp_tex_h, int imp_offset_w,
                           int imp_offset_h)
{
    if (!reference || !imp || w <= 0 || h <= 0 || imp_tex_w - imp_offset_w < w || imp_tex_h - imp_offset_h < h)
        return 1;
    int sz = w*h;
    int df = 0;
            
    for (int i=0;i<w;i++)
    {
        for (int j = 0;j<h;j++)
        {
            uint index = 4*(i*h + j);
            uint imp_w = (imp_tex_w - 1 - (i + imp_offset_w));
            uint imp_h = (imp_tex_h - 1 - (j + imp_offset_h));
            uint imp_index = 4*((imp_w)*(imp_tex_h) + imp_h);
            int t1 = get_type(reference_raw[index], reference_raw[index+1]);
            int t2 = get_type(imp[imp_index], imp[imp_index+1]);
            df += (t1 != t2);
            /*if (i % 8 == 0 && j % 8 == 0)
            {
                if (t2 != 3)
                    debug("[%d]",t2);
                else 
                    debug("[ ]",t2);
            }*/
        }
        //if (i % 8 == 0)
        //    debugnl();
    }
    logerr("difference calculated %f", (float)df/sz);
    return (float)df/sz;
}

double DummyMetric::get(GrovePacked &g)
{
    return 1;
}