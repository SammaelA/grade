#include "metric.h"
#include "common_utils/utility.h"
//#define STB_IMAGE_IMPLEMENTATION
#include "third_party/stb_image.h"
#include "core/tree.h"
#include "tinyEngine/postfx.h"
unsigned char reference_raw[4*256*256];
double CompressionMetric::get(GrovePacked &g)
{
    int origins = 0;
    int instances = g.instancedBranches.size();
    for (int j = 0; j < g.instancedCatalogue.levels(); j++)
    {
        origins += g.instancedCatalogue.get_level_size(j);
    }
    float comp = (float)instances / origins;//bullshit
    logerr("instanced branches %d/%d = %f", instances, origins, comp);
    return comp;
}
ImpostorMetric::ImpostorMetric(TreeSilhouette tree_sil):
reference(textureManager.empty())
{
    w = (int)Quality::ULTRALOW;
    h = (int)Quality::ULTRALOW;
    //reference_raw = safe_new<unsigned char>(4*w*h, "metric_reference_raw");
    int type = 3;
    for (int i=0;i<w;i++)
    {
        for (int j = 0;j<h;j++)
        {
            type = 3;
            uint index = 4*(i*h + j);
            float x = (float)j/w;
            float y = 1 - (float)i/h;
            if (y <= tree_sil.trunk_height)
            {
                float d = y/tree_sil.trunk_height;
                float trunk_r = d*tree_sil.trunk_up_r + (1-d)*tree_sil.trunk_down_r;
                if (abs(0.5 - x) < trunk_r)
                    type = 1;
            }
            float ellips = pow(abs(x-0.5)/tree_sil.crown_width_r, tree_sil.crown_ellipsoid_power) + 
            pow(abs(y-tree_sil.crown_center_height)/tree_sil.crown_height_r, tree_sil.crown_ellipsoid_power); 
            if (ellips < 1)
                type = 2;

            if (type == 1)
            {
                reference_raw[index] = 255;
                reference_raw[index+1] = 0;
                reference_raw[index+2] = 0;
                reference_raw[index+3] = 255;
            }
            else if (type == 2)
            {
                reference_raw[index] = 0;
                reference_raw[index+1] = 255;
                reference_raw[index+2] = 0;
                reference_raw[index+3] = 255;
            }
            else if (type == 3)
            {
                reference_raw[index] = 0;
                reference_raw[index+1] = 0;
                reference_raw[index+2] = 0;
                reference_raw[index+3] = 0;
            }
            //if (i % 8 == 0 && j % 8 == 0)
            //    debug("%d ",(int)(reference_raw[index+3] > 0));
        }
        //if (i % 8 == 0)
        //    debugnl();
    }
}
ImpostorMetric::ImpostorMetric(Texture &_reference_image):
reference(_reference_image)
{
    if (reference.type != GL_TEXTURE_2D || !reference.is_valid())
    {
        logerr("only a valid 2d RGBA texture could be a reference for impostor metric");
        //reference_raw = nullptr;
        w = 0;
        h = 0;
    }
    else
    {
        w = (int)Quality::ULTRALOW;
        h = (int)Quality::ULTRALOW;
        GLuint fbo;
        glGenFramebuffers(1, &fbo);
        Texture ref_resized = textureManager.create_unnamed(w,h);
        PostFx copy = PostFx("copy.fs");
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ref_resized.texture, 0);
        copy.use();
        copy.get_shader().texture("tex",reference.texture);
        copy.render();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        //reference_raw = safe_new<unsigned char>(4*w*h, "metric_reference_raw");

        glBindTexture(GL_TEXTURE_2D, ref_resized.texture);

        glGetTexImage(GL_TEXTURE_2D,
                    0,
                    GL_RGBA,
                    GL_UNSIGNED_BYTE,
                    reference_raw);
        glBindTexture(GL_TEXTURE_2D, 0);
        //textureManager.save_bmp(ref_resized,"imp");
        glDeleteTextures(1, &ref_resized.texture);
    }
    
    //int channels = 0;
    //reference_raw = stbi_load("resources/textures/reference_tree_test.png",&w,&h,&channels,0);
    //logerr("loaded channels = %d",channels);
}

ImpostorMetric::~ImpostorMetric()
{
    //stbi_image_free(reference_raw);
    //safe_delete<unsigned char>(reference_raw, "metric_reference_raw");
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
        //textureManager.save_bmp(imp,"imp");
        
        for (int i=0;i<w;i+=8)
        {
            for (int j = 0;j<h;j+=8)
            {
                uint imp_index = 4*( i*(h) + j);
                int t1 = get_type(reference_raw[imp_index], reference_raw[imp_index+1]);
                //if (t1 != 3)
                //    debug("%d|",reference_raw[imp_index + 3]);
                //else 
                //    debug(" |");
            }
            //debugnl();
        }
        //textureManager.save_bmp_raw(reference_raw,w,h,4,"imp");
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
        textureManager.save_bmp_raw(imp_raw,iw,ih,4,"imp");
        //textureManager.save_bmp(imp,"imp");
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
    int t1_pixels = 0, t2_pixels = 0, filled_pixels = 0;
    int t1_hits = 0, t2_hits = 0, filled_hits = 0;
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
            
            t1_pixels += (t1 == 1) + (t2 == 1);
            t2_pixels += (t1 == 2) + (t2 == 2);
            filled_pixels += (t1 != 3) + (t2 != 3);

            t1_hits += (t1 == 1) && (t2 == 1);
            t2_hits += (t1 == 2) && (t2 == 2);
            filled_hits += (t1 != 3) && (t2 != 3);

            df += (t1 != t2);
            if (t1 != 3 && t2 != 3)
            {
                imp[imp_index] = 255;
                imp[imp_index+1] = 255;
                imp[imp_index+2] = 255;
                imp[imp_index+3] = 255;
            }
            else if (t1 != 3 && t2 == 3)
            {
                imp[imp_index] = 0;
                imp[imp_index+1] = 255;
                imp[imp_index+2] = 0;
                imp[imp_index+3] = 255;
            }
            else if (t1 == 3 && t2 != 3)
            {
                imp[imp_index] = 0;
                imp[imp_index+1] = 0;
                imp[imp_index+2] = 255;
                imp[imp_index+3] = 255;
            }
            else if (t1 == 3 && t2 == 3)
            {
                imp[imp_index] = 0;
                imp[imp_index+1] = 0;
                imp[imp_index+2] = 0;
                imp[imp_index+3] = 0;
            }
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
    float d_t1 = 2.0f*t1_hits/t1_pixels;
    float d_t2 = 2.0f*t2_hits/t2_pixels;
    float d_fill = 2.0f*filled_hits/filled_pixels;

    float q_t1 = 0;
    float q_t2 = 0;
    float q_fill = 1;
    float norm_q = 3*(q_t1 + q_t2 + q_fill);
    float new_diff = 1 - (q_t1*d_t1 + q_t2*d_t2+q_fill*d_fill)/norm_q;
    //logerr("difference calculated %f %f %f = %f", d_t1, d_t2, d_fill, new_diff);
    //logerr("difference calculated %f", (float)df/sz);
    return new_diff;
    return (float)df/sz;
}

double DummyMetric::get(GrovePacked &g)
{
    return 1;
}