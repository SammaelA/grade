#pragma once
#include "preprocessing.h"
#include "tinyEngine/engine.h"

namespace upg
{
    Texture get_mask(Texture &image)
    {
        //TODO: image to mask
        logerr("preprocessing: get_mask is not implemented");
        return {}; 
    }

    ReferenceView preprocess_get_reference_view(const Block &view_blk)
    {
        ReferenceView rv;

        Block *synt_reference = view_blk.get_block("synthetic reference");
        if (synt_reference)
        {
            //our reference is a set of procedural generator's parameters
            //we should create model from these parameters and render it
            //with given camera (it contains in view_blk too)
            logerr("preprocessing: loading from synthetic reference is not implemented");
        }

        std::string mask_dir = view_blk.get_string("mask", "");
        std::string image_dir = view_blk.get_string("image", "");

        if (mask_dir == "")
        {
            if (image_dir == "")
                logerr("preprocessing: each view block must have mask or image path");
            Texture image = engine::textureManager->load_unnamed_tex(image_dir, 1);
            rv.mask = get_mask(image);
        }
        else
        {
            rv.mask = engine::textureManager->load_unnamed_tex(mask_dir, 1);           
        }
    }
}