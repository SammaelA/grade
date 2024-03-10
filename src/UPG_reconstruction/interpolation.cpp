#include "interpolation.h"
#include "density_field_process.h"


float 
interpolation::bilinear(const float &tx, const float &ty, const float &c00, const float &c10, const float &c01, const float &c11)
{
    float a = c00 * (1.f - tx) + c10 * tx;
    float b = c01 * (1.f - tx) + c11 * tx;
    
    return a * (1.f - ty) + b * ty;
}


//?     Guess that grid is [0, 1]x[0, 1]x[0, 1] cube
float 
interpolation::trilinear(const std::vector<float>& sdf_model, glm::vec3& p, size_t dim_size)
{
    float gx = 0, gy = 0, gz = 0, tx = 0, ty = 0, tz = 0; 
    unsigned gxi = 0, gyi = 0, gzi = 0; 

    gx = p.x * dim_size;
    gxi = (int)gx;
    tx = gx - gxi;

    gy = p.y * dim_size;
    gyi = (int)gy;
    ty = gy - gyi;

    gz = p.z * dim_size;
    gzi = (int)gz;
    tz = gz - gzi;

    const float& c000 = sdf_model[df::get_index(dim_size, gxi, gyi, gzi)]; 
    const float& c100 = sdf_model[df::get_index(dim_size,gxi + 1, gyi, gzi)]; 
    const float& c010 = sdf_model[df::get_index(dim_size,gxi, gyi + 1, gzi)]; 
    const float& c110 = sdf_model[df::get_index(dim_size,gxi + 1, gyi + 1, gzi)]; 
    const float& c001 = sdf_model[df::get_index(dim_size,gxi, gyi, gzi + 1)]; 
    const float& c101 = sdf_model[df::get_index(dim_size,gxi + 1, gyi, gzi + 1)]; 
    const float& c011 = sdf_model[df::get_index(dim_size,gxi, gyi + 1, gzi + 1)]; 
    const float& c111 = sdf_model[df::get_index(dim_size,gxi + 1, gyi + 1, gzi + 1)]; 

    float e = interpolation::bilinear(tx, ty, c000, c100, c010, c110);
    float f = interpolation::bilinear(tx, ty, c001, c101, c011, c111);

    return e * (1 - tz) + f * tz;
}