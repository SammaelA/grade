#include "helpers.h"

using namespace glm;

glm::vec3 canonical_bbox()
{
    return glm::vec3(100,50,50);
}
bool get_dedicated_bbox(Branch *branch, BBox &bbox)
{
    if (!branch || branch->segments.empty())
        return false;
    vec3 a(0, 0, 0);
    vec3 b(0, 0, 0);
    vec3 c;
    a = normalize(branch->joints.back().pos - branch->joints.front().pos);
    for (Joint &j : branch->joints)
    {
        for (Branch *br : j.childBranches)
        {
            b += br->joints.back().pos - br->joints.front().pos;
        }
    }
    if (length(cross(a, b)) < 0.01)
        b = vec3(0, 1, 0);
    if (length(cross(a, b)) < 0.01)
        b = vec3(0, 0, 1);
    b = normalize(b - dot(a, b) * a);
    c = cross(a, b);

    bbox = BillboardCloudRaw::get_bbox(branch, a, b, c);
    return true;
}
void voxelize_original_branch(Branch *b, LightVoxelsCube *light, int level_to, float scale)
{
    glm::vec3 second_vec = glm::vec3(b->plane_coef.x,b->plane_coef.y,b->plane_coef.z);
    for (Segment &s : b->segments)
    {
        glm::vec3 dir = s.end - s.begin;
        glm::vec3 n = glm::normalize(glm::cross(dir, second_vec));
        float len = length(dir);
        float v = (1/3.0)*PI*(SQR(s.rel_r_begin) + s.rel_r_begin*s.rel_r_end + SQR(s.rel_r_end));
        float R = 0.5*scale*(s.rel_r_begin + s.rel_r_end);
        int samples = MIN(5 + 0.25*v, 50);

        //approximate partial cone with cylinder
        //uniformly distributed points in cylinder
        for (int i = 0; i<samples; i++)
        {
            float phi = urand(0, 2*PI);
            float h = urand(0, 1);
            float r = R*sqrtf(urand(0, 1));
            glm::vec3 nr = glm::rotate(glm::mat4(1),phi,dir)*glm::vec4(n,0);
            glm::vec3 sample = s.begin + h*dir + r*nr;
            light->set_occluder_trilinear(sample, 25/samples);
        } 
    }
    if (b->level < level_to)
    {
        for (Joint &j : b->joints)
        {
            for (Branch *br : j.childBranches)
            {
                voxelize_original_branch(br, light, level_to, scale);
            }
        }
    }
}

void set_occlusion(Branch *b, LightVoxelsCube *light)
{
    for (Joint &j : b->joints)
    {
        if (j.leaf)
            light->set_occluder_trilinear(j.pos, 1);
        for (Branch *br : j.childBranches)
        {
            set_occlusion(br, light);
        }
    }
}