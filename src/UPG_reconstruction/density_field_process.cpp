#include "density_field_process.h"

glm::vec3 df::closest_point_triangle(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c)
{
    //implementation taken from Embree library
    const glm::vec3 ab = b - a;
    const glm::vec3 ac = c - a;
    const glm::vec3 ap = p - a;

    const float d1 = dot(ab, ap);
    const float d2 = dot(ac, ap);
    if (d1 <= 0.f && d2 <= 0.f) return a; //#1

    const glm::vec3 bp = p - b;
    const float d3 = dot(ab, bp);
    const float d4 = dot(ac, bp);
    if (d3 >= 0.f && d4 <= d3) return b; //#2

    const glm::vec3 cp = p - c;
    const float d5 = dot(ab, cp);
    const float d6 = dot(ac, cp);
    if (d6 >= 0.f && d5 <= d6) return c; //#3

    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
    {
        const float v = d1 / (d1 - d3);
        return a + v * ab; //#4
    }
        
    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
    {
        const float v = d2 / (d2 - d6);
        return a + v * ac; //#5
    }
        
    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
    {
        const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + v * (c - b); //#6
    }

    const float denom = 1.f / (va + vb + vc);
    const float v = vb * denom;
    const float w = vc * denom;
    return a + v * ab + w * ac; //#0
}

std::vector<float>
df::create_density_field(const std::vector<float>& model, const VoxelGrid& grid)
{
    std::srand(std::time(nullptr));

    std::vector<float> density_field(grid.dimension * grid.dimension * grid.dimension, 0);
    
    glm::vec3 gridSize = grid.bounds[1] - grid.bounds[0];
    glm::vec3 voxelSize = gridSize / (float)grid.dimension;

    for (int xi = 0; xi < grid.dimension; xi++)
    {
        for (int yi = 0; yi < grid.dimension; yi++)
        {
            for (int zi = 0; zi < grid.dimension; zi++)
            {
                glm::vec3 voxelCoord = glm::vec3(xi, yi, zi) * voxelSize + grid.bounds[0];

                //  sampling every voxel and testing its points if they are inside a mesh

                size_t inside_points_count = 0;
                
                for (int point_ind = 0; point_ind < grid.samplePointsCount; point_ind++)
                {
                    glm::vec3 samplingPoint = glm::vec3(
                        0 + (voxelSize.x - 0) * rand() / (float)RAND_MAX,
                        0 + (voxelSize.y - 0) * rand() / (float)RAND_MAX,
                        0 + (voxelSize.z - 0) * rand() / (float)RAND_MAX
                    );

                    //  check if point is inside
                    
                    glm::vec3 shortest_distance_vec_to_tr;
                    glm::vec3 norm_to_nearest_tr;
                    
                    float shortest_distance = 10000.f;

                    for (int i = 0; i < model.size(); i += 3*FLOAT_PER_VERTEX)
                    {
                        //  Calculate angle between triangle norm and vector(sample point - triangle center) 
                        
                        glm::vec3 point_on_tr = df::closest_point_triangle(
                            samplingPoint, 
                            glm::vec3(model[i+0*FLOAT_PER_VERTEX], model[i+0*FLOAT_PER_VERTEX+1], model[i+0*FLOAT_PER_VERTEX+2]),
                            glm::vec3(model[i+1*FLOAT_PER_VERTEX], model[i+1*FLOAT_PER_VERTEX+1], model[i+1*FLOAT_PER_VERTEX+2]),
                            glm::vec3(model[i+2*FLOAT_PER_VERTEX], model[i+2*FLOAT_PER_VERTEX+1], model[i+2*FLOAT_PER_VERTEX+2])
                        );

                        float distance = glm::length(voxelCoord + samplingPoint - point_on_tr);

                        if (distance < shortest_distance)
                        {
                            shortest_distance = distance;
                            shortest_distance_vec_to_tr = voxelCoord + samplingPoint - point_on_tr;
                            norm_to_nearest_tr = glm::vec3(
                                model[i + 0 * FLOAT_PER_VERTEX + 3],
                                model[i + 0 * FLOAT_PER_VERTEX + 3 + 1],
                                model[i + 0 * FLOAT_PER_VERTEX + 3 + 2]
                            );
                        }
                    }

                    if (glm::dot(shortest_distance_vec_to_tr, norm_to_nearest_tr) < 0)
                    {
                        inside_points_count++;
                    }
                }

                //  Calc density in each voxel
                density_field[(zi * grid.dimension + yi) * grid.dimension + xi] = inside_points_count / (float)grid.samplePointsCount;
            }
        }
    }

    return density_field;
}

void
df::erase(std::vector<float>& density, const VoxelGrid& grid, const float& trashold)
{
    for (size_t i = 0; i < grid.dimension * grid.dimension * grid.dimension; i++)
    {
        if (density[i] < trashold)
        {
            density[i] = 0.f;
        }
    }
}

std::vector<float> 
df::create_sdf(const std::vector<float>& density, const VoxelGrid& grid)
{
    std::vector<glm::vec3> voxel_border_center_coords;
    glm::vec3 gridSize = grid.bounds[1] - grid.bounds[0];
    glm::vec3 voxelSize = gridSize / (float)grid.dimension;

    std::vector<float> sdf_model(grid.dimension * grid.dimension * grid.dimension, 0);

    //  Collect all voxels' centers that are border of future sdf model
    for (int x = 0; x < grid.dimension; x++)
    {
        for (int y = 0; y < grid.dimension; y++)
        {
            for (int z = 0; z < grid.dimension; z++)
            {
                if (density[(z * grid.dimension + y) * grid.dimension + x] > 0 && 
                    density[(z * grid.dimension + y) * grid.dimension + x] < 1)
                {
                    glm::vec3 voxelCenterCoord = glm::vec3(x, y, z) * voxelSize + grid.bounds[0] + voxelSize / 2.f;
                    voxel_border_center_coords.push_back(voxelCenterCoord);
                }
            }   
        }
    }

    //  Calculate distance between current voxel's center and the nearest border voxel's center from voxel_border_center_coords
    for (int x = 0; x < grid.dimension; x++)
    {
        for (int y = 0; y < grid.dimension; y++)
        {
            for (int z = 0; z < grid.dimension; z++)
            {
                if (density[(z * grid.dimension + y) * grid.dimension + x] > 0 && density[(z * grid.dimension + y) * grid.dimension + x] < 1)
                {
                    sdf_model[(z * grid.dimension + y) * grid.dimension + x] = 0;
                }
                else 
                {
                    float shortest_distance = 10000.f, distance = 0;
                    glm::vec3 voxelCenterCoord = glm::vec3(x, y, z) * voxelSize + grid.bounds[0] + voxelSize / 2.f;

                    for (const glm::vec3 &borderVoxelCenterCoord: voxel_border_center_coords)
                    {
                        distance = glm::distance(voxelCenterCoord, borderVoxelCenterCoord);
                        shortest_distance = (distance < shortest_distance) ? distance : shortest_distance;
                    }

                    sdf_model[(z * grid.dimension + y) * grid.dimension + x] = 
                        (density[(z * grid.dimension + y) * grid.dimension + x] > 0) ? -1 * shortest_distance : shortest_distance;
                }
            }   
        }
    }

    return sdf_model;
}

std::vector<float>
df::pipeline(const std::vector<float>& model)
{
    std::vector<float> sdf_model;
    //  First step
    glm::vec3 bounds[2] = {glm::vec3(-1.2f,-1.2f,-1.2f), glm::vec3(1.2f, 1.2f, 1.2f)};
    VoxelGrid grid(64, 5, bounds);

    std::vector<float> density = df::create_density_field(model, grid);

    //  Second
    df::erase(density, grid, 0.5);

    //  Third
    sdf_model = df::create_sdf(density, grid);

    return sdf_model;
}